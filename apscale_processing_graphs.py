#!/usr/bin/env python3

"""
A script to generate processing graphs for apscale runs. Must be started from the
parent directory of an apscale project directory.

Requirements: apscale project directory that contains a Project_report.xlsx file
with sheets named 3_PE merging, 4_primer_trimming, and 5_quality_filtering,
as well as the folders 7_otu_clustering, 8_denoising, and 9_lulu_filtering
with all result files.

By Chris Hempel (christopher.hempel@kaust.edu.sa) on Jun 7 2023
"""

import os
import datetime
import pandas as pd
import plotly.express as px
import argparse
import warnings

# Define that warnings are not printed to console
warnings.filterwarnings("ignore")


# Funtion to print datetime and text
def time_print(text):
    datetime_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{datetime_now}  ---  " + text)


# Define arguments
parser = argparse.ArgumentParser(
    description="""A script to generate processing graphs for apscale runs.""",
)
parser.add_argument(
    "-p",
    "--project_dir",
    help="Directory containing apscale results to generate reports for.",
    required=True,
)
parser.add_argument(
    "-f",
    "--graph_format",
    help="Graph format, either png or svg.",
    required=True,
    choices=["png", "svg"],
)
parser.add_argument(
    "-S",
    "--scaling_factor",
    help="Scaling factor for graph width. Manual trial and error in 0.2 increments might be required (default: 1).",
    default=1,
    type=float,
)
args = parser.parse_args()

# Set arguments
project_dir = args.project_dir
graph_format = args.graph_format
scaling_factor = args.scaling_factor
project_name = os.path.basename(project_dir)

# Import files
time_print("Importing files...")
report_file = os.path.join(project_dir, "Project_report.xlsx")
esv_postlulu_file = os.path.join(
    project_dir,
    "9_lulu_filtering",
    "denoising",
    f"{project_dir}_ESV_table_filtered.parquet.snappy",
)
esv_prelulu_file = os.path.join(
    project_dir, "8_denoising", f"{project_name}_ESV_table.parquet.snappy"
)
otu_postlulu_file = os.path.join(
    project_dir,
    "9_lulu_filtering",
    "otu_clustering",
    f"{project_name}_OTU_table_filtered.parquet.snappy",
)
otu_prelulu_file = os.path.join(
    project_dir,
    "7_otu_clustering",
    f"{project_name}_OTU_table.parquet.snappy",
)
report_sheet_dict = pd.read_excel(report_file, sheet_name=None)
esv_postlulu_df = pd.read_parquet(esv_postlulu_file, engine="fastparquet")
time_print("1/4 files imported...")
esv_prelulu_df = pd.read_parquet(esv_prelulu_file, engine="fastparquet")
time_print("2/4 files imported...")
otu_postlulu_df = pd.read_parquet(otu_postlulu_file, engine="fastparquet")
time_print("3/4 files imported...")
otu_prelulu_df = pd.read_parquet(otu_prelulu_file, engine="fastparquet")
time_print("Import done. Generating graphs...")

# ESV table processing
esv_postlulu_df = esv_postlulu_df.drop("Seq", axis=1).set_index("ID")
esv_postlulu_df[esv_postlulu_df > 1] = 1
esv_postlulu_sums = esv_postlulu_df.sum()
esv_postlulu_sums.name = "ESVs postlulu"
esv_postlulu_sums.index = (
    esv_postlulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)
esv_prelulu_df = esv_prelulu_df.drop("Seq", axis=1).set_index("ID")
esv_prelulu_df[esv_prelulu_df > 1] = 1
esv_prelulu_sums = esv_prelulu_df.sum()
esv_prelulu_sums.name = "ESVs prelulu"
esv_prelulu_sums.index = (
    esv_prelulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)
# OTU table processing
otu_postlulu_df = otu_postlulu_df.drop("Seq", axis=1).set_index("ID")
otu_postlulu_df[otu_postlulu_df > 1] = 1
otu_postlulu_sums = otu_postlulu_df.sum()
otu_postlulu_sums.name = "OTUs postlulu"
otu_postlulu_sums.index = (
    otu_postlulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)
otu_prelulu_df = otu_prelulu_df.drop("Seq", axis=1).set_index("ID")
otu_prelulu_df[otu_prelulu_df > 1] = 1
otu_prelulu_sums = otu_prelulu_df.sum()
otu_prelulu_sums.name = "OTUs prelulu"
otu_prelulu_sums.index = (
    otu_prelulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)

# Graphs
# Set graph width based on number of samples (note: doesn't work consistently)
graph_width = 400 + len(report_sheet_dict["3_PE merging"]) * 25 * scaling_factor

# PE merging
perc_kept_pe = pd.Series(
    report_sheet_dict["3_PE merging"]["merged reads"]
    / report_sheet_dict["3_PE merging"]["processed reads"]
    * 100
)
perc_kept_pe.index = (
    report_sheet_dict["3_PE merging"]["File"]
    .str.replace("_PE", "")
    .str.replace(".fastq.gz", "")
)  # type: ignore
perc_kept_pe = perc_kept_pe.fillna(0)
pe_graph = px.bar(
    perc_kept_pe.sort_values(),
    range_y=[0, 100],
    labels={"value": "Reads kept [%]", "File": "Samples"},
    title="Percentage of reads kept during PE merging",
    width=graph_width,
)
pe_graph.add_hline(
    y=(sum(perc_kept_pe) / len(perc_kept_pe)),
    line_width=1,
    line_dash="dash",
)
pe_graph.update_layout(showlegend=False)
pe_graph.update_xaxes(tickangle=55)

# Primer trimming
perc_kept_trim = pd.Series(
    report_sheet_dict["4_primer_trimming"]["trimmed reads"]
    / report_sheet_dict["4_primer_trimming"]["processed reads"]
    * 100
)
perc_kept_trim.index = (
    report_sheet_dict["4_primer_trimming"]["File"]
    .str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace(".fastq.gz", "")
)  # type: ignore
perc_kept_trim = perc_kept_trim.fillna(0)
trim_graph = px.bar(
    perc_kept_trim.sort_values(),
    range_y=[0, 100],
    labels={"value": "Reads kept [%]", "File": "Samples"},
    title="Percentage of reads kept during trimming",
    width=graph_width,
)
trim_graph.add_hline(
    y=(sum(perc_kept_trim) / len(perc_kept_trim)),
    line_width=1,
    line_dash="dash",
)
trim_graph.update_layout(showlegend=False)
trim_graph.update_xaxes(tickangle=55)

# Quality filtering
perc_kept_qf = (
    report_sheet_dict["5_quality_filtering"]["passed reads"]
    / report_sheet_dict["5_quality_filtering"]["processed reads"]
    * 100
)
perc_kept_qf.index = (
    report_sheet_dict["5_quality_filtering"]["File"]
    .str.replace("_PE", "")
    .str.replace("_trimmed_filtered.fasta.gz", "")
    .str.replace("_trimmed.fasta.gz", "")
    .str.replace(".fasta.gz", "")
    .str.replace("_filtered.fasta.gz", "")
)  # type: ignore
perc_kept_qf = perc_kept_qf.fillna(0)
qf_graph = px.bar(
    perc_kept_qf.sort_values(),
    range_y=[0, 100],
    labels={"value": "Reads kept [%]", "File": "Samples"},
    title="Percentage of reads kept during quality filtering",
    width=graph_width,
)
qf_graph.add_hline(
    y=(sum(perc_kept_qf) / len(perc_kept_qf)),
    line_width=1,
    line_dash="dash",
)
qf_graph.update_layout(showlegend=False)
qf_graph.update_xaxes(tickangle=55)

# Raw number of reads
num_reads_raw = report_sheet_dict["3_PE merging"]["processed reads"]
num_reads_raw.index = (
    report_sheet_dict["3_PE merging"]["File"]
    .str.replace("_PE.fastq.gz", "")
    .str.replace(".fastq.gz", "")
)  # type: ignore
ymax_reads = max(num_reads_raw)
rawreads_graph = px.bar(
    num_reads_raw.sort_values(),
    labels={"value": "Number of reads", "File": "Samples"},
    title="Raw number of reads",
    range_y=[0, ymax_reads],
    width=graph_width,
)
rawreads_graph.add_hline(
    y=(sum(num_reads_raw) / len(num_reads_raw)), line_width=1, line_dash="dash"
)
rawreads_graph.update_layout(showlegend=False)
rawreads_graph.update_xaxes(tickangle=55)

# Number of reads after PE merging and quality filtering
num_reads_filtered = report_sheet_dict["5_quality_filtering"]["passed reads"]
num_reads_filtered.index = (
    report_sheet_dict["5_quality_filtering"]["File"]
    .str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace(".fasta.gz", "")
)  # type: ignore
filteredreads_graph = px.bar(
    num_reads_filtered.sort_values(),
    labels={"value": "Number of reads", "File": "Samples"},
    title="Number of reads after PE merging, trimming and quality filtering",
    range_y=[0, ymax_reads],
    width=graph_width,
)
filteredreads_graph.add_hline(
    y=(sum(num_reads_filtered) / len(num_reads_filtered)),
    line_width=1,
    line_dash="dash",
)
filteredreads_graph.update_layout(showlegend=False)
filteredreads_graph.update_xaxes(tickangle=55)

# Number of ESVs per sample before LULU
ymax_esvs = esv_prelulu_sums.max()
prelulu_graph_esvs = px.bar(
    esv_prelulu_sums.sort_values(),
    labels={"value": "Number of ESVs", "index": "Samples"},
    title="Number of ESVs per sample before LULU",
    range_y=[0, ymax_esvs],
    width=graph_width,
)
prelulu_graph_esvs.add_hline(
    y=(sum(esv_prelulu_sums) / len(esv_prelulu_sums)),
    line_width=1,
    line_dash="dash",
)
prelulu_graph_esvs.update_layout(showlegend=False)
prelulu_graph_esvs.update_xaxes(tickangle=55)

# Number of ESVs per sample after LULU
postlulu_graph_esvs = px.bar(
    esv_postlulu_sums.sort_values(),
    labels={"value": "Number of ESVs", "index": "Samples"},
    title="Number of ESVs per sample after LULU",
    range_y=[0, ymax_esvs],
    width=graph_width,
)
postlulu_graph_esvs.add_hline(
    y=(sum(esv_postlulu_sums) / len(esv_postlulu_sums)),
    line_width=1,
    line_dash="dash",
)
postlulu_graph_esvs.update_layout(showlegend=False)
postlulu_graph_esvs.update_xaxes(tickangle=55)

# Number of OTUs per sample before LULU
ymax_otus = otu_prelulu_sums.max()
prelulu_graph_otus = px.bar(
    otu_prelulu_sums.sort_values(),
    labels={"value": "Number of OTUs", "index": "Samples"},
    title="Number of OTUs per sample before LULU",
    range_y=[0, ymax_otus],
    width=graph_width,
)
prelulu_graph_otus.add_hline(
    y=(sum(otu_prelulu_sums) / len(otu_prelulu_sums)),
    line_width=1,
    line_dash="dash",
)
prelulu_graph_otus.update_layout(showlegend=False)
prelulu_graph_otus.update_xaxes(tickangle=55)

# Number of OTUs per sample after LULU
postlulu_graph_otus = px.bar(
    otu_postlulu_sums.sort_values(),
    labels={"value": "Number of OTUs", "index": "Samples"},
    title="Number of OTUs per sample after LULU",
    range_y=[0, ymax_otus],
    width=graph_width,
)
postlulu_graph_otus.add_hline(
    y=(sum(otu_postlulu_sums) / len(otu_postlulu_sums)),
    line_width=1,
    line_dash="dash",
)
postlulu_graph_otus.update_layout(showlegend=False)
postlulu_graph_otus.update_xaxes(tickangle=55)

# Number of reads vs number of ESVs
reads_esvs_graph = px.scatter(
    pd.concat([num_reads_filtered, esv_postlulu_sums], axis=1).reset_index(),
    x="passed reads",
    y="ESVs postlulu",
    trendline="ols",
    trendline_color_override="black",
    hover_name="index",
    title=project_name,
    width=480,
    labels={
        "ESVs postlulu": "Number of ESVs after LULU",
        "passed reads": "Number of filtered reads",
    },
)

# Number of reads vs number of OTUs
reads_otus_graph = px.scatter(
    pd.concat([num_reads_filtered, otu_postlulu_sums], axis=1).reset_index(),
    x="passed reads",
    y="OTUs postlulu",
    trendline="ols",
    trendline_color_override="black",
    hover_name="index",
    title=project_name,
    width=480,
    labels={
        "OTUs postlulu": "Number of OTUs after LULU",
        "passed reads": "Number of filtered reads",
    },
)


# Save graphs
## Make outdir for project_dir if it doesn't already exist
outdir = os.path.join(project_dir, "processing_graphs")
os.makedirs(outdir, exist_ok=True)

pe_graph.write_image(
    os.path.join(
        outdir,
        f"{project_name}_1_pe_merging.{graph_format}",
    )
)
trim_graph.write_image(
    os.path.join(
        outdir,
        f"{project_name}_2_trimming.{graph_format}",
    )
)
qf_graph.write_image(
    os.path.join(
        outdir,
        f"{project_name}_3_qualityFiltering.{graph_format}",
    )
)
rawreads_graph.write_image(
    os.path.join(
        outdir,
        f"{project_name}_4_numberRawreads.{graph_format}",
    )
)
filteredreads_graph.write_image(
    os.path.join(
        outdir,
        f"{project_name}_5_numberFilteredreads.{graph_format}",
    )
)
prelulu_graph_esvs.write_image(
    os.path.join(
        outdir,
        f"{project_name}_6_prelulu_esvs.{graph_format}",
    )
)
postlulu_graph_esvs.write_image(
    os.path.join(
        outdir,
        f"{project_name}_7_postlulu_esvs.{graph_format}",
    )
)
prelulu_graph_otus.write_image(
    os.path.join(
        outdir,
        f"{project_name}_8_prelulu_otus.{graph_format}",
    )
)
postlulu_graph_otus.write_image(
    os.path.join(
        outdir,
        f"{project_name}_9_postlulu_otus.{graph_format}",
    )
)
reads_esvs_graph.write_image(
    os.path.join(
        outdir,
        f"{project_name}_10_reads_vs_esvs.{graph_format}",
    )
)
reads_otus_graph.write_image(
    os.path.join(
        outdir,
        f"{project_name}_11_reads_vs_otus.{graph_format}",
    )
)
