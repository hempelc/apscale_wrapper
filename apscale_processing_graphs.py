#!/usr/bin/env python3

"""
A script to generate processing graphs for apscale runs.

Requirements: apscale project directory that contains a Project_report.xlsx file
with sheets named 3_PE merging, 4_primer_trimming, and 5_quality_filtering,
as well as the folders 4_primer_trimming, 7_otu_clustering, 8_denoising, and 9_lulu_filtering with all result files.

Credit: some graphs are modified versions of the graphs generated by the apscale GUI

By Chris Hempel (christopher.hempel@kaust.edu.sa) on Jun 7 2023
"""

import os
import datetime
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import argparse
import warnings
import dash_bio
import subprocess
from statistics import mean, median, stdev


# Define that warnings are not printed to console
warnings.filterwarnings("ignore")


# Funtion to print datetime and text
def time_print(text):
    print(datetime.datetime.now().strftime("%H:%M:%S"), ": ", text, sep="")


# Function for statistics calculation
def calculate_read_stats(lst):
    return [
        min(lst),
        max(lst),
        round(mean(lst), 2),
        round(median(lst), 2),
        round(stdev(lst), 2),
    ]


# Function to format a df for krona
def krona_formatting(df):
    if database_format == "bold":
        ranks = ["phylum", "class", "order", "family", "genus", "species"]
    else:
        ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    # Sum samples
    sample_sums = (
        df.drop(columns=["ID", "Seq", "lowest_taxon", "lowest_rank"] + ranks)
        .sum(axis=1)
        .rename("Sum")
    )
    krona_df = pd.concat([sample_sums, df[ranks]], axis=1)
    # Fix taxonomy formatting
    ## Turn all non-taxa names into NaN
    krona_df = (
        krona_df.replace("Taxonomy unreliable - multiple matching taxa", np.nan)
        .replace(
            "Taxonomy unreliable - percentage similarity threshold for rank not met",
            np.nan,
        )
        .replace(
            "Taxonomy unreliable - bitscore and alignment length threshold not met",
            np.nan,
        )
        .replace("Not available in database", np.nan)
        .replace("Unknown in PR2 database", np.nan)
        .replace("Unknown in BOLD database", np.nan)
        .replace("_", " ", regex=True)
    )
    ## If entire taxonomy is NaN, replace with "Taxonomy unreliable"
    for index, row in krona_df.iterrows():
        if database_format == "bold":
            if pd.isna(row["phylum"]):
                krona_df.loc[index, "phylum":] = "Taxonomy unreliable"
        elif pd.isna(row["domain"]):
            krona_df.loc[index, "domain":] = "Taxonomy unreliable"
    ## Fill NaNs with last tax entry
    krona_df = krona_df.fillna(method="ffill", axis=1)
    # Aggregate taxa
    krona_df_agg = (
        krona_df.groupby(krona_df.columns[1:].tolist())["Sum"].sum().reset_index()
    )
    # Move the Sum column to the front
    column_order = [krona_df_agg.columns.tolist()[-1]] + krona_df_agg.columns.tolist()[
        :-1
    ]
    krona_df_agg = krona_df_agg[column_order]
    return krona_df_agg.loc[krona_df_agg["Sum"] != 0]


# Define arguments
parser = argparse.ArgumentParser(
    description="A script to generate processing graphs for apscale runs.",
)
parser.add_argument(
    "--project_dir",
    help="Directory containing apscale results to generate reports for.",
    required=True,
)
parser.add_argument(
    "--graph_format",
    help="Graph format, either HTML, png, svg. HTML is recommended because it creates interactive plots (default: html).",
    default="html",
    choices=["png", "svg", "html"],
)
parser.add_argument(
    "--min_length",
    metavar="NNN",
    type=int,
    help="Minimum limit of expected amplicon length (will be plotted in length graph).",
    required=True,
)
parser.add_argument(
    "--max_length",
    metavar="NNN",
    type=int,
    help="Maximum limit of expected amplicon length (will be plotted in length graph).",
    required=True,
)
parser.add_argument(
    "--blast",
    help="Has BLAST been run on the OTU and ESV tables? Information required for krona graphs.",
    choices=["True", "False"],
)
parser.add_argument(
    "--database_format",
    help="Format of the database used for BLAST. Currently available formats are: midori2, pr2, silva, bold. Information required for krona graphs.",
    choices=["midori2", "pr2", "silva", "bold"],
)
parser.add_argument(
    "--remove_negative_controls",
    help="Have negatve controls been removed? Information required for krona graphs.",
    choices=["True", "False"],
)
parser.add_argument(
    "--scaling_factor",
    help="Scaling factor for graph width. Manual trial and error in 0.2 increments might be required (default: 1).",
    default=1,
    type=float,
)

args = parser.parse_args()

# Set arguments
project_dir = args.project_dir
graph_format = args.graph_format
min_length = args.min_length
max_length = args.max_length
scaling_factor = args.scaling_factor
blast = args.blast
database_format = args.database_format
project_name = os.path.basename(project_dir)
if args.remove_negative_controls == "True":
    microdecon_suffix = "_microdecon-filtered"
else:
    microdecon_suffix = ""

# Start of pipeline
time_print("Generating apscale processing graphs...")

# Make outdir for project_dir if it doesn't already exist
outdir = os.path.join(project_dir, "0_statistics_and_graphs")
os.makedirs(outdir, exist_ok=True)

# Import files
time_print("Importing files. This can take a while...")
report_file = os.path.join(project_dir, "Project_report.xlsx")
esv_postlulu_file = os.path.join(
    project_dir,
    "9_lulu_filtering",
    "denoising",
    f"{project_name}_ESV_table_filtered.parquet.snappy",
)
esv_prelulu_file = os.path.join(
    project_dir, "8_denoising", f"{project_name}_ESV_table.parquet.snappy"
)
esv_final_file = os.path.join(
    project_dir,
    "9_lulu_filtering",
    "denoising",
    f"{project_name}_ESV_table_filtered{microdecon_suffix}_with_taxonomy.csv",
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
otu_final_file = os.path.join(
    project_dir,
    "9_lulu_filtering",
    "otu_clustering",
    f"{project_name}_OTU_table_filtered{microdecon_suffix}_with_taxonomy.csv",
)

report_sheet_dict = pd.read_excel(report_file, sheet_name=None)
esv_postlulu_df = pd.read_parquet(esv_postlulu_file, engine="fastparquet")
time_print("1/6 files imported...")
esv_prelulu_df = pd.read_parquet(esv_prelulu_file, engine="fastparquet")
time_print("2/6 files imported...")
esv_final_df = pd.read_csv(esv_final_file)
time_print("3/6 files imported...")
otu_postlulu_df = pd.read_parquet(otu_postlulu_file, engine="fastparquet")
time_print("4/6 files imported...")
otu_prelulu_df = pd.read_parquet(otu_prelulu_file, engine="fastparquet")
time_print("5/6 files imported...")
otu_final_df = pd.read_csv(otu_final_file)
time_print("6/6 files imported. Import done. Generating graphs...")

# ESV table processing
esv_postlulu_df_mod = esv_postlulu_df.drop("Seq", axis=1).set_index("ID")
esv_postlulu_df_mod[esv_postlulu_df_mod > 1] = 1
esv_postlulu_sums = esv_postlulu_df_mod.sum()
esv_postlulu_sums.name = "count"
esv_postlulu_sums.index = (
    esv_postlulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)
esv_postlulu_sums = esv_postlulu_sums.reset_index()

esv_prelulu_df = esv_prelulu_df.drop("Seq", axis=1).set_index("ID")
esv_prelulu_df[esv_prelulu_df > 1] = 1
esv_prelulu_sums = esv_prelulu_df.sum()
esv_prelulu_sums.name = "count"
esv_prelulu_sums.index = (
    esv_prelulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)
esv_prelulu_sums = esv_prelulu_sums.reset_index()

# OTU table processing
otu_postlulu_df_mod = otu_postlulu_df.drop("Seq", axis=1).set_index("ID")
otu_postlulu_df_mod[otu_postlulu_df_mod > 1] = 1
otu_postlulu_sums = otu_postlulu_df_mod.sum()
otu_postlulu_sums.name = "count"
otu_postlulu_sums.index = (
    otu_postlulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)
otu_postlulu_sums = otu_postlulu_sums.reset_index()

otu_prelulu_df = otu_prelulu_df.drop("Seq", axis=1).set_index("ID")
otu_prelulu_df[otu_prelulu_df > 1] = 1
otu_prelulu_sums = otu_prelulu_df.sum()
otu_prelulu_sums.name = "count"
otu_prelulu_sums.index = (
    otu_prelulu_sums.index.str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
)
otu_prelulu_sums = otu_prelulu_sums.reset_index()

# Stats table processing
## Gather stats
samples = (
    report_sheet_dict["3_PE merging"]["File"].str.replace("_PE.fastq.gz", "").tolist()
)
raw_reads = report_sheet_dict["3_PE merging"]["processed reads"].values.tolist()
merged_reads = report_sheet_dict["3_PE merging"]["merged reads"].values.tolist()
trimmed_reads = report_sheet_dict["4_primer_trimming"]["trimmed reads"].values.tolist()
filtered_reads = report_sheet_dict["5_quality_filtering"][
    "passed reads"
].values.tolist()
mapped_reads_OTUs = [sum(otu_postlulu_df[sample].values.tolist()) for sample in samples]
mapped_reads_ESVs = [sum(esv_postlulu_df[sample].values.tolist()) for sample in samples]
## Make stats df
df_stats = pd.DataFrame()
df_stats["Raw reads"] = raw_reads + calculate_read_stats(raw_reads)
df_stats["Merged reads"] = merged_reads + calculate_read_stats(merged_reads)
df_stats["Trimmed reads"] = trimmed_reads + calculate_read_stats(trimmed_reads)
df_stats["Filtered reads"] = filtered_reads + calculate_read_stats(filtered_reads)
df_stats["Mapped reads (OTUs)"] = mapped_reads_OTUs + calculate_read_stats(
    mapped_reads_OTUs
)
df_stats["Mapped reads (ESVs)"] = mapped_reads_ESVs + calculate_read_stats(
    mapped_reads_ESVs
)
df_stats.index = samples + [
    "_minimum",
    "_maximum",
    "_average",
    "_median",
    "_deviation",
]
## Save stats df
out_xlsx = os.path.join(
    outdir,
    f"{project_name}_summary_stats.xlsx",
)
df_stats.to_excel(out_xlsx)

# Graphs
# Set graph width based on number of samples (note: doesn't work consistently)
graph_width = 400 + len(report_sheet_dict["3_PE merging"]) * 25 * scaling_factor

# PE merging
perc_kept_pe = pd.DataFrame(
    {
        "value": report_sheet_dict["3_PE merging"]["merged reads"]
        / report_sheet_dict["3_PE merging"]["processed reads"]
        * 100
    }
)
perc_kept_pe["sample"] = (
    report_sheet_dict["3_PE merging"]["File"]
    .str.replace("_PE", "")
    .str.replace(".fastq.gz", "")
)  # type: ignore
perc_kept_pe = perc_kept_pe.fillna(0)
pe_graph = px.bar(
    perc_kept_pe.sort_values("value"),
    x="sample",
    y="value",
    range_y=[0, 100],
    labels={"value": "Reads kept [%]", "sample": "Sample"},
    title="Percentage of reads kept during PE merging",
    width=graph_width,
    # hover_name="sample",
    # hover_data=["value"],
)
pe_graph.add_hline(
    y=(sum(perc_kept_pe["value"]) / len(perc_kept_pe)),
    line_width=1,
    line_dash="dash",
)
pe_graph.update_layout(showlegend=False)
pe_graph.update_xaxes(tickangle=55)
## Save
if graph_format == "html":
    pe_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_1_pe_merging.{graph_format}",
        )
    )
else:
    pe_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_1_pe_merging.{graph_format}",
        )
    )
time_print("PE graph generated.")

# Primer trimming
perc_kept_trim = pd.DataFrame(
    {
        "value": report_sheet_dict["4_primer_trimming"]["trimmed reads"]
        / report_sheet_dict["4_primer_trimming"]["processed reads"]
        * 100
    }
)
perc_kept_trim["sample"] = (
    report_sheet_dict["4_primer_trimming"]["File"]
    .str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace(".fastq.gz", "")
)  # type: ignore
perc_kept_trim = perc_kept_trim.fillna(0)
trim_graph = px.bar(
    perc_kept_trim.sort_values("value"),
    x="sample",
    y="value",
    range_y=[0, 100],
    labels={"value": "Reads kept [%]", "sample": "Samples"},
    title="Percentage of reads kept during trimming",
    width=graph_width,
)
trim_graph.add_hline(
    y=(sum(perc_kept_trim["value"]) / len(perc_kept_trim)),
    line_width=1,
    line_dash="dash",
)
trim_graph.update_layout(showlegend=False)
trim_graph.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    trim_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_2_trimming.{graph_format}",
        )
    )
else:
    trim_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_2_trimming.{graph_format}",
        )
    )
time_print("Trimming graph generated.")

# Read length distribution
time_print(
    "Reading in sequence files to generate read lengths graph. This can take a while..."
)
trimmed_seqs_dir = os.path.join(
    project_dir,
    "4_primer_trimming",
    "data",
)
# Make list of all samples
trimmed_seqs_files = [
    os.path.join(trimmed_seqs_dir, sample) for sample in os.listdir(trimmed_seqs_dir)
]

# Use the commandline to generate read length distribution file
## Construct the command
readlength_command = [
    "gunzip",
    "-c",
    *trimmed_seqs_files,
    "|",
    "awk",
    "'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'",
]
readlength_command = " ".join(readlength_command)

# Run the command and capture the output
readlength_result = subprocess.run(
    readlength_command, shell=True, stdout=subprocess.PIPE, text=True
)

# Convert the stdout to a Pandas DataFrame
readlength_result_lines = readlength_result.stdout.strip().split("\n")
readlength_result_data = [line.split() for line in readlength_result_lines]
readlength_df = pd.DataFrame(readlength_result_data, columns=["ReadLength", "Count"])
readlength_df = readlength_df.astype(int)

# Plot
length_graph = px.bar(
    readlength_df,
    x="ReadLength",
    y="Count",
    labels={"ReadLength": "Sequence length", "Count": "Frequency"},
    title="PE-merged read lengths before quality filtering",
)
length_graph.add_vline(
    x=min_length,
    line_width=1,
    line_dash="dash",
)
length_graph.add_vline(
    x=max_length,
    line_width=1,
    line_dash="dash",
)

# Save
if graph_format == "html":
    length_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_3_read_lengths.{graph_format}",
        )
    )
else:
    length_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_3_read_lengths.{graph_format}",
        )
    )
time_print("Read lengths graph generated.")

# Quality filtering
perc_kept_qf = pd.DataFrame(
    {
        "value": report_sheet_dict["5_quality_filtering"]["passed reads"]
        / report_sheet_dict["5_quality_filtering"]["processed reads"]
        * 100
    }
)
perc_kept_qf["sample"] = (
    report_sheet_dict["5_quality_filtering"]["File"]
    .str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace(".fasta.gz", "")
)  # type: ignore
perc_kept_qf = perc_kept_qf.fillna(0)
qf_graph = px.bar(
    perc_kept_qf.sort_values("value"),
    x="sample",
    y="value",
    range_y=[0, 100],
    labels={"value": "Reads kept [%]", "sample": "Sample"},
    title="Percentage of reads kept during quality filtering",
    width=graph_width,
)
qf_graph.add_hline(
    y=(sum(perc_kept_qf["value"]) / len(perc_kept_qf)),
    line_width=1,
    line_dash="dash",
)
qf_graph.update_layout(showlegend=False)
qf_graph.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    qf_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_4_qualityFiltering.{graph_format}",
        )
    )
else:
    qf_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_4_qualityFiltering.{graph_format}",
        )
    )
time_print("Quality filtering graph generated.")

# Dereplication
perc_kept_derep = pd.DataFrame(
    {
        "value": report_sheet_dict["6_dereplication"]["unique sequences"]
        / report_sheet_dict["6_dereplication"]["processed sequences"]
        * 100
    }
)
perc_kept_derep["sample"] = (
    report_sheet_dict["6_dereplication"]["File"]
    .str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace("_dereplicated", "")
    .str.replace(".fasta.gz", "")
)  # type: ignore
perc_kept_derep = perc_kept_derep.fillna(0)
derep_graph = px.bar(
    perc_kept_derep.sort_values("value"),
    x="sample",
    y="value",
    range_y=[0, 100],
    labels={"value": "Unique reads [%]", "sample": "Sample"},
    title="Number of unique reads per sample",
    width=graph_width,
)
derep_graph.add_hline(
    y=(sum(perc_kept_derep["value"]) / len(perc_kept_derep)),
    line_width=1,
    line_dash="dash",
)
derep_graph.update_layout(showlegend=False)
derep_graph.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    derep_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_5_dereplication.{graph_format}",
        )
    )
else:
    derep_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_5_dereplication.{graph_format}",
        )
    )
time_print("Dereplication graph generated.")

# Raw number of reads
num_reads_raw = pd.DataFrame(
    {"value": report_sheet_dict["3_PE merging"]["processed reads"]}
)
num_reads_raw["sample"] = (
    report_sheet_dict["3_PE merging"]["File"]
    .str.replace("_PE.fastq.gz", "")
    .str.replace(".fastq.gz", "")
)  # type: ignore
ymax_reads = max(num_reads_raw["value"])
rawreads_graph = px.bar(
    num_reads_raw.sort_values("value"),
    labels={"value": "Number of reads", "sample": "Sample"},
    title="Number of reads before PE merging, trimming and quality filtering",
    range_y=[0, ymax_reads],
    x="sample",
    y="value",
    width=graph_width,
)
rawreads_graph.add_hline(
    y=(sum(num_reads_raw["value"]) / len(num_reads_raw)), line_width=1, line_dash="dash"
)
rawreads_graph.update_layout(showlegend=False)
rawreads_graph.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    rawreads_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_6_numberRawreads.{graph_format}",
        )
    )
else:
    rawreads_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_6_numberRawreads.{graph_format}",
        )
    )
time_print("Raw reads graph generated.")

# Number of reads after PE merging and quality filtering
num_reads_filtered = pd.DataFrame(
    {"value": report_sheet_dict["5_quality_filtering"]["passed reads"]}
)
num_reads_filtered["sample"] = (
    report_sheet_dict["5_quality_filtering"]["File"]
    .str.replace("_PE", "")
    .str.replace("_trimmed", "")
    .str.replace("_filtered", "")
    .str.replace(".fasta.gz", "")
)  # type: ignore
filteredreads_graph = px.bar(
    num_reads_filtered.sort_values("value"),
    labels={"value": "Number of reads", "sample": "Sample"},
    title="Number of reads after PE merging, trimming and quality filtering",
    range_y=[0, ymax_reads],
    width=graph_width,
    x="sample",
    y="value",
)
filteredreads_graph.add_hline(
    y=(sum(num_reads_filtered["value"]) / len(num_reads_filtered)),
    line_width=1,
    line_dash="dash",
)
filteredreads_graph.update_layout(showlegend=False)
filteredreads_graph.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    filteredreads_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_7_numberFilteredreads.{graph_format}",
        )
    )
else:
    filteredreads_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_7_numberFilteredreads.{graph_format}",
        )
    )
time_print("Filtered reads graph generated.")

# Number of ESVs per sample before LULU
ymax_esvs = esv_prelulu_sums["count"].max()
prelulu_graph_esvs = px.bar(
    esv_prelulu_sums.sort_values("count"),
    labels={"count": "Number of ESVs", "index": "Sample"},
    title="Number of ESVs per sample before LULU filtering",
    range_y=[0, ymax_esvs],
    width=graph_width,
    y="count",
    x="index",
)
prelulu_graph_esvs.add_hline(
    y=(sum(esv_prelulu_sums["count"]) / len(esv_prelulu_sums)),
    line_width=1,
    line_dash="dash",
)
prelulu_graph_esvs.update_layout(showlegend=False)
prelulu_graph_esvs.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    prelulu_graph_esvs.write_html(
        os.path.join(
            outdir,
            f"{project_name}_8_prelulu_esvs.{graph_format}",
        )
    )
else:
    prelulu_graph_esvs.write_image(
        os.path.join(
            outdir,
            f"{project_name}_8_prelulu_esvs.{graph_format}",
        )
    )
time_print("Number of ESVs before LULU graph generated.")

# Number of ESVs per sample after LULU
postlulu_graph_esvs = px.bar(
    esv_postlulu_sums.sort_values("count"),
    labels={"count": "Number of ESVs", "index": "Sample"},
    title="Number of ESVs per sample after LULU filtering",
    range_y=[0, ymax_esvs],
    width=graph_width,
    y="count",
    x="index",
)
postlulu_graph_esvs.add_hline(
    y=(sum(esv_postlulu_sums["count"]) / len(esv_postlulu_sums)),
    line_width=1,
    line_dash="dash",
)
postlulu_graph_esvs.update_layout(showlegend=False)
postlulu_graph_esvs.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    postlulu_graph_esvs.write_html(
        os.path.join(
            outdir,
            f"{project_name}_9_postlulu_esvs.{graph_format}",
        )
    )
else:
    postlulu_graph_esvs.write_image(
        os.path.join(
            outdir,
            f"{project_name}_9_postlulu_esvs.{graph_format}",
        )
    )
time_print("Number of ESVs after LULU graph generated.")

# Number of OTUs per sample before LULU
ymax_otus = otu_prelulu_sums["count"].max()
prelulu_graph_otus = px.bar(
    otu_prelulu_sums.sort_values("count"),
    labels={"count": "Number of OTUs", "index": "Sample"},
    title="Number of OTUs per sample before LULU filtering",
    range_y=[0, ymax_otus],
    width=graph_width,
    y="count",
    x="index",
)
prelulu_graph_otus.add_hline(
    y=(sum(otu_prelulu_sums["count"]) / len(otu_prelulu_sums)),
    line_width=1,
    line_dash="dash",
)
prelulu_graph_otus.update_layout(showlegend=False)
prelulu_graph_otus.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    prelulu_graph_otus.write_html(
        os.path.join(
            outdir,
            f"{project_name}_10_prelulu_otus.{graph_format}",
        )
    )
else:
    prelulu_graph_otus.write_image(
        os.path.join(
            outdir,
            f"{project_name}_10_prelulu_otus.{graph_format}",
        )
    )
time_print("Number of OTUs before LULU graph generated.")

# Number of OTUs per sample after LULU
postlulu_graph_otus = px.bar(
    otu_postlulu_sums.sort_values("count"),
    labels={"count": "Number of OTUs", "index": "Sample"},
    title="Number of OTUs per sample after LULU filtering",
    range_y=[0, ymax_otus],
    width=graph_width,
    y="count",
    x="index",
)
postlulu_graph_otus.add_hline(
    y=(sum(otu_postlulu_sums["count"]) / len(otu_postlulu_sums)),
    line_width=1,
    line_dash="dash",
)
postlulu_graph_otus.update_layout(showlegend=False)
postlulu_graph_otus.update_xaxes(tickangle=55)
# Save
if graph_format == "html":
    postlulu_graph_otus.write_html(
        os.path.join(
            outdir,
            f"{project_name}_11_postlulu_otus.{graph_format}",
        )
    )
else:
    postlulu_graph_otus.write_image(
        os.path.join(
            outdir,
            f"{project_name}_11_postlulu_otus.{graph_format}",
        )
    )
time_print("Number of OTUs after LULU graph generated.")

# LULU pre post overview
ESVs = len(esv_prelulu_df)
ESVs_filtered = len(esv_postlulu_df)
OTUs = len(otu_prelulu_df)
OTUs_filtered = len(otu_postlulu_df)
lulu_pre_post_graph = go.Figure()
x_values = ["ESVs", "ESVs LULU-filtered", "OTUs", "OTUs LULU-filtered"]
y_values = [ESVs, ESVs_filtered, OTUs, OTUs_filtered]
text = [ESVs, ESVs_filtered, OTUs, OTUs_filtered]
lulu_pre_post_graph.add_trace(go.Bar(x=x_values, y=y_values, text=text))
lulu_pre_post_graph.update_layout(
    width=1000,
    height=800,
    title="Total number of ESVs and OTUs before and after LULU filtering",
)
lulu_pre_post_graph.update_traces(textposition="outside")
lulu_pre_post_graph.update_yaxes(title="OTUs/ESVs")
# Save
if graph_format == "html":
    lulu_pre_post_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_12_lulu_filtering_stats.{graph_format}",
        )
    )
else:
    lulu_pre_post_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_12_lulu_filtering_stats.{graph_format}",
        )
    )
time_print("Pre- vs. post-LULU comparison graph generated.")

# Number of reads vs number of ESVs
reads_esvs_graph = px.scatter(
    pd.concat([num_reads_filtered, esv_postlulu_sums["count"]], axis=1).reset_index(),
    x="value",
    y="count",
    trendline="ols",
    trendline_color_override="black",
    hover_name="sample",
    title="Filtered reads vs. ESVs after LULU filtering",
    width=480,
    labels={
        "count": "Number of ESVs after LULU filtering",
        "value": "Number of filtered reads",
    },
)
# Save
if graph_format == "html":
    reads_esvs_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_13_reads_vs_esvs.{graph_format}",
        )
    )
else:
    reads_esvs_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_13_reads_vs_esvs.{graph_format}",
        )
    )
time_print("Number of reads vs. ESVs graph generated.")


# Number of reads vs number of OTUs
reads_otus_graph = px.scatter(
    pd.concat([num_reads_filtered, otu_postlulu_sums["count"]], axis=1).reset_index(),
    x="value",
    y="count",
    trendline="ols",
    trendline_color_override="black",
    hover_name="sample",
    title="Filtered reads vs. OTUs after LULU filtering",
    width=480,
    labels={
        "count": "Number of OTUs after LULU filtering",
        "value": "Number of filtered reads",
    },
)
# Save
if graph_format == "html":
    reads_otus_graph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_14_reads_vs_otus.{graph_format}",
        )
    )
else:
    reads_otus_graph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_14_reads_vs_otus.{graph_format}",
        )
    )
time_print("Number of reads vs. OTUs graph generated.")

# Add Accumulation curves - ESVs, OTUs, ESVs species, OTUs species


# Lineplot ESVs
esv_linegraph = go.Figure()
for sample in samples:
    y_values = df_stats.loc[sample].values.tolist()[:-2] + [
        df_stats.loc[sample].values.tolist()[-1]
    ]
    x_values = df_stats.columns.tolist()[:-2] + [df_stats.columns.tolist()[-1]]
    esv_linegraph.add_trace(go.Scatter(x=x_values, y=y_values, name=sample))
esv_linegraph.update_layout(
    template="simple_white",
    width=1000,
    height=800,
    title="ESVs - Reads per sample for each processing step",
)
esv_linegraph.update_yaxes(title="Reads")
esv_linegraph = esv_linegraph.update_layout(showlegend=False)
# Save
if graph_format == "html":
    esv_linegraph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_15_linegraph_esvs.{graph_format}",
        )
    )
else:
    esv_linegraph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_15_linegraph_esvs.{graph_format}",
        )
    )
time_print("Lineplot ESVs graph generated.")

# Lineplot OTUs
otu_linegraph = go.Figure()
for sample in samples:
    y_values = df_stats.loc[sample].values.tolist()[:-1]
    x_values = df_stats.columns.tolist()[:-1]
    otu_linegraph.add_trace(go.Scatter(x=x_values, y=y_values, name=sample))
otu_linegraph.update_layout(
    template="simple_white",
    width=1000,
    height=800,
    title="OTUs - Reads per sample for each processing step",
)
otu_linegraph.update_yaxes(title="Reads")
otu_linegraph = otu_linegraph.update_layout(showlegend=False)
# Save
if graph_format == "html":
    otu_linegraph.write_html(
        os.path.join(
            outdir,
            f"{project_name}_16_linegraph_otus.{graph_format}",
        )
    )
else:
    otu_linegraph.write_image(
        os.path.join(
            outdir,
            f"{project_name}_16_linegraph_otus.{graph_format}",
        )
    )
time_print("Lineplot OTUs graph generated.")

# Boxplot
boxplot = go.Figure()
for category in df_stats.columns.tolist():
    y_values = df_stats.loc[samples][category].values.tolist()
    boxplot.add_trace(go.Box(y=y_values, name=category))
boxplot.update_layout(
    width=1000,
    height=800,
    title="Read number overview for each processing step",
    showlegend=False,
)
boxplot.update_yaxes(title="Reads")
# Save
if graph_format == "html":
    boxplot.write_html(
        os.path.join(
            outdir,
            f"{project_name}_17_boxplot_read_summary.{graph_format}",
        )
    )
else:
    boxplot.write_image(
        os.path.join(
            outdir,
            f"{project_name}_17_boxplot_read_summary.{graph_format}",
        )
    )
time_print("Boxplot generated.")

# ESV clustergram
## Format
ID_list_esv = esv_postlulu_df["ID"].values.tolist()
esv_postlulu_samples_df = esv_postlulu_df.drop(columns=["ID", "Seq"])
sample_list = esv_postlulu_samples_df.columns.tolist()
## Take log
esv_postlulu_samples_log_df = np.where(
    esv_postlulu_samples_df != 0, np.log(esv_postlulu_samples_df), 0
)
## Define colours
colors = px.colors.sample_colorscale("plasma", [n / 5 for n in range(5)])
## Create clustergram
time_print("Generating clustergram for ESVs. This can take a while...")
esv_clustergram = dash_bio.Clustergram(
    data=esv_postlulu_samples_log_df,
    column_labels=list(esv_postlulu_samples_df.columns.values),
    row_labels=list(esv_postlulu_samples_df.index),
    height=800,
    width=min(graph_width, 3000),
    hidden_labels="row",
    color_map=[
        [0.0, colors[0]],
        [0.25, colors[1]],
        [0.5, colors[2]],
        [0.75, colors[3]],
        [1.0, colors[4]],
    ],
    paper_bg_color="white",
)
esv_clustergram.update_layout(
    title="Clustergram with log-transformed ESV abundances",
)
# Save
if graph_format == "html":
    esv_clustergram.write_html(
        os.path.join(
            outdir,
            f"{project_name}_18_esv_clustergram.{graph_format}",
        )
    )
else:
    esv_clustergram.write_image(
        os.path.join(
            outdir,
            f"{project_name}_18_esv_clustergram.{graph_format}",
        )
    )
time_print("Clustergram generated for ESVs.")

# OTU clustergram
## Format
ID_list_otu = otu_postlulu_df["ID"].values.tolist()
otu_postlulu_samples_df = otu_postlulu_df.drop(columns=["ID", "Seq"])
sample_list = otu_postlulu_samples_df.columns.tolist()
## Take log
otu_postlulu_samples_log_df = np.where(
    otu_postlulu_samples_df != 0, np.log(otu_postlulu_samples_df), 0
)
## Define colours
colors = px.colors.sample_colorscale("plasma", [n / 5 for n in range(5)])
## Create clustergram
time_print("Generating clustergram for OTUs. This can take a while...")
otu_clustergram = dash_bio.Clustergram(
    data=otu_postlulu_samples_log_df,
    column_labels=list(otu_postlulu_samples_df.columns.values),
    row_labels=list(otu_postlulu_samples_df.index),
    height=800,
    width=min(graph_width, 3000),
    hidden_labels="row",
    color_map=[
        [0.0, colors[0]],
        [0.25, colors[1]],
        [0.5, colors[2]],
        [0.75, colors[3]],
        [1.0, colors[4]],
    ],
    paper_bg_color="white",
)
otu_clustergram.update_layout(
    title="Clustergram with log-transformed OTU abundances",
)
# Save
if graph_format == "html":
    otu_clustergram.write_html(
        os.path.join(
            outdir,
            f"{project_name}_19_otu_clustergram.{graph_format}",
        )
    )
else:
    otu_clustergram.write_image(
        os.path.join(
            outdir,
            f"{project_name}_19_otu_clustergram.{graph_format}",
        )
    )
time_print("Clustergram generated for OTUs.")

# Kronagraphs
if blast == "True":  # Requirement as we need taxonomic information for Kronagraphs
    time_print("Generating kronagraphs...")
    # Format dfs for Krona
    esv_krona_df = krona_formatting(esv_final_df)
    otu_krona_df = krona_formatting(otu_final_df)
    # Save so that krona can be run in command line
    esv_krona_df.to_csv(
        os.path.join(
            outdir,
            f"{project_name}_ESVs_krona-formatted.csv",
        ),
        header=False,
        index=False,
        sep="\t",
    )
    otu_krona_df.to_csv(
        os.path.join(
            outdir,
            f"{project_name}_OTUs_krona-formatted.csv",
        ),
        header=False,
        index=False,
        sep="\t",
    )
    # Use the commandline to run krona on both
    ## Construct the krona command for ESVs
    krona_command_esvs = " ".join(
        [
            "ktImportText",
            "-o",
            os.path.join(outdir, f"{project_name}_20_esv_krona.html"),
            os.path.join(outdir, f"{project_name}_ESVs_krona-formatted.csv"),
        ]
    )
    ## Run the command
    subprocess.call(krona_command_esvs, shell=True)
    ## Construct the krona command for OTUs
    krona_command_otus = " ".join(
        [
            "ktImportText",
            "-o",
            os.path.join(outdir, f"{project_name}_20_otu_krona.html"),
            os.path.join(outdir, f"{project_name}_OTUs_krona-formatted.csv"),
        ]
    )
    ## Run the command
    subprocess.call(krona_command_otus, shell=True)
    # Remove formatted files
    os.remove(os.path.join(outdir, f"{project_name}_ESVs_krona-formatted.csv"))
    os.remove(os.path.join(outdir, f"{project_name}_OTUs_krona-formatted.csv"))

time_print("Finished graph generation.")
