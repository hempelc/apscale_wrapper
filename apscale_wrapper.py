#!/usr/bin/env python3

"""
A wrapper to run apscale on forward and reverse reads, generate processing QC graphs, and assign taxonomy to generated taxonomic units.

Requires the submodules settings_generator.py and XXX in path.

By Chris Hempel (christopher.hempel@kaust.edu.sa) on 15 Aug 2023
"""

import datetime
import pandas as pd
import argparse
import warnings
import subprocess
from pathlib import Path

# Define that warnings are not printed to console
warnings.filterwarnings("ignore")


# Funtion to print datetime and text
def time_print(text):
    datetime_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"{datetime_now}  ---  " + text)


# Function to generate an apscale Settings.xlsx file
def generateSettings(**kwargs):
    with pd.ExcelWriter(
        Path(f"{project_name}_apscale").joinpath("Settings.xlsx"),
        mode="w",
        engine="openpyxl",
    ) as writer:
        # ## write the 0_general_settings sheet
        df_0 = pd.DataFrame(
            [[cores, 9]],
            columns=["cores to use", "compression level"],
        )

        df_0.to_excel(writer, sheet_name="0_general_settings", index=False)

        ## write the 3_PE_merging sheet
        df_3 = pd.DataFrame(
            [[25, 199, 5]], columns=["maxdiffpct", "maxdiffs", "minovlen"]
        )

        df_3.to_excel(writer, sheet_name="3_PE_merging", index=False)

        ## write the 4_primer_trimming sheet
        df_4 = pd.DataFrame(
            [[forward_primer, reverse_primer, "False"]],
            columns=["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", "anchoring"],
        )

        df_4.to_excel(writer, sheet_name="4_primer_trimming", index=False)

        ## write the 5_quality_filtering sheet
        df_5 = pd.DataFrame(
            [[maxEE, min_length, max_length]],
            columns=["maxEE", "min length", "max length"],
        )

        df_5.to_excel(writer, sheet_name="5_quality_filtering", index=False)

        ## write the 6_dereplication_pooling sheet
        df_6 = pd.DataFrame([[4]], columns=["min size to pool"])

        df_6.to_excel(writer, sheet_name="6_dereplication_pooling", index=False)

        ## write the 7_otu_clustering sheet
        df_7 = pd.DataFrame([[otu_perc, "True"]], columns=["pct id", "to excel"])

        df_7.to_excel(writer, sheet_name="7_otu_clustering", index=False)

        ## write the 8_denoising sheet
        df_8 = pd.DataFrame([[2, 8, "True"]], columns=["alpha", "minsize", "to excel"])

        df_8.to_excel(writer, sheet_name="8_denoising", index=False)

        ## write the 8_denoising sheet
        df_9 = pd.DataFrame(
            [[84, 95, 1, "True"]],
            columns=[
                "minimum similarity",
                "minimum relative cooccurence",
                "minimum ratio",
                "to excel",
            ],
        )

        df_9.to_excel(writer, sheet_name="9_lulu_filtering", index=False)


# Define arguments
parser = argparse.ArgumentParser(
    description="""A wrapper to run apscale on forward and reverse reads, generate processing QC graphs,
    and assign taxonomy to generated taxonomic units.""",
)
parser.add_argument(
    "-s",
    "--sequence_dir",
    metavar="PATH/TO/SEQUENCE_DIR",
    help="""Path to directory containing demultiplexed forward and reverse sequences in .fastq.gz format.
    Note: file names must end with R1 and R2 to identify read pairs.""",
    required=True,
)
parser.add_argument(
    "-p",
    "--project_name",
    help="Name of the apscale project to be generated.",
    required=True,
)
parser.add_argument(
    "-f",
    "--forward_primer",
    metavar="PRIMER_SEQUENCE",
    help="Forward primer sequence to trim.",
    required=True,
)
parser.add_argument(
    "-r",
    "--reverse_primer",
    metavar="PRIMER_SEQUENCE",
    help="Reverse primer sequence to trim.",
    required=True,
)
parser.add_argument(
    "-m",
    "--min_length",
    metavar="NNN",
    type=int,
    help="Minimum limit of expected amplicon length (used for length filtering).",
    required=True,
)
parser.add_argument(
    "-M",
    "--max_length",
    metavar="NNN",
    type=int,
    help="Maximum limit of expected amplicon length (used for length filtering).",
    required=True,
)
parser.add_argument(
    "-o",
    "--otu_perc",
    metavar="NN",
    default=97,
    type=int,
    help="OTU identify treshold for clustering (default=97).",
)
parser.add_argument(
    "-e",
    "--maxEE",
    metavar="N",
    default=2,
    type=int,
    help="maxEE (maximum estimated error) value used for quality filtering (default=2).",
)
parser.add_argument(
    "-n",
    "--cores",
    metavar="N",
    default=2,
    type=int,
    help="Number of cores to use (default=2).",
)
parser.add_argument(
    "-g",
    "--graph_format",
    help="Format for processing report graphs, either png or svg.",
    required=True,
    choices=["png", "svg"],
)
parser.add_argument(
    "-S",
    "--scaling_factor",
    help="Scaling factor for graph width. Manual trial and error in 0.2 increments might be required (default: 1).",
    default=1,
    type=int,
)
args = parser.parse_args()

# Set arguments
sequence_dir = args.sequence_dir
project_name = args.project_name
forward_primer = args.forward_primer
reverse_primer = args.reverse_primer
min_length = args.min_length
max_length = args.max_length
graph_format = args.graph_format
otu_perc = args.otu_perc
maxEE = args.maxEE
cores = args.cores
scaling_factor = args.scaling_factor

### Start of pipeline
time_print("Starting apscale wrapper.")

# Create an apscale directory using bash
time_print("Creating apscale directory...")
subprocess.run(["apscale", "--create_project", project_name])

# Generate symlinks to demultiplexed reads
time_print("Generating symlinks to demultiplexed reads...")
subprocess.run(
    f"ln -s $(realpath {sequence_dir})/* {project_name}_apscale/2_demultiplexing/data",
    shell=True,
)

# Generate a Settings.xlsx file with given parameters
time_print("Generating apscale settings file...")
generateSettings(
    project_name=project_name,
    forward_primer=forward_primer,
    reverse_primer=reverse_primer,
    min_length=min_length,
    max_length=max_length,
    otu_perc=otu_perc,
    maxEE=maxEE,
    cores=cores,
)

# Run apscale
time_print("Starting apscale...")
subprocess.run(["apscale", "--run_apscale", f"{project_name}_apscale"])
time_print("Apscale done.")


# Generate processing graphs using separate script
time_print("Generating apscale processing graphs...")
subprocess.run(
    [
        "apscale_processing_graphs.py",
        "--project_name",
        f"{project_name}",
        "--graph_format",
        f"{graph_format}",
        "--scaling_factor",
        f"{scaling_factor}",
    ]
)

time_print("Apscale wrapper done.")
