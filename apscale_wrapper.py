#!/usr/bin/env python3

"""
A wrapper to run apscale on forward and reverse reads and to generate 
processing QC graphs.

Requires the submodules apscale_processing_graphs.py and apscale_subtract_negatives.py in PATH.

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
    print(datetime.datetime.now().strftime("%H:%M:%S"), "---", text)


# Function to generate an apscale Settings.xlsx file
def generateSettings(**kwargs):
    with pd.ExcelWriter(
        Path(f"{args.project_name}_apscale").joinpath("Settings.xlsx"),
        mode="w",
        engine="openpyxl",
    ) as writer:
        # ## write the 0_general_settings sheet
        df_0 = pd.DataFrame(
            [[args.cores, 9]],
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
            [[args.forward_primer, args.reverse_primer, "False"]],
            columns=["P5 Primer (5' - 3')", "P7 Primer (5' - 3')", "anchoring"],
        )

        df_4.to_excel(writer, sheet_name="4_primer_trimming", index=False)

        ## write the 5_quality_filtering sheet
        df_5 = pd.DataFrame(
            [[args.maxEE, args.min_length, args.max_length]],
            columns=["maxEE", "min length", "max length"],
        )

        df_5.to_excel(writer, sheet_name="5_quality_filtering", index=False)

        ## write the 6_dereplication_pooling sheet
        df_6 = pd.DataFrame([[4]], columns=["min size to pool"])

        df_6.to_excel(writer, sheet_name="6_dereplication_pooling", index=False)

        ## write the 7_otu_clustering sheet
        df_7 = pd.DataFrame(
            [
                [
                    args.clusteringtool,
                    args.otu_perc,
                    args.swarm_distance,
                    args.prior_denoising,
                    args.coi,
                    "True",
                ]
            ],
            columns=[
                "clustering tool",
                "vsearch pct id",
                "swarm distance",
                "prior denoise",
                "coi",
                "to excel",
            ],
        )

        df_7.to_excel(writer, sheet_name="7_otu_clustering", index=False)

        ## write the 8_denoising sheet
        df_8 = pd.DataFrame(
            [[2, 8, args.coi, "True"]], columns=["alpha", "minsize", "coi", "to excel"]
        )

        df_8.to_excel(writer, sheet_name="8_denoising", index=False)

        ## write the 9_lulu_filtering sheet
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
    description="""A wrapper to run apscale on forward and reverse 
	reads and to generate processing QC graphs.
	Requires the submodule apscale_processing_graphs.py in PATH.""",
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
    help="OTU identify treshold for clustering with vsearch. Only used when --clusteringtool=vsearch (default=97).",
)
parser.add_argument(
    "-c",
    "--coi",
    choices=["False", "True"],
    default="False",
    help="""Are you processing COI data? If yes, the fact that COI is a coding gene can be used to
    improve denoising (default=False).""",
)
parser.add_argument(
    "-d",
    "--swarm_distance",
    default=1,
    metavar="N",
    type=int,
    help="Distance used for swarm. Overwritten to 13 if --coi=True. Only used when --clusteringtool=swarm (default: 1).",
)
parser.add_argument(
    "-t",
    "--clusteringtool",
    default="vsearch",
    choices=["vsearch", "swarm"],
    help="Tool used for OTU clustering (default=vsearch).",
)
parser.add_argument(
    "-D",
    "--prior_denoising",
    choices=["False", "True"],
    default="False",
    help="Set to True if you want to denoise reads prior to OTU clustering (default: False).",
)
parser.add_argument(
    "-e",
    "--maxEE",
    metavar="N",
    default=2,
    type=int,
    help="maxEE (maximum estimated error) value used for quality filtering (default: 2).",
)
parser.add_argument(
    "-n",
    "--cores",
    metavar="N",
    default=2,
    type=int,
    help="Number of cores to use (default: 2).",
)
parser.add_argument(
    "-g",
    "--graph_format",
    help="Graph format, either HTML, png, svg. HTML is recommended because it creates interactive plots (default: html).",
    default="html",
    choices=["png", "svg", "html"],
)
parser.add_argument(
    "-R",
    "--remove_negatives",
    help="Do you want to remove reads in negative controls from the samples using the R package microDecon? (default: False)",
    default="False",
    choices=["True", "False"],
)
parser.add_argument(
    "-N",
    "--negatives",
    help="Required if --removes_negatives=True. List the names of all negative controls (without _R1/2.fq.gz), separated by commas without spaces.",
    metavar="control1,control2,control3",
    type=str,
)
parser.add_argument(
    "-S",
    "--scaling_factor",
    help="Scaling factor for graph width. Manual trial and error in 0.2 increments might be required (default: 1.0).",
    default=1.0,
    metavar="N.N",
    type=float,
)


# Define a custom validation function to enforce the requirement for --remove_negatives
def validate_args(args):
    if args.remove_negatives == "True" and not args.negatives:
        parser.error("--negatives is required when --remove_negatives=True.")


# Register the custom validation function to be called after parsing arguments
parser.set_defaults(func=validate_args)

# Parse arguments
args = parser.parse_args()

# Call the custom validation function to check the requirements
args.func(args)

### Start of pipeline
time_print("Starting apscale wrapper.")

# Create an apscale directory using bash
time_print("Creating apscale directory...")
subprocess.run(["apscale", "--create_project", args.project_name])

# Create an empty Project_report.xlsx file
## Create an ExcelWriter object using the openpyxl engine
excel_writer = pd.ExcelWriter(
    f"{args.project_name}_apscale/Project_report.xlsx", engine="openpyxl"
)
## Write an empty DataFrame to the Excel file
pd.DataFrame().to_excel(excel_writer, sheet_name="Sheet1", index=False)
# Save the Excel file
excel_writer.book.save(f"{args.project_name}_apscale/Project_report.xlsx")

# Generate symlinks to demultiplexed reads
time_print("Generating symlinks to demultiplexed reads...")
subprocess.run(
    f'ln -s "$(realpath "{args.sequence_dir}")"/* {args.project_name}_apscale/2_demultiplexing/data',
    shell=True,
)

# Generate a Settings.xlsx file with given parameters
time_print("Generating apscale settings file...")
generateSettings(
    project_name=args.project_name,
    forward_primer=args.forward_primer,
    reverse_primer=args.reverse_primer,
    min_length=args.min_length,
    max_length=args.max_length,
    otu_perc=args.otu_perc,
    maxEE=args.maxEE,
    cores=args.cores,
    coi=args.coi,
    clusteringtool=args.clusteringtool,
    swarm_distance=args.swarm_distance,
    prior_denoising=args.prior_denoising,
)

# Run apscale
time_print("Starting apscale...")
subprocess.run(["apscale", "--run_apscale", f"{args.project_name}_apscale"])
time_print("Apscale done.")

if args.remove_negatives == "True":
    subprocess.run(
        [
            "apscale_subtract_negatives.py",
            "--project_dir",
            f"{args.project_name}_apscale",
            "--negatives",
            f"{args.negatives}",
        ]
    )

# Generate processing graphs using separate script
time_print("Generating apscale processing graphs...")
subprocess.run(
    [
        "apscale_processing_graphs.py",
        "--project_dir",
        f"{args.project_name}_apscale",
        "--graph_format",
        f"{args.graph_format}",
        "--scaling_factor",
        f"{args.scaling_factor}",
    ]
)

time_print("Apscale wrapper done.")
