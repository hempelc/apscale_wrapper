#!/usr/bin/env python3

"""
A submodule for the apscale wrapper to filter reads in negative controls from
samples with microDecon once apscale is done.

By Chris Hempel (christopher.hempel@kaust.edu.sa) on 1 Nov 2023
"""

import datetime
import os
import pandas as pd
import argparse
import warnings
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri

# Define that warnings are not printed to console
warnings.filterwarnings("ignore")


# Funtion to print datetime and text
def time_print(text):
    print(datetime.datetime.now().strftime("%H:%M:%S"), "---", text)


# Define a custom argument type for a list of strings
def list_of_strings(arg):
    return arg.split(",")


# Define arguments
parser = argparse.ArgumentParser(
    description="""A submodule for the apscale wrapper to filter reads in negative controls from
    samples with microDecon once apscale is done.""",
)
parser.add_argument(
    "--negative_controls",
    help="List the names of all negative controls (without _R1/2.fq.gz), separated by commas without spaces.",
    metavar="control1,control2,control3",
    required=True,
    type=list_of_strings,
)
parser.add_argument(
    "--project_dir",
    help="Directory containing apscale results.",
    required=True,
)

# Parse argument
args = parser.parse_args()

# Set project_name argument
project_name = os.path.basename(args.project_dir)

### Start of pipeline
time_print(
    "Removing negative control reads from samples with the R package microDecon..."
)

# Path to ESV/OTU post-LULU files
otu_postlulu_file = os.path.join(
    args.project_dir,
    "9_lulu_filtering",
    "otu_clustering",
    f"{project_name}_OTU_table_filtered.parquet.snappy",
)
esv_postlulu_file = os.path.join(
    args.project_dir,
    "9_lulu_filtering",
    "denoising",
    f"{project_name}_ESV_table_filtered.parquet.snappy",
)

esv_postlulu_df = pd.read_parquet(esv_postlulu_file, engine="fastparquet")
otu_postlulu_df = pd.read_parquet(otu_postlulu_file, engine="fastparquet")

# Test if microdecon is installed:
if not rpackages.isinstalled("microDecon"):
    time_print(
        "The required R package microDecon is not installed. Installing now. This can take a while..."
    )
    # Import R utils to install other packages
    utils = rpackages.importr("utils")
    # Set mirror
    utils.chooseCRANmirror(ind=1)
    if not rpackages.isinstalled("devtools"):
        time_print(
            "The R package devtools is required to install microDecon. Installing devtools. This can take a while..."
        )
        utils.install_packages("devtools")
    devtools = rpackages.importr("devtools")
    devtools.install_github("donaldtmcknight/microDecon")
    time_print(
        "Installation of microDecon finished. Removing negative control reads from samples..."
    )

# Load microDecon
robjects.r('library("microDecon")')
# Define the decon function
microDecon = robjects.r["decon"]
# Activate that pandas can be converted to R
pandas2ri.activate()


# Define function to process dfs
def remove_negs_from_df(df, unit, negative_controls):
    # Identify samples and neg controls
    true_samples = list(df.columns.difference(negative_controls + ["ID", "Seq"]))
    samples = negative_controls + true_samples

    # Generate sample df and format for microDecon
    df_samples = df[["ID"] + samples]

    # Separating other information
    df_other = df.drop(columns=samples)

    # Remove negatives
    decon_results = microDecon(
        df_samples,
        numb_blanks=len(negative_controls),
        numb_ind=robjects.IntVector([len(true_samples)]),
        taxa=False,
    )

    # Convert the result to a Python object and save filtered df
    df_samples_decon = robjects.conversion.rpy2py(decon_results.rx2("decon.table"))

    # Merge reads and other data
    df_decon = df_samples_decon.merge(df_other, on="ID")

    # Restore row order
    df_decon["NumericValues"] = df_decon["ID"].str.replace(f"{unit}_", "").astype(int)
    df_decon = df_decon.sort_values(by="NumericValues")
    df_decon = df_decon.drop(columns=["NumericValues"])

    # Delete negative control column
    df_decon = df_decon.drop(df_decon.columns[1], axis=1)
    return df_decon


# Process dfs
otu_postlulu_df_microdeconFiltered = remove_negs_from_df(
    otu_postlulu_df, "OTU", args.negative_controls
)
esv_postlulu_df_microdeconFiltered = remove_negs_from_df(
    esv_postlulu_df, "ESV", args.negative_controls
)

# Export dfs
otu_postlulu_df_microdeconFiltered.to_csv(
    otu_postlulu_file.replace(".parquet.snappy", "_microdecon-filtered.csv"),
    index=False,
)
esv_postlulu_df_microdeconFiltered.to_csv(
    esv_postlulu_file.replace(".parquet.snappy", "_microdecon-filtered.csv"),
    index=False,
)


# Export cleaned sequences as FASTA file
# Function to write FASTA file
def write_fasta(df, filename):
    with open(filename, "w") as file:
        for index, row in df.iterrows():
            file.write(f'>{row["ID"]}\n{row["Seq"]}\n')


# Provide the desired filenames for the FASTA files
fasta_filename_otus = os.path.join(
    args.project_dir,
    "9_lulu_filtering",
    "otu_clustering",
    f"{project_name}_OTUs_filtered_microdecon-filtered.fasta",
)
fasta_filename_esvs = os.path.join(
    args.project_dir,
    "9_lulu_filtering",
    "otu_clustering",
    f"{project_name}_ESVs_filtered-microdecon-filtered.fasta",
)

# Write the FASTA files
write_fasta(otu_postlulu_df_microdeconFiltered, fasta_filename_otus)
write_fasta(esv_postlulu_df_microdeconFiltered, fasta_filename_esvs)
