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
    print(datetime.datetime.now().strftime("%H:%M:%S"), ": ", text, sep="")


# Define a custom argument type for a list of strings
def list_of_strings(arg):
    return arg.split(",")


# Export cleaned sequences as FASTA file
# Function to write FASTA file
def write_fasta(df, filename):
    with open(filename, "w") as file:
        for index, row in df.iterrows():
            file.write(f'>{row["ID"]}\n{row["Seq"]}\n')


# Define function to process dfs
def remove_negs_from_df(df, unit, negative_controls):
    # Identify negative controls that do and don't contain 0 reads (microDecon gives an error if used with negative controls with 0 reads)
    negative_controls_keep = [neg for neg in negative_controls if df[neg].sum() > 0]
    negative_controls_drop = [
        neg for neg in negative_controls if neg not in negative_controls_keep
    ]

    # If all negative controls are removed because they don't contain any reads, return the df without neg controls
    if not negative_controls_keep:
        return df.drop(columns=negative_controls_drop)

    # Identify true samples and all samples
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    true_samples = list(
        df.columns.difference(negative_controls + ["ID", "Seq", "total_reads"] + ranks)
    )
    samples = negative_controls_keep + true_samples

    # Generate sample df and format for microDecon
    df_samples = df[["ID"] + samples]

    # Separating other information and dropping negative controls with 0 reads
    df_other = df.drop(columns=samples + negative_controls_drop)

    # If only 1 OTU/ESV in the blank(s) contains reads, microDecon will fail.
    # In that case, we instead just subtract the summed reads from the 1 OTU/ESV from the samples.
    non_zero_rows = df_samples[negative_controls_keep][
        (df_samples[negative_controls_keep] != 0).any(axis=1)
    ]
    if len(non_zero_rows) == 1:
        # Get the row sum
        non_zero_rows_sum = df_samples[negative_controls_keep].sum(axis=1)

        # Subtract from samples
        df_samples_decon = df_samples[true_samples].sub(non_zero_rows_sum, axis=0)

        # Turn negative values to 0
        df_samples_decon = df_samples_decon.applymap(lambda x: max(0, x))

        # Add ID column
        df_samples_decon = pd.DataFrame(
            {"ID": df["ID"], **df_samples_decon.to_dict("list")}
        )

    # Otherwise we use microDecon
    else:
        # Remove negatives
        decon_results = microDecon(
            df_samples,
            numb_blanks=len(negative_controls_keep),
            numb_ind=robjects.IntVector([len(true_samples)]),
            taxa=False,
            thresh=1,  # Turns off to drop OTUs that show up in only 70% of samples
            prop_thresh=0,  # Turns off to drop OTUs with a total abundance of < 0.005%
        )

        # Convert the result to a Python object and save filtered df
        df_samples_decon = robjects.conversion.rpy2py(decon_results.rx2("decon.table"))

        # Delete negative control column
        df_samples_decon = df_samples_decon.drop(df_samples_decon.columns[1], axis=1)

    # Merge reads and other data
    df_decon = df_samples_decon.merge(df_other, on="ID")

    # Restore row order
    df_decon["NumericValues"] = df_decon["ID"].str.replace(f"{unit}_", "").astype(int)
    df_decon = df_decon.sort_values(by="NumericValues")
    df_decon = df_decon.drop(columns=["NumericValues"])

    # Drop ESVs with 0 reads
    df_decon = df_decon[df_decon.drop(columns=["ID", "Seq"]).sum(axis=1) != 0]

    return df_decon


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

# Start of pipeline
time_print(
    "Removing negative control reads from samples with the R package microDecon..."
)

# Define path name variables
path_to_otu_clustering = os.path.join(
    project_name, "9_lulu_filtering", "otu_clustering"
)
path_to_denoising = os.path.join(project_name, "9_lulu_filtering", "denoising")

# Path to ESV/OTU post-LULU files
otu_postlulu_file = os.path.join(
    path_to_otu_clustering,
    f"{project_name}_OTU_table-with_filtered_taxonomy.csv",
)
esv_postlulu_file = os.path.join(
    path_to_denoising,
    f"{project_name}_ESV_table-with_filtered_taxonomy.csv",
)

esv_postlulu_df = pd.read_csv(esv_postlulu_file)
otu_postlulu_df = pd.read_csv(otu_postlulu_file)

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

# Process dfs
otu_postlulu_df_microdeconFiltered = remove_negs_from_df(
    otu_postlulu_df, "OTU", args.negative_controls
)
esv_postlulu_df_microdeconFiltered = remove_negs_from_df(
    esv_postlulu_df, "ESV", args.negative_controls
)

# Export dfs
otu_postlulu_df_microdeconFiltered.to_csv(
    otu_postlulu_file.replace(".csv", "-without_NegControls.csv"),
    index=False,
)
esv_postlulu_df_microdeconFiltered.to_csv(
    esv_postlulu_file.replace(".csv", "-without_NegControls.csv"),
    index=False,
)

# Provide the desired filenames for the FASTA files
fasta_filename_otus = os.path.join(
    path_to_otu_clustering,
    f"{project_name}-OTU_sequences-without_NegControls.fasta",
)
fasta_filename_esvs = os.path.join(
    path_to_denoising,
    f"{project_name}-ESV_sequences-without_NegControls.fasta",
)

# Write the FASTA files
write_fasta(otu_postlulu_df_microdeconFiltered, fasta_filename_otus)
write_fasta(esv_postlulu_df_microdeconFiltered, fasta_filename_esvs)

time_print("Negative control reads removed.")
