#!/usr/bin/env python3

"""
A submodule for the abscale wrapper to run BLAST on generated FASTA files and filter BLAST hits.

By Chris Hempel (christopher.hempel@kaust.edu.sa) on 20 Jan 2022
"""

import datetime
import pandas as pd
import argparse
import warnings
import subprocess
import copy
import os
import numpy as np


# Define that warnings are not printed to console
warnings.filterwarnings("ignore")


# Funtion to print datetime and text
def time_print(text):
    print(datetime.datetime.now().strftime("%H:%M:%S"), ": ", text, sep="")


# Function to filter bitscores
def bitscore_cutoff(x, percentage):
    min_bitscore = x.max() - x.max() * percentage
    return x[x >= min_bitscore]


# Define a custom argument type for a list of strings
def list_of_strings(arg):
    return arg.split(",")


# Define a custom argument type for a list of integers
def list_of_integers(arg):
    return [int(value) for value in arg.split(",")]


# Function to replace unknown taxonomy in the PR2 database
def replace_pr2_tax(taxon):
    if taxon.endswith("X") or taxon.endswith("sp."):
        return "Unknown in PR2 database"
    else:
        return taxon


# Function to add suffix to filename despite format
def add_suffix(filename, suffix="_no_cutoff"):
    root, extension = os.path.splitext(filename)
    return f"{root}{suffix}{extension}"


# Define function to process tax dfs
def post_processing(df):
    # Keep only relevant columns and put species to last column
    df_tmp = df[["qseqid"] + ranks]

    ## Make a df mask: group dfs, check if ranks have more than one taxon, and if yes, True, else False
    lca_mask = df_tmp.groupby(["qseqid"]).transform(lambda x: len(set(x)) != 1)

    ## Replace ranks in df with "NA" based on mask
    df_tmp = df_tmp.mask(lca_mask, "NA")

    # Add qseqid and pident info
    df_tmp["qseqid"] = df["qseqid"]
    df_tmp["pident"] = df["pident"]

    ## Drop duplicate rows == aggregate taxonomic info
    df_processed = df_tmp.drop_duplicates()

    # Per ID and identical taxon match, keep the max pident
    idx_pident = (
        df_processed.groupby(["qseqid"] + ranks)["pident"].transform(max)
        == df_processed["pident"]
    )
    df_processed = df_processed[idx_pident]
    df_processed = df_processed.rename(columns={"pident": "percentage_similarity"})

    # Change column name, fill NaNs, replace _ in species names, and save the df
    df_processed = df_processed.rename(columns={"qseqid": "ID"}).fillna(
        "Taxonomy unreliable"
    )
    df_processed["species"] = df_processed["species"].str.replace("_", " ")
    return df_processed


# Function to determine the cutoff rank of the processed taxonomy df
def determine_cutoff_rank(id_value):
    if id_value == "No match":
        return "No match"
    tax = "phylum"
    if id_value >= 75:
        tax = "class"
    if id_value >= 80:
        tax = "order"
    if id_value >= 85:
        tax = "family"
    if id_value >= 90:
        tax = "genus"
    if id_value >= 95:
        tax = "species"
    if id_value >= 98:
        tax = "none"
    return tax


# Define arguments
parser = argparse.ArgumentParser(
    description="A submodule for the abscale wrapper to run BLAST on generated FASTA files and filter BLAST hits."
)
parser.add_argument(
    "--fastafile",
    help="Input FASTA file.",
    required=True,
    metavar="file.fasta",
)
parser.add_argument(
    "--database",
    help="BLAST database.",
    metavar="/PATH/TO/DATABASE",
    required=True,
)
parser.add_argument(
    "--database_format",
    help="Format of the database. Currently available formats are: midori2, pr2.",
    choices=["midori2", "pr2"],
    required=True,
)
parser.add_argument(
    "--evalue",
    help="E-value for BLAST (default: 1e-05).",
    metavar="1e[exponent]",
    default="1e-05",
)
parser.add_argument(
    "--filter_mode",
    choices=["soft", "strict"],
    help="""Filter mode.

        soft:
        Keeps the best hit (highest bitscore) for each sequence. If multiple hits have the same highest bitscore, an LCA approach is applied (assigns the taxonomy to each sequence based on all taxonomic ranks that are identical in the remaining hits of each sequence).

        strict:
        Performs 3 steps:
        (1) bitscore filtering - keeps all hits with a bitscore >= --bitscore_treshold, an alignment length >= --length, and within --bitscore_percentage of the best bitscore per sequence.
        (2) similarity cutoff - only keeps the taxonomy of hits up to a certain rank, depending on the hits' blast percentage identity and cutoff values given in argument --cutoff_pidents.
        (3) LCA approach - assigns the taxonomy to each sequence based on all taxonomic ranks that are identical in the remaining hits of each sequence.
        """,
    default="strict",
)
parser.add_argument(
    "--bitscore_percentage",
    metavar="%",
    default=2.0,
    type=float,
    help=(
        """Percentage threshold (in %%) for bitscore filter when choosing
        filter_mode option "strict" (default=2.0)."""
    ),
)
parser.add_argument(
    "--alignment_length",
    metavar="NNN",
    default=100,
    type=int,
    help=(
        "Alignment length threshold to perform bitscore filtering on when "
        'choosing filter_mode option "strict" (default=100).'
    ),
)
parser.add_argument(
    "--bitscore_threshold",
    metavar="NNN",
    default=150,
    type=int,
    help=(
        """Bitscore threshold to perform bitscore filtering on when choosing
        filter_mode option "strict" (default=150)."""
    ),
)
parser.add_argument(
    "--cutoff_pidents",
    metavar="N,N,N,N,N,N",
    default=[98, 95, 90, 85, 80, 75],
    type=list_of_integers,
    help=(
        """Similarity cutoff per hit based on BLAST pident values when choosing
        filter_mode option "strict". cutoff pident values have to be divided by commas
        without spaces, in the order species, genus, family, order,
        class, phylum. Domain is always retained. Taxonomy is only kept for a rank if the BLAST hit's
        pident is >= the respective cutoff (default=98,95,90,85,80,75)."""
    ),
)
parser.add_argument(
    "--outfile",
    required=True,
    metavar="outfile.csv",
    help="Name of output file in .csv format.",
)
parser.add_argument(
    "--cores",
    metavar="N",
    default="2",
    help="Number of cores to use (default: 2).",
)

# Set arguments
args = parser.parse_args()
bitscore_percentage = args.bitscore_percentage / 100

# Start of pipeline
time_print(f"Running BLAST on {args.fastafile} with database {args.database}...")

# Run BLAST on FASTA file
subprocess.run(
    [
        "blastn",
        "-query",
        args.fastafile,
        "-db",
        args.database,
        "-out",
        "apscale_wrapper_blast_output.tsv",
        "-outfmt",
        "6 qseqid sseqid pident length bitscore",
        "-evalue",
        args.evalue,
        "-num_threads",
        args.cores,
    ]
)

# Load in and format BLAST output
time_print("Formatting BLAST output...")
df = pd.read_table(
    "apscale_wrapper_blast_output.tsv",
    header=None,
    names=["qseqid", "sseqid", "pident", "length", "bitscore"],
    dtype={"qseqid": str, "bitscore": float, "pident": float, "length": int},
)
# Clean up
os.remove("apscale_wrapper_blast_output.tsv")
if args.database_format == "midori2":
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    # Remove any of "phylum", "class", "order", "family", "genus", and "species" followed by _ as well as _ followed by a number and all the extra information before the domain
    df["taxonomy"] = (
        df["sseqid"]
        .str.replace(r"(phylum|class|order|family|genus|species)_", "", regex=True)
        .str.replace(r"_\d+", "", regex=True)
        .str.replace(r".*root;", "", regex=True)
    )
    # Split the taxonomy column by semicolon and expand into new columns
    df[ranks] = df["taxonomy"].str.split(";", expand=True)
    # Only keep desired columns and ranks and fill missing values with "NA"
    df = df.drop(["sseqid", "taxonomy"], axis=1).fillna("NA")

# TO BE ADDED AT SOME POINT
# elif args.database_format == "ncbi_nt":
#     # Drop rows containing "Unknown" = taxid could not be translated
#     df = df[~df[ranks].apply(lambda row: row.str.contains("Unknown")).any(axis=1)]

elif args.database_format == "pr2":
    ranks_all = [
        "domain",
        "supergroup",
        "phylum",
        "subdivision",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    # Split the taxonomy column by semicolon and expand into new columns
    df[ranks_all] = df["sseqid"].str.rstrip(";").str.split(";", expand=True)
    # Only keep desired columns and ranks and fill missing values with "NA"
    df = df.drop(["sseqid", "supergroup", "subdivision"], axis=1).fillna("NA")

    # Replace unknown taxonomy in the PR2 database
    df[ranks] = df[ranks].apply(replace_pr2_tax)


if args.filter_mode == "soft":
    time_print(
        "Grouping IDs and filtering hits based on the highest bitscore of each ID..."
    )
    idx = df.groupby(["qseqid"])["bitscore"].transform(max) == df["bitscore"]
    df = df[idx]

elif args.filter_mode == "strict":
    time_print("Filtering hits based on bitscore and length...")
    df.loc[
        (df["length"] < args.alignment_length)
        & (df["bitscore"] < args.bitscore_threshold),
        ranks,
    ] = "NA"

    time_print(
        "Grouping IDs and filtering hits based on bitscore cutoff for each ID..."
    )
    idx = (
        df.groupby(["qseqid"])["bitscore"].transform(
            bitscore_cutoff, percentage=bitscore_percentage
        )
        == df["bitscore"]
    )
    df = df[idx]

    # Copy df to keep version without cutoffs
    df_no_cutoffs = copy.deepcopy(df)

    time_print("Applying similarity cutoffs and LCA filter...")
    # Process the original df
    df.loc[df["pident"] < args.cutoff_pidents[0], "species"] = "NA"
    df.loc[df["pident"] < args.cutoff_pidents[1], "genus"] = "NA"
    df.loc[df["pident"] < args.cutoff_pidents[2], "family"] = "NA"
    df.loc[
        df["pident"] < args.cutoff_pidents[3], ["order", "suborder", "infraorder"]
    ] = "NA"
    df.loc[df["pident"] < args.cutoff_pidents[4], ["class", "subclass"]] = "NA"
    df.loc[df["pident"] < args.cutoff_pidents[5], ["phylum", "subphylum"]] = "NA"

    # Process the copied df so that a column for cutoff ranks is used instead of the actual cutoff
    df_no_cutoffs = post_processing(df_no_cutoffs)
    df_no_cutoffs["cutoff_rank"] = df_no_cutoffs["percentage_similarity"].apply(
        determine_cutoff_rank
    )
    # Save copied df
    df_no_cutoffs.to_csv(add_suffix(args.outfile), index=False)

# Process and save df
df = post_processing(df)
df = df.drop("percentage_similarity", axis=1)
df.to_csv(args.outfile, index=False)

time_print("BLAST filtering done.")
