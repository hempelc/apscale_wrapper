#!/usr/bin/env python3

"""
A submodule for the abscale wrapper to run SINTAX on generated FASTA files.

By Chris Hempel (christopher.hempel@kaust.edu.sa) on 29 Apr 2024
"""

import datetime
import argparse
import warnings
import subprocess
import os
import re
import pandas as pd

# Define that warnings are not printed to console
warnings.filterwarnings("ignore")


# Funtion to print datetime and text
def time_print(text):
    print(datetime.datetime.now().strftime("%H:%M:%S"), ": ", text, sep="")


# Function to replace unknown taxonomy in the PR2 database
def replace_pr2_tax(taxon_series):
    return taxon_series.apply(
        lambda taxon: (
            "Unknown in PR2 database"
            if taxon.endswith("X") or taxon.endswith("sp.") or taxon.endswith("_sp")
            else taxon
        )
    )


# Function to add suffix to filename despite format
def add_suffix(filename, suffix="_no_cutoff"):
    root, extension = os.path.splitext(filename)
    return f"{root}{suffix}{extension}"


# Function to get lowest identified taxon and rank per row
def lowest_taxon_and_rank(row):
    exceptions = [
        "Taxonomy unreliable - confidence threshold not met",
        "Unknown in PR2 database",
        "Unknown in BOLD database",
        "Unknown in SILVA database",
        "Unknown in MIDORI2 database",
    ]
    lowest_taxon = "Taxonomy unreliable"
    lowest_rank = "Taxonomy unreliable"

    # Iterate over columns in reverse order
    for col in reversed(ranks):
        value = re.sub(r"\(\d\.\d{2}\)", "", row[col])
        if value == "No match in database":
            lowest_taxon = value
            lowest_rank = value
            break
        elif pd.notna(value) and value not in exceptions:
            lowest_taxon = value
            lowest_rank = col
            break  # Break on the first valid entry found

    return pd.Series([lowest_taxon, lowest_rank], index=["lowest_taxon", "lowest_rank"])


# Define function to process tax dfs
def post_processing(df):
    # Remove underscores
    df[ranks] = df[ranks].replace(r"_", " ", regex=True)
    # Make a df without confidence values (for df_no_cutoff)
    df_ranks = df[ranks].replace(r"\(\d\.\d{2}\)", "", regex=True)
    lowest_columns = df_ranks.apply(lowest_taxon_and_rank, axis=1)
    df = pd.concat([df, lowest_columns], axis=1)
    return df


# Function to replace unknown MIDORI2 ranks
def replace_if_match(taxon):
    if bool(re.search(r"(phylum|class|order|family|genus|species)_", taxon)):
        return "Unknown in MIDORI2 database"
    return taxon


def tax_formatting(df, tax_col, ranks):
    # Taxonomy formatting
    if args.database_format == "silva":
        # Annotate entries with no DB match
        df = df.fillna(",".join(["No match in database" for _ in range(len(ranks))]))
        # Split the tax_col column by comma and expand into new columns
        df[ranks] = df[tax_col].str.split(",", expand=True)
        # Only keep desired columns and ranks and fill missing values
        df = df.drop([tax_col], axis=1).fillna(
            "Taxonomy unreliable - confidence threshold not met"
        )
        # Replace taxa containing Not_available
        df[ranks] = df[ranks].replace("Not_available", "Unknown in SILVA database")

    elif args.database_format == "midori2":
        # Annotate entries with no DB match
        df = df.fillna(",".join(["No match in database" for _ in range(len(ranks))]))
        # Split the tax_col column by comma and expand into new columns
        df[ranks] = df[tax_col].str.split(",", expand=True)
        # Fill missing values with threshold note
        df[ranks] = df[ranks].fillna(
            "Taxonomy unreliable - confidence threshold not met"
        )
        # If any taxon starts with "phylum", "class", "order", "family", "genus", or "species" followed by _, replace the taxa with "Unknown in MIDORI2 database"
        for rank in ranks:
            df[rank] = df[rank].apply(replace_if_match)
        # Only keep desired columns and ranks
        df = df.drop([tax_col], axis=1)

    elif args.database_format == "bold":
        # Annotate entries with no DB match
        df = df.fillna(",".join(["No match in database" for _ in range(len(ranks))]))
        # Split the tax_col column by comma and expand into new columns
        df[ranks] = (
            df[tax_col]
            .str.split(",", expand=True)
            .fillna("Taxonomy unreliable - confidence threshold not met")
        )
        # Replace species names containing " sp. "" or ending with "sp"
        mask = (
            df["species"].str.endswith("_sp")
            | df["species"].str.contains("_sp\._", regex=True)
            | df["species"].str.contains("_cf\._", regex=True)
        )
        df.loc[mask, "species"] = "Unknown in BOLD database"
        # Replace taxa containing incertae_sedis
        df[ranks] = df[ranks].applymap(
            lambda x: "Unknown in BOLD database" if "incertae_sedis" in x else x
        )
        # Replace 'Unknown_in_BOLD_database' by 'Unknown in BOLD database' in the entire DataFrame
        df = df.replace("Unknown_in_BOLD_database", "Unknown in BOLD database")
        # Only keep desired columns and ranks and fill missing values with "NA"
        df = df.drop([tax_col], axis=1)

    # TO BE ADDED AT SOME POINT
    # elif args.database_format == "ncbi_nt":
    #     # Drop rows containing "Unknown" = taxid could not be translated
    #     df = df[~df[ranks].apply(lambda row: row.str.contains("Unknown")).any(axis=1)]

    elif args.database_format == "pr2":
        ranks_pr2 = [
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
        # Annotate entries with no DB match
        df = df.fillna(
            ",".join(["No match in database" for _ in range(len(ranks_pr2))])
        )
        # Split the taxonomy column by comma and expand into new columns
        df[ranks_pr2] = (
            df[tax_col]
            .str.replace(":nucl", "")
            .str.replace(":plas", "")
            .str.replace(":apic", "")
            .str.replace(":chrom", "")
            .str.replace(":mito", "")
            .str.split(",", expand=True)
        )
        # Only keep desired columns and ranks and fill missing values with "NA"
        df = df.drop([tax_col, "supergroup", "subdivision"], axis=1).fillna(
            "Taxonomy unreliable - confidence threshold not met"
        )
        # Replace unknown taxonomy in the PR2 database
        df[ranks] = df[ranks].apply(replace_pr2_tax)
    return df


# Define arguments
parser = argparse.ArgumentParser(
    description="A submodule for the abscale wrapper to run SINTAX on generated FASTA files."
)
parser.add_argument(
    "--fastafile",
    help="Input FASTA file.",
    required=True,
    metavar="file.fasta",
)
parser.add_argument(
    "--sintax_database",
    help="SINTAX database.",
    metavar="/PATH/TO/DATABASE",
    required=True,
)
parser.add_argument(
    "--database_format",
    help="Format of the database. Note: the SILVA and BOLD databases have to have a specific format.",
    choices=["midori2", "pr2", "silva", "bold"],
    required=True,
)
parser.add_argument(
    "--sintax_confidence_cutoff",
    metavar="N.N",
    default="0.75",
    help=("""Confidence value cutoff level, ranging from 0-1 (default=0.75)."""),
)
parser.add_argument(
    "--outfile",
    required=True,
    metavar="OUTFILENAME.csv",
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
if args.database_format == "bold":
    ranks = ["phylum", "class", "order", "family", "genus", "species"]
else:
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# Start of pipeline
time_print(
    f"Running SINTAX on {args.fastafile} with database {args.sintax_database}..."
)

# Run SINTAX on FASTA file
sintaxout = os.path.join(
    os.path.dirname(args.fastafile), "apscale_wrapper_sintax_output.tsv"
)
subprocess.run(
    [
        "vsearch",
        "--sintax",
        args.fastafile,
        "--db",
        args.sintax_database,
        "--tabbedout",
        sintaxout,
        "--sintax_cutoff",
        args.sintax_confidence_cutoff,
        "--threads",
        args.cores,
    ]
)

# Load in SINTAX output
time_print("Formatting SINTAX output...")
df = pd.read_table(
    sintaxout,
    header=None,
    delim_whitespace=True,  # Ensure consistent parsing
    names=["ID", "tax_with_scores", "spacer", "tax"],
    usecols=["ID", "tax_with_scores", "tax"],
)

# Revert SINTAX-specific format
df["tax"] = (
    df["tax"]
    .str.replace(r"(d|k|p|c|o|f|g|s):", "", regex=True)
    .str.replace(r"_\d+", "", regex=True)
)
df["tax_with_scores"] = (
    df["tax_with_scores"]
    .str.replace(r"(d|k|p|c|o|f|g|s):", "", regex=True)
    .str.replace(r"_\d+", "", regex=True)
)

# Sort by ID number
df["OTU_ESV_number"] = df["ID"].str.extract(r"_(\d+)", expand=False).astype(int)
df = df.sort_values(by="OTU_ESV_number")
df = df.drop(columns=["OTU_ESV_number"])

# Separate both dfs
df_no_cutoffs = df[["ID", "tax_with_scores"]]
df_with_cutoffs = df[["ID", "tax"]]

# Format taxonomy
df_no_cutoffs = tax_formatting(df_no_cutoffs, "tax_with_scores", ranks)
df_with_cutoffs = tax_formatting(df_with_cutoffs, "tax", ranks)

# Process and save df
df_no_cutoffs = post_processing(df_no_cutoffs)
df_with_cutoffs = post_processing(df_with_cutoffs)

# Save
df_no_cutoffs.to_csv(add_suffix(args.outfile), index=False)
df_with_cutoffs.to_csv(args.outfile, index=False)

time_print("SINTAX formatting done.")
