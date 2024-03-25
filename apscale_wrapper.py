#!/usr/bin/env python3

"""
A wrapper to run apscale on forward and reverse reads and to generate 
processing QC graphs.

Requires the submodules apscale_processing_graphs.py, apscale_blast.py, and apscale_remove_negatives.py in PATH.

By Chris Hempel (christopher.hempel@kaust.edu.sa) on 15 Aug 2023
"""

import datetime
import pandas as pd
import argparse
import warnings
import subprocess
from pathlib import Path
import os
import sys


# Define that warnings are not printed to console
warnings.filterwarnings("ignore")


# Funtion to print datetime and text
def time_print(text):
    timetext = datetime.datetime.now().strftime("%H:%M:%S") + ": " + text
    print(timetext)
    log.write(timetext + '\n')  # Write to the log file
    log.flush()  # Flush the buffer to ensure immediate writing to the file


# Define a custom validation function to enforce the requirement for --remove_negative_controls and --run_blast
def validate_args(args):
    if args.remove_negative_controls == "True" and not args.negative_controls:
        parser.error(
            "--negative_controls is required when --remove_negative_controls=True."
        )
    if args.run_blast == "True" and (not args.database or not args.database_format):
        parser.error(
            "--database and --database_format are required when --run_blast=True."
        )


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


# Function to BLAST FASTA files using the apscale blast submodule
def blasting(fastafile, outfile, **kwargs):
    proc=subprocess.run(
        [
            "apscale_blast.py",
            "--fastafile",
            fastafile,
            "--database",
            args.database,
            "--database_format",
            args.database_format,
            "--evalue",
            args.evalue,
            "--filter_mode",
            args.filter_mode,
            "--bitscore_percentage",
            args.bitscore_percentage,
            "--alignment_length",
            args.alignment_length,
            "--bitscore_threshold",
            args.bitscore_threshold,
            "--cutoff_pidents",
            args.cutoff_pidents,
            "--outfile",
            outfile,
            "--cores",
            args.cores,
        ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    for line in proc.stdout:
        sys.stdout.write(str(line))
        log.write(str(line))
    
# Define arguments
parser = argparse.ArgumentParser(
    description="""A wrapper to run apscale on forward and reverse 
	reads, filter Negative control reads, BLAST generated sequences, and generate processing QC graphs.
	Requires the submodules apscale_processing_graphs.py, apscale_blast.py, and apscale_remove_negatives.py
    in PATH.""",
)
parser.add_argument(
    "--sequence_dir",
    metavar="PATH/TO/SEQUENCE_DIR",
    help="""Path to directory containing demultiplexed forward and reverse sequences in .fastq.gz format.
    Note: file names must end with R1 and R2 to identify read pairs.""",
    required=True,
)
parser.add_argument(
    "--project_name",
    help="Name of the apscale project to be generated.",
    required=True,
)
parser.add_argument(
    "--forward_primer",
    metavar="PRIMER_SEQUENCE",
    help="Forward primer sequence to trim.",
    required=True,
)
parser.add_argument(
    "--reverse_primer",
    metavar="PRIMER_SEQUENCE",
    help="Reverse primer sequence to trim.",
    required=True,
)
parser.add_argument(
    "--min_length",
    metavar="NNN",
    type=int,
    help="Minimum limit of expected amplicon length (used for length filtering).",
    required=True,
)
parser.add_argument(
    "--max_length",
    metavar="NNN",
    type=int,
    help="Maximum limit of expected amplicon length (used for length filtering).",
    required=True,
)
parser.add_argument(
    "--otu_perc",
    metavar="NN",
    default=97,
    type=int,
    help="OTU identify treshold for clustering with vsearch. Only used when --clusteringtool=vsearch (default=97).",
)
parser.add_argument(
    "--coi",
    choices=["False", "True"],
    default="False",
    help="""Are you processing COI data? If yes, the fact that COI is a coding gene can be used to
    improve clustering and denoising (default=False).""",
)
parser.add_argument(
    "--swarm_distance",
    default=1,
    metavar="N",
    type=int,
    help="Distance used for swarm. Overwritten to 13 if --coi=True. Only used when --clusteringtool=swarm (default: 1).",
)
parser.add_argument(
    "--clusteringtool",
    default="vsearch",
    choices=["vsearch", "swarm"],
    help="Tool used for OTU clustering (default=vsearch).",
)
parser.add_argument(
    "--prior_denoising",
    choices=["False", "True"],
    default="False",
    help="Set to True if you want to denoise reads prior to OTU clustering (default: False).",
)
parser.add_argument(
    "--maxEE",
    metavar="N",
    default=2,
    type=int,
    help="maxEE (maximum estimated error) value used for quality filtering (default: 2).",
)
parser.add_argument(
    "--graph_format",
    help="Graph format, either HTML, png, svg. HTML is recommended because it creates interactive plots (default: html).",
    default="html",
    choices=["png", "svg", "html"],
)
parser.add_argument(
    "--scaling_factor",
    help="Scaling factor for graph width. Manual trial and error in 0.2 increments might be required (default: 1.0).",
    default=1.0,
    metavar="N.N",
    type=float,
)
parser.add_argument(
    "--remove_negative_controls",
    help="Do you want to remove reads in negative controls from the samples using the R package microDecon? (default: False)",
    default="False",
    choices=["True", "False"],
)
parser.add_argument(
    "--negative_controls",
    help="Required if --remove_negative_controls=True. List the names of all negative controls (without _R1/2.fq.gz), separated by commas without spaces.",
    metavar="control1,control2,control3",
)
parser.add_argument(
    "--run_blast",
    help="Do you want to run BLAST on the generated ESVs and OTUs? (default: False)",
    default="False",
    choices=["True", "False"],
)
parser.add_argument(
    "--database",
    help="Required if --run_blast=True. BLAST database.",
    metavar="/PATH/TO/DATABASE",
)
parser.add_argument(
    "--database_format",
    help="Required if --run_blast=True. Format of the database. Currently available formats are: midori2, pr2, silva. Note: the SILVA database has to have a specific format.",
    choices=["midori2", "pr2", "silva"],
)
parser.add_argument(
    "--evalue",
    help="Used if --run_blast=True. E-value for BLAST (default: 1e-05).",
    metavar="1e[exponent]",
    default="1e-05",
)
parser.add_argument(
    "--filter_mode",
    choices=["soft", "strict"],
    help="""Used if --run_blast=True. Filter mode.

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
    default="2.0",
    help=(
        """Used if --run_blast=True. Percentage threshold (in %%) for bitscore filter when choosing
        filter_mode option "strict" (default=2.0)."""
    ),
)
parser.add_argument(
    "--alignment_length",
    metavar="NNN",
    default="100",
    help=(
        "Used if --run_blast=True. Alignment length threshold to perform bitscore filtering on when "
        'choosing filter_mode option "strict" (default=100).'
    ),
)
parser.add_argument(
    "--bitscore_threshold",
    metavar="NNN",
    default="150",
    help=(
        """Used if --run_blast=True. Bitscore threshold to perform bitscore filtering on when choosing
        filter_mode option "strict" (default=150)."""
    ),
)
parser.add_argument(
    "--cutoff_pidents",
    metavar="N,N,N,N,N,N",
    default="98,95,90,85,80,75",
    help=(
        """Used if --run_blast=True. Similarity cutoff per hit based on BLAST pident values when choosing
        filter_mode option "strict". cutoff pident values have to be divided by commas
        without spaces, in the order species, genus, family, order,
        class, phylum. Domain is always retained. Taxonomy is only kept for a rank if the BLAST hit's
        pident is >= the respective cutoff (default=98,95,90,85,80,75)."""
    ),
)
parser.add_argument(
    "--cores",
    metavar="N",
    default="2",
    help="Number of cores to use (default: 2).",
)

# Register the custom validation function to be called after parsing arguments
parser.set_defaults(func=validate_args)

# Parse arguments
args = parser.parse_args()

# Call the custom validation function to check the requirements
args.func(args)

# Overwrite the parameter for --swarm_distance if --coi=="True"
if args.coi == "True":
    args.swarm_distance = 13
    
# Log options
## Define the log file name
log_file = f"{args.project_name}_apscale_wrapper.log"
# Open the log file in append mode
log = open(log_file, 'a')


### Start of pipeline
time_print("Starting apscale wrapper.")

# Create an apscale directory using bash
time_print("Creating apscale directory...")
proc=subprocess.run(["apscale", "--create_project", args.project_name],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in proc.stdout:
    sys.stdout.write(str(line))
    log.write(str(line))

# Create an empty Project_report.xlsx file
## Create an ExcelWriter object using the openpyxl engine
excel_writer = pd.ExcelWriter(
    os.path.join(f"{args.project_name}_apscale", "Project_report.xlsx"),
    engine="openpyxl",
)
## Write an empty DataFrame to the Excel file
pd.DataFrame().to_excel(excel_writer, sheet_name="Sheet1", index=False)
# Save the Excel file
excel_writer.book.save(
    os.path.join(f"{args.project_name}_apscale", "Project_report.xlsx")
)

# Generate symlinks to demultiplexed reads
time_print("Generating symlinks to demultiplexed reads...")
# sequence_files = os.path.join(os.path.realpath(args.sequence_dir), "*")
target_directory = os.path.join(
    f"{args.project_name}_apscale", "2_demultiplexing", "data"
)
# TO DO: involve os.path.join for sequence files to make it universal for systems
# subprocess.run(
#     ["ln", "-s", sequence_files, target_directory],
#     shell=True,
# )
proc=subprocess.run(
    f'ln -s "$(realpath "{args.sequence_dir}")"/* {target_directory}',
    shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
)
for line in proc.stdout:
    sys.stdout.write(str(line))
    log.write(str(line))

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
proc=subprocess.run(["apscale", "--run_apscale", f"{args.project_name}_apscale"],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
for line in proc.stdout:
    sys.stdout.write(str(line))
    log.write(str(line))

time_print("Apscale done.")

if args.remove_negative_controls == "True":
    microdecon_suffix = "_microdecon-filtered"
    proc=subprocess.run(
        [
            "apscale_remove_negatives.py",
            "--project_dir",
            f"{args.project_name}_apscale",
            "--negative_controls",
            f"{args.negative_controls}",
        ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    for line in proc.stdout:
        sys.stdout.write(str(line))
        log.write(str(line))
    
else:
    microdecon_suffix = ""

# Run BLAST if specified
if args.run_blast == "True":
    # Define file names
    fastafile_otus = os.path.join(
        f"{args.project_name}_apscale",
        "9_lulu_filtering",
        "otu_clustering",
        f"{args.project_name}_apscale_OTUs_filtered{microdecon_suffix}.fasta",
    )
    fastafile_esvs = os.path.join(
        f"{args.project_name}_apscale",
        "9_lulu_filtering",
        "denoising",
        f"{args.project_name}_apscale_ESVs_filtered{microdecon_suffix}.fasta",
    )
    blastoutFile_otus = os.path.join(
        f"{args.project_name}_apscale",
        "9_lulu_filtering",
        "otu_clustering",
        f"{args.project_name}_apscale_OTUs{microdecon_suffix}_blasted.csv",
    )
    blastoutFile_esvs = os.path.join(
        f"{args.project_name}_apscale",
        "9_lulu_filtering",
        "denoising",
        f"{args.project_name}_apscale_ESVs{microdecon_suffix}_blasted.csv",
    )
    otu_outfile = os.path.join(
        f"{args.project_name}_apscale",
        "9_lulu_filtering",
        "otu_clustering",
        f"{args.project_name}_apscale_OTU_table_filtered{microdecon_suffix}_with_taxonomy.csv",
    )
    esv_outfile = os.path.join(
        f"{args.project_name}_apscale",
        "9_lulu_filtering",
        "denoising",
        f"{args.project_name}_apscale_ESV_table_filtered{microdecon_suffix}_with_taxonomy.csv",
    )
    if args.filter_mode == "strict":
        blastoutFile_otus_noCutoff = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "otu_clustering",
            f"{args.project_name}_apscale_OTUs{microdecon_suffix}_blasted_no_cutoff.csv",
        )
        blastoutFile_esvs_noCutoff = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "denoising",
            f"{args.project_name}_apscale_ESVs{microdecon_suffix}_blasted_no_cutoff.csv",
        )
        otu_outfile_noCutoff = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "otu_clustering",
            f"{args.project_name}_apscale_OTU_table_filtered{microdecon_suffix}_with_taxonomy_no_cutoff.csv",
        )
        esv_outfile_noCutoff = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "denoising",
            f"{args.project_name}_apscale_ESV_table_filtered{microdecon_suffix}_with_taxonomy_no_cutoff.csv",
        )
    if args.remove_negative_controls == "True":
        otu_table_file = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "otu_clustering",
            f"{args.project_name}_apscale_OTU_table_filtered{microdecon_suffix}.csv",
        )
        esv_table_file = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "denoising",
            f"{args.project_name}_apscale_ESV_table_filtered{microdecon_suffix}.csv",
        )
    else:
        otu_table_file = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "otu_clustering",
            f"{args.project_name}_apscale_OTU_table_filtered{microdecon_suffix}.parquet.snappy",
        )
        esv_table_file = os.path.join(
            f"{args.project_name}_apscale",
            "9_lulu_filtering",
            "denoising",
            f"{args.project_name}_apscale_ESV_table_filtered{microdecon_suffix}.parquet.snappy",
        )

    # Run BLAST command
    blasting(
        fastafile=fastafile_otus,
        outfile=blastoutFile_otus,
        database=args.database,
        database_format=args.database_format,
        evalue=args.evalue,
        filter_mode=args.filter_mode,
        bitscore_percentage=args.bitscore_percentage,
        alignment_length=args.alignment_length,
        bitscore_threshold=args.bitscore_threshold,
        cutoff_pidents=args.cutoff_pidents,
        cores=args.cores,
    )
    blasting(
        fastafile=fastafile_esvs,
        outfile=blastoutFile_esvs,
        database=args.database,
        database_format=args.database_format,
        evalue=args.evalue,
        filter_mode=args.filter_mode,
        bitscore_percentage=args.bitscore_percentage,
        alignment_length=args.alignment_length,
        bitscore_threshold=args.bitscore_threshold,
        cutoff_pidents=args.cutoff_pidents,
        cores=args.cores,
    )

    time_print("OTUs and ESVs BLASTed. Merging taxonomy with OTU/ESV tables...")
    # Read in OTU and ESV tables
    if args.remove_negative_controls == "True":
        otu_table = pd.read_csv(otu_table_file)
        esv_table = pd.read_csv(esv_table_file)
    else:
        otu_table = pd.read_parquet(otu_table_file, engine="fastparquet")
        esv_table = pd.read_parquet(esv_table_file, engine="fastparquet")
    # Read in BLAST results
    blastout_otus = pd.read_csv(blastoutFile_otus)
    blastout_esvs = pd.read_csv(blastoutFile_esvs)
    # Clean up
    os.remove(blastoutFile_otus)
    os.remove(blastoutFile_esvs)
    # Merge tables
    otu_table_with_tax = pd.merge(
        otu_table,
        blastout_otus,
        on="ID",
        how="outer",
    ).fillna("No match in database")
    esv_table_with_tax = pd.merge(
        esv_table,
        blastout_esvs,
        on="ID",
        how="outer",
    ).fillna("No match in database")
    # Save
    otu_table_with_tax.to_csv(otu_outfile, index=False)
    esv_table_with_tax.to_csv(esv_outfile, index=False)
    if args.filter_mode == "strict":
        # Read in BLAST results with no cutoff
        blastout_otus_noCutoff = pd.read_csv(blastoutFile_otus_noCutoff)
        blastout_esvs_noCutoff = pd.read_csv(blastoutFile_esvs_noCutoff)
        # Clean up
        os.remove(blastoutFile_otus_noCutoff)
        os.remove(blastoutFile_esvs_noCutoff)
        # Merge tables
        otu_table_with_tax_noCutoff = pd.merge(
            otu_table,
            blastout_otus_noCutoff,
            on="ID",
            how="outer",
        ).fillna("No match in database")
        esv_table_with_tax_noCutoff = pd.merge(
            esv_table,
            blastout_esvs_noCutoff,
            on="ID",
            how="outer",
        ).fillna("No match in database")
        # Save
        otu_table_with_tax_noCutoff.to_csv(otu_outfile_noCutoff, index=False)
        esv_table_with_tax_noCutoff.to_csv(esv_outfile_noCutoff, index=False)

# Generate processing graphs using separate script
proc=subprocess.run(
    [
        "apscale_processing_graphs.py",
        "--project_dir",
        f"{args.project_name}_apscale",
        "--graph_format",
        f"{args.graph_format}",
        "--scaling_factor",
        f"{args.scaling_factor}",
        "--blast",
        f"{args.run_blast}",
        "--remove_negative_controls",
        args.remove_negative_controls,
    ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
)
for line in proc.stdout:
    sys.stdout.write(str(line))
    log.write(str(line))


# Generate detailed report
# TO DO

time_print("Apscale wrapper done.")

# Close the log file
log.close()