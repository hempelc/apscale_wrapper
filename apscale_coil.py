# !!! NOT FUNCTIONAL !!!

import pandas as pd
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects


# Activate conversion between R dataframes and pandas dataframes
pandas2ri.activate()

# Define a list of required libraries
required_libraries = ["seqinr", "readxl", "openxlsx", "dplyr", "coil", "tools", "rlang"]

# Install required libraries if not already installed
# Import R utils to install other packages
utils = rpackages.importr("utils")
# Set mirror
utils.chooseCRANmirror(ind=1)
for lib in required_libraries:
    if not rpackages.isinstalled(lib):
        print(f"The required R package {lib} is not installed. Installing now. This can take a while...")
        utils.install_packages(StrVector([lib]))

# Load required libraries
for lib in required_libraries:
    r.library(lib)

# Options
## Give path to ESV or OTU table with taxonomy
file = '/Users/simplexdna/Desktop/spongebob-coi_apscale_ESV_table_filtered_microdecon-filtered_with_taxonomy.csv'
## If you want to visually inspect sequences flagged by coil, set this to FALSE
# Otherwise, if TRUE, flagged seqs will be dropped automatically
auto_drop = True
# ESV or OTU?
unit="ESV"

# Import data
df = pd.read_csv(file)

# Function to identify the genetic code of the lowest rank that is recognized by coil's function which_translate_table
def genetic_code_lowest_rank(row):
    exceptions = ["Taxonomy unreliable - multiple matching taxa",
        "Taxonomy unreliable - percentage similarity threshold for rank not met",
        "Taxonomy unreliable - bitscore and alignment length threshold not met",
        "Taxonomy unreliable",
        "No match in database",
        "Unknown in PR2 database",
        "Unknown in BOLD database"]
    which_trans_table = robjects.r["which_trans_table"]
    result = "Not identifiable"
    col_names = reversed(row.index)
    for col in col_names:
        value = row[col]
        if pd.isna(value):
            continue  # Skip processing if value is NA
        if value in exceptions:
            continue  # Skip processing if value is NA
        try:
            result = which_trans_table("Chromadoridae")
            # SCRIPT BREAKS HERE BECAUSE WARNINGS FROM which_trans_table CAN'T BE CAUGHT
            if not warning_occurred:
            # Break out of the loop if processing succeeds
                break
        except Exception as e:
            # Continue to the next column if an error occurs
            print(f"Error for taxon {value}: {str(e)}")
            continue
    return result

# Apply the process_row function to each allowed rank of the dataframe
df['genetic_code'] = df[["phylum", "class", "order", "family"]].apply(lambda row: genetic_code_lowest_rank(row), axis=1)

# Create flag if coil can't process sequence (= not reliable)
# (likely because it's no animal sequenced due to no genetic code being found)
df['not_reliable'] = df['genetic_code'] == "Not identifiable"

# Mutate unidentifiable genetic codes 
df['genetic_code'] = df['genetic_code'].apply(lambda x: 0 if x == "Not identifiable" else x)

##### Formatting done

##### Start coil

# For primers mlCO1intF and jgHCO2198R, start position is 343/115 (nt/aa) and end is 657/218 (nt/aa)
# See coil vignette for more details
nt_start = 346
nt_end = 657
aa_start = 116
aa_end = 218

#Subsetting the PHHM to target length
coil = importr("coil")
meta_nt_phmm = coil.subsetPHMM(nt_coi_PHMM, start = nt_start, end = nt_end)
meta_aa_phmm = coil.subsetPHMM(aa_coi_PHMM, start = aa_start, end = aa_end)

# TO DO: CUT DOWN DF TO SAVE SPACE
df_reliable = df[~df['not_reliable']]
df_unreliable = df[df['not_reliable']]

# Running coil
# Note: I turned on triple translate just in case the PHMM subsetting is not 100% accurate
# Runs longer but is safer
coil_result_df = r.flatten_coi5p(
  lapply(df_reliable['ID'], function(i):
    r.coi5p_pipe(df_reliable.loc[i, 'Seq'], 
                 name = df_reliable.loc[i, 'ID'], 
                 trans_table = df_reliable.loc[i, 'genetic_code'],
                 nt_PHMM = meta_nt_phmm,
                 aa_PHMM = meta_aa_phmm,
                 triple_translate=True)
  )
)

# Concatenating df and coil output
df_with_coil = pd.concat([df_reliable, pandas2ri.ri2py(coil_result_df)], axis=1)

# Create coil flag columns that flags sequences that contain either an indel or stop codon
df_with_coil['coil_flag'] = (df_with_coil['indel_likely'] | df_with_coil['stop_codons']).map(lambda x: "CONTAINS INDEL OR STOP CODON" if x else None)

# Concatenate dfs
df_unreliable = df_unreliable.rename(columns=lambda x: x if x in df_with_coil.columns else None)
concatenated_df = pd.concat([df_with_coil, df_unreliable])

# Sort df by ID names
concatenated_df = concatenated_df.astype({'ID': str}).sort_values(by='ID')

if auto_drop:
    # Set outfile name
    outfile = file.replace(f".{fileformat}", "_coil_filtered.xlsx")
    # Drop rows that are flagged
    num_dropped = concatenated_df['coil_flag'].notna().sum()
    print(f"Dropping {num_dropped} flagged sequences.")
    concatenated_df = concatenated_df[concatenated_df['coil_flag'].isna()]
    # Drop coil-related columns
    concatenated_df = concatenated_df.drop(columns=['genetic_code', 'not_reliable', 'indel_likely', 'stop_codons', 'coil_flag'])
    concatenated_df.to_excel(outfile, index=False)
else:
    # Set outfile name
    outfile = file.replace(f".{fileformat}", "_with_coilflag.xlsx")
    # Cut down df
    concatenated_df = concatenated_df.drop(columns=['genetic_code', 'not_reliable', 'indel_likely', 'stop_codons'])
    concatenated_df.to_excel(outfile, index=False)
