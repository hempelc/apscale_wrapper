# Script to identify indels and stop codons in COI sequences using coil
# Requires an ESV or OTU table with taxonomy results
# Note: the script leaves non-animal sequences untouched, only animal sequences are checked

# Define a list of required libraries
required_libraries <- c("seqinr", "readxl", "openxlsx", "dplyr", "coil", "tools")

# Check if libraries are installed, if not, install them
for (lib in required_libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    message(paste("The required R package", lib, "is not installed. Installing now. This can take a while..."))
    chooseCRANmirror(ind=1)
    install.packages(lib, dependencies = TRUE)
  }
}

# Load required libraries
library(seqinr)
library(readxl)
library(openxlsx)
library(dplyr)
library(coil)
library(tools)

# Get command-line arguments (order: file, unit)
args <- commandArgs(trailingOnly = TRUE)

# Options
## Set path to ESV or OTU table with taxonomy
file = as.character(args[1])
# Set unit (ESV or OTU)
unit = as.character(args[2])

# Import data
df <- read.csv(file)


# Function to identify the genetic code of the lowest rank that is recognized by coil's function which_translate_table
genetic_code_lowest_rank <- function(row) {
  exceptions = c("Taxonomy unreliable - multiple matching taxa",
                "Taxonomy unreliable - percentage similarity threshold for rank not met",
                "Taxonomy unreliable - bitscore and alignment length threshold not met",
                "Taxonomy unreliable",
                "No match in database",
                "Unknown in PR2 database",
                "Unknown in BOLD database")
  result <- "Not identifiable"
  for (col in rev(names(row))) {
    value <- row[[col]]
    if (is.na(value) || value %in% exceptions) {
      next  # Skip processing if value is in exceptions or NA
    }
    tryCatch({
      result <- which_trans_table(value)
      break  # Break out of the loop if processing succeeds
    }, error = function(e) {
      # Continue to the next column if an error occurs
      # Deactivating warning and error messages for this function as they could clutter the stdout
      #warning(paste("Error for taxon", value, ":", e$message))
    }, warning = function(w) {
      # Continue to the next column if a warning occurs
      #message(paste("Warning for taxon", value, ":", w$message))
    })
  }
  return(result)
}

# Apply the process_row function to each allowed rank of the dataframe
df$genetic_code <- apply(df[c("phylum", "class", "order", "family")], 1, genetic_code_lowest_rank)

# Create flag if coil can't process sequence (= not reliable)
# (likely because it's no animal sequenced due to no genetic code being found)
df <- df %>%
  mutate(not_reliable = genetic_code == "Not identifiable")

# Mutate unidentifiable genetic codes 
df <- df %>%
  mutate(genetic_code = ifelse(genetic_code == "Not identifiable", 0, genetic_code))

# Genetic code needs to be in format integer
df$genetic_code <- as.integer(df$genetic_code)

##### Formatting done

##### Start coil

# For primers mlCO1intF and jgHCO2198R, start position is 343/115 (nt/aa) and end is 657/218 (nt/aa)
# See coil vignette for more details
nt_start = 346
nt_end = 657
aa_start = 116
aa_end = 218

#Subsetting the PHHM to target length
meta_nt_phmm = subsetPHMM(nt_coi_PHMM, start = nt_start, end = nt_end)
meta_aa_phmm = subsetPHMM(aa_coi_PHMM, start = aa_start, end = aa_end)

# TO DO: CUT DOWN DF TO SAVE SPACE
df_reliable = df[!df$not_reliable, ]
df_unreliable = df[df$not_reliable, ]

# Running coil
# Note: I turned on triple translate just in case the PHMM subsetting is not 100% accurate
# Runs longer but is safer
coil_result_df = flatten_coi5p(
  lapply(1:length(df_reliable$ID), function(i){
    coi5p_pipe(df_reliable$Seq[i], 
               name = df_reliable$ID[i], 
               trans_table = df_reliable$genetic_code[i],
               nt_PHMM = meta_nt_phmm,
               aa_PHMM = meta_aa_phmm,
               triple_translate=TRUE)
  })
)

# Concatenating df and coil output
df_with_coil = cbind(df_reliable, coil_result_df)

# Create coil flag columns that flags sequences that contain either an indel or stop codon
df_with_coil$coil_flag <- df_with_coil$indel_likely | df_with_coil$stop_codons
df_with_coil$coil_flag <- ifelse(df_with_coil$coil_flag == TRUE, "CONTAINS INDEL OR STOP CODON", NA)

# Concatenate dfs
df_unreliable <- df_unreliable %>%
  rename_all(~ ifelse(. %in% colnames(df_with_coil), ., NA))
concatenated_df <- bind_rows(df_with_coil, df_unreliable)

# Sort df by ID names
concatenated_df <- concatenated_df %>%
  arrange(as.numeric(gsub(paste0(unit, "_"), "", ID)))

# Set outfile name
outfile <- sub(paste0("\\.csv$"), "_coil_filtered.csv", file)
# Drop rows that are flagged
num_dropped = sum(!is.na(concatenated_df$coil_flag))
warning("R package coil: ", num_dropped, " sequences have been identified as NUMTs and were removed.")
concatenated_df = subset(concatenated_df, is.na(coil_flag))
# Drop coil-related columns
concatenated_df = select(concatenated_df, -genetic_code:-coil_flag)
write.csv(concatenated_df, file = outfile, row.names = FALSE)
