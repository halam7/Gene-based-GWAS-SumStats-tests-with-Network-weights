install.packages("data.table", type = "source")
library(data.table)
library(dplyr)
library(tidyverse)


bim <- fread(here("data", "processed" ,"1000G_EUR_PHASE3_PLINK_MERGED", "1000G_EUR_merged.bim"), header = FALSE)
setnames(bim, c("CHR", "SNP", "GD", "POS", "A1_bim", "A2_bim"))

process_sumstats <- function(sumstats_filename, bim, output_prefix) {
  library(data.table)
  library(here)
  
  # Build full path to input file
  sumstats_path <- here("data", "processed", "sumstats_zscores", sumstats_filename)
  sumstats <- fread(sumstats_path)
  
  # Check for "SNP" column
  if (!"SNP" %in% colnames(sumstats)) {
    stop("Error: Sumstats file must contain a 'SNP' column.")
  }
  
  # Merge with BIM file
  merged <- merge(sumstats, bim, by = "SNP")
  
  # Allele match or flip
  merged$allele_match <- with(merged,
                              (A1 == A1_bim & A2 == A2_bim) |
                                (A1 == A2_bim & A2 == A1_bim)
  )
  
  # Mismatches
  mismatches <- merged[!merged$allele_match, ]
  
  # Flip Z-scores where alleles are reversed
  merged$flipped <- with(merged, A1 == A2_bim & A2 == A1_bim)
  merged$Z[merged$flipped] <- -merged$Z[merged$flipped]
  
  # Summary output
  cat("Total SNPs checked:", nrow(merged), "\n")
  cat("Matching alleles:", sum(merged$allele_match), "\n")
  cat("Mismatches:", nrow(mismatches), "\n")
  cat("Flips:", sum(merged$flipped), "\n")
  
  # Write cleaned summary stats to processed folder
  output_path <- here("data", "processed", "sumstats_cleaned", paste0(output_prefix, "_sumstats_clean.txt"))
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(merged, output_path, sep = "\t")
}

ADHD <- process_sumstats("ADHD_sumstats_Z.txt", bim, output_prefix = "ADHD")
AN <- process_sumstats("AN_sumstats_Z.txt", bim, output_prefix = "AN")
ASD <- process_sumstats("ASD_sumstats_Z.txt", bim, output_prefix = "ASD")
AXD <- process_sumstats("AXD_sumstats_Z.txt", bim,output_prefix = "AXD")
BD <- process_sumstats("BD_sumstats_Z.txt", bim, output_prefix = "BD")
CP <- process_sumstats("CP_sumstats_Z.txt", bim, output_prefix = "CP")
EA <- process_sumstats("EA_sumstats_Z.txt", bim, output_prefix = "EA")
MDD <- process_sumstats("MDD_sumstats_Z.txt", bim, output_prefix = "MDD")
NSM <- process_sumstats("NSM_sumstats_Z.txt", bim, output_prefix = "NSM")
OCD <- process_sumstats("OCD_sumstats_Z.txt", bim, output_prefix = "OCD")
SCZ <- process_sumstats("SCZ_sumstats_Z.txt", bim, output_prefix = "SCZ")
SMK <- process_sumstats("SmkInit_Z.txt", bim, output_prefix = "SMK")




