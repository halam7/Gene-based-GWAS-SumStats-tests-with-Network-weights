args <- commandArgs(trailingOnly = TRUE)
sumstats_path <- args[1]
chunk_id <- as.integer(args[2])  # new argument

.libPaths("/home/halam7/project/R/4.4")

library(mkatr)
library(ieugwasr)
library(genetics.binaRies)
library(here)
library(data.table)
library(R.utils)
library(callr)

source(here("scripts", "16_safesats_penalty_genes_parallel.R"))

setid <- fread(here("data", "processed", "setid", "SKAT_setID.txt"), header = FALSE)
colnames(setid) <- c("Gene", "SNP")
bfile <- here("data", "processed", "1000G_EUR_PHASE3_PLINK_MERGED", "1000G_EUR_merged")

# Load gene list for this chunk
gene_list <- readLines(here("data", "processed", "gene_chunks", paste0("chunk_", chunk_id, ".txt")))

run_sats_analysis(sumstats_path, setid, bfile, gene_list, chunk_id)

cat("Running chunk", chunk_id, "on trait file:", sumstats_path, "\n")