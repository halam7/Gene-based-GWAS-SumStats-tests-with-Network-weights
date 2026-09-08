## A function that takes in a sumstats_weights.text file and runs all the
## 5 weights and outputs the results.
library(here)
library(data.table)
#install.packages("R.utils")
#install.packages("Matrix")

setid <- fread(here("data", "processed", "setid", "SKAT_setID.txt"), header = FALSE)
colnames(setid) <- c("Gene", "SNP")

bfile <- here("data", "processed", "1000G_EUR_PHASE3_PLINK_MERGED", "1000G_EUR_merged")
sumstats_dir <- here("data", "processed", "sumstats_weights")

run_sats_analysis <- function(sumstats_path, setid, bfile, gene_list, chunk_id) {
  library(mkatr)
  library(data.table)
  library(ieugwasr)
  library(genetics.binaRies)
  library(here)
  library(callr)   # for process-level timeout
  library(Matrix)
  library(MASS)
  
  # Safe wrapper: runs sats in a clean R session with a hard timeout
  safe_sats <- function(Z, R, W, timeout = 120, label = "sats") {
    tryCatch({
      callr::r(
        function(Z, R, W) {
          # Run in a clean session
          library(mkatr)
          mkatr::sats(Z = Z, R = R, W = W)
        },
        args = list(Z, R, W),
        timeout = timeout
      )
    }, error = function(e) {
      message(sprintf("  %s failed: %s", label, conditionMessage(e)))
      NULL
    })
  }
  
  # Load sumstats
  trait <- sub("_weights\\.txt$", "", basename(sumstats_path))
  sumstats <- fread(sumstats_path)
  
  # Load setid and filter
  setid_filtered <- setid[setid$SNP %in% sumstats$SNP, ]
  setid_wide <- split(setid_filtered$SNP, setid_filtered$Gene)
  
  # Filter to only genes in this chunk
  setid_wide <- setid_wide[names(setid_wide) %in% gene_list]
  
  # Define constants
  plink_bin <- "plink"
  twas_dir <- here("data", "raw", "twas", "CMC.BRAIN.RNASEQ")
  
  results_list <- list()
  
  for (gene in names(setid_wide)) {
    tryCatch({
      snps <- setid_wide[[gene]]
      snps_in_data <- intersect(snps, sumstats$SNP)
      
      #Z <- sumstats$Z[sumstats$SNP %in% snps_in_data]
      #names(Z) <- sumstats$SNP[sumstats$SNP %in% snps_in_data]
      #if (!is.numeric(Z) || is.null(names(Z))) next
      
      LD <- ieugwasr::ld_matrix_local(
        variants = snps_in_data,
        bfile = bfile,
        plink_bin = plink_bin,
        with_alleles = FALSE
      )
      #LD <- LD[names(Z), names(Z)]
      snps_order <- snps_in_data
      LD <- LD[snps_order, snps_order]
      
      
      LD <- as.matrix(Matrix::nearPD(LD)$mat) 
      
      
      Z <- MASS::mvrnorm(n = 1, mu = rep(0, length(snps_order)), Sigma = LD)
      names(Z) <- snps_order
      
      W_uniform <- rep(1, length(Z))
      result_uniform <- safe_sats(Z = Z, R = LD, W = W_uniform, timeout = 60, label = paste0("U-", gene))
      
      W_maf <- sumstats$MAF_weight[match(names(Z), sumstats$SNP)]
      if (all(is.na(W_maf)) || sum(abs(W_maf), na.rm = TRUE) == 0) {
        result_maf <- NULL
        W_maf <- rep(0, length(Z))
      } else {
        W_maf <- W_maf / sum(abs(W_maf), na.rm = TRUE)
        W_maf <- W_maf + 10^-5
        W_maf <- W_maf / sum(abs(W_maf), na.rm = TRUE)
        result_maf <- safe_sats(Z = Z, R = LD, W = W_maf, timeout = 120, label = paste0("MAF-", gene))
      }
      
      W_net_degree <- sumstats$degree[match(names(Z), sumstats$SNP)]
      if (all(is.na(W_net_degree)) || sum(abs(W_net_degree), na.rm = TRUE) == 0) {
        result_deg <- NULL
        W_net_degree <- rep(0, length(Z))
      } else {
        W_net_degree <- W_net_degree / sum(abs(W_net_degree), na.rm = TRUE)
        W_net_degree <- W_net_degree + 10^-5
        W_net_degree <- W_net_degree / sum(abs(W_net_degree), na.rm = TRUE)
        result_deg <- safe_sats(Z = Z, R = LD, W = W_net_degree, timeout = 120, label = paste0("DEG-", gene))
      }
      
      W_net_bet <- sumstats$between[match(names(Z), sumstats$SNP)]
      if (all(is.na(W_net_bet)) || sum(abs(W_net_bet), na.rm = TRUE) == 0) {
        result_bet <- NULL
        W_net_bet <- rep(0, length(Z))
      } else {
        W_net_bet <- W_net_bet / sum(abs(W_net_bet), na.rm = TRUE)
        W_net_bet <- W_net_bet + 10^-5
        W_net_bet <- W_net_bet / sum(abs(W_net_bet), na.rm = TRUE)
        result_bet <- safe_sats(Z = Z, R = LD, W = W_net_bet, timeout = 120, label = paste0("BET-", gene))
      }
      
      T_A <- T_S2 <- T_S <- NA
      N_twas <- NA
      result_twas <- NULL
      twas_file <- file.path(twas_dir, paste0("CMC.", gene, ".wgt.RDat"))
      if (file.exists(twas_file)) {
        load(twas_file)
        if (exists("wgt.matrix") && !is.null(wgt.matrix)) {
          twas_weights <- as.data.table(wgt.matrix, keep.rownames = "SNP")
          common_snps <- intersect(twas_weights$SNP, names(Z))
          if (length(common_snps) >= 1) {
            N_twas <- length(common_snps)
            ZT <- Z[common_snps]
            LDT <- LD[common_snps, common_snps]
            W_twas_blup <- twas_weights[match(common_snps, twas_weights$SNP), blup]
            W_twas_blup[is.na(W_twas_blup)] <- 0
            if (sum(abs(W_twas_blup)) > 0 && all(is.finite(W_twas_blup))) {
              W_twas_blup <- W_twas_blup / sum(abs(W_twas_blup))
              W_twas_blup <- W_twas_blup + 10^-5
              W_twas_blup <- W_twas_blup / sum(abs(W_twas_blup))
              result_twas <- safe_sats(Z = ZT, R = LDT, W = W_twas_blup, timeout = 120, label = paste0("TWAS-", gene))
              if (!is.null(result_twas)) {
                T_A  <- result_twas$p.value["A"]
                T_S2 <- result_twas$p.value["S2"]
                T_S  <- result_twas$p.value["S"]
              }
            }
          }
        }
      }
      
      N_uniform <- length(Z)
      N_maf     <- sum(!is.na(W_maf))
      N_deg     <- sum(abs(W_net_degree) > 0, na.rm = TRUE)
      N_bet     <- sum(abs(W_net_bet) > 0, na.rm = TRUE)
      
      result_row <- data.frame(
        Gene     = gene,
        N_uniform = N_uniform,
        N_maf     = N_maf,
        N_twas    = N_twas,
        N_deg     = N_deg,
        N_bet     = N_bet,
        
        U_A  = if (!is.null(result_uniform)) result_uniform$p.value["A"]  else NA,
        U_S2 = if (!is.null(result_uniform)) result_uniform$p.value["S2"] else NA,
        U_S  = if (!is.null(result_uniform)) result_uniform$p.value["S"]  else NA,
        
        M_A  = if (!is.null(result_maf)) result_maf$p.value["A"]  else NA,
        M_S2 = if (!is.null(result_maf)) result_maf$p.value["S2"] else NA,
        M_S  = if (!is.null(result_maf)) result_maf$p.value["S"]  else NA,
        
        T_A  = if (!is.null(result_twas)) result_twas$p.value["A"]  else NA,
        T_S2 = if (!is.null(result_twas)) result_twas$p.value["S2"] else NA,
        T_S  = if (!is.null(result_twas)) result_twas$p.value["S"]  else NA,
        
        D_A  = if (!is.null(result_deg)) result_deg$p.value["A"]  else NA,
        D_S2 = if (!is.null(result_deg)) result_deg$p.value["S2"] else NA,
        D_S  = if (!is.null(result_deg)) result_deg$p.value["S"]  else NA,
        
        B_A  = if (!is.null(result_bet)) result_bet$p.value["A"]  else NA,
        B_S2 = if (!is.null(result_bet)) result_bet$p.value["S2"] else NA,
        B_S  = if (!is.null(result_bet)) result_bet$p.value["S"]  else NA
      )
      
      # Write immediately
      out_file <- file.path(here("data", "processed", "results_simulation2"), paste0(trait, "_chunk_", chunk_id, ".txt"))
      fwrite(result_row,
             file = out_file,
             sep = "\t",
             append = file.exists(out_file),
             col.names = !file.exists(out_file))
      
      cat("Saved result for gene:", gene, "\n")
      
    }, error = function(e) {
      message("Error in gene: ", gene)
      message("  ", conditionMessage(e))
    })
  }
  
  cat("Finished chunk", chunk_id, "for trait", trait, "\n")
  return(NULL)
}

