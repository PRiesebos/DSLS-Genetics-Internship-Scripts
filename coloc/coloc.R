#' coloc.r, R script used to run colocalization per gene for GWAS and eQTL variants
#' Author: Peter Riesebos
#' input: tsv files containing both the GWAS and eQTL variants per gene
#' output: colocalization results, table per gene, containing gene names, PPH values, distances and  top found SNP/p-value

library(coloc)

# Define paths
coloc_input_files_path <- "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/coloc_files_newest/"
coloc_output_path <- "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/coloc_output/"

# Get the list of files
coloc_files <- list.files(coloc_input_files_path, pattern = ".*\\.tsv$", full.names = TRUE)
coloc_files <- sort(coloc_files)  # Sort for reproducibility

# Ensure output directory exists
if (!dir.exists(coloc_output_path)) {
  dir.create(coloc_output_path, recursive = TRUE)
}

# Initialize a list to store files with errors
error_files <- list()

# Process each file
for (file in coloc_files) {
  tryCatch({
    # Read input file
    eqtl_gwas_file <- read.table(file, header = TRUE, sep = "\t")
    
    # Extract gene name, symbol, and TSS directly from the file
    gennaam <- eqtl_gwas_file$Gene[1]
    gensymbol <- eqtl_gwas_file$GeneSymbol[1]
    TSS <- eqtl_gwas_file$GenePos[1]
    
    # Prepare eQTL and GWAS data
    eqtl_data <- list(
      beta = eqtl_gwas_file$MetaBeta,
      snp = eqtl_gwas_file$SNP,
      position = eqtl_gwas_file$SNPPos,
      type = "quant",
      pvalues = eqtl_gwas_file$MetaP,
      N = eqtl_gwas_file$MetaPN,
      MAF = eqtl_gwas_file$SNPEffectAlleleFreq
    )
    
    gwas_data <- list(
      beta = eqtl_gwas_file$beta,
      snp = eqtl_gwas_file$rsid,
      position = eqtl_gwas_file$base_pair_location,
      type = "cc",
      pvalues = eqtl_gwas_file$p_value,
      N = 398668,
      s = 45106 / 398668,
      MAF = eqtl_gwas_file$effect_allele_frequency
    )
    
    # Perform coloc analysis
    coloc_result <- coloc.abf(dataset1 = eqtl_data, dataset2 = gwas_data)
    
    # Extract posterior probabilities
    PP3 <- coloc_result$summary["PP.H3.abf"]
    PP4 <- coloc_result$summary["PP.H4.abf"]
    
    # Identify top GWAS SNP
    top_gwas_index <- which.min(gwas_data$pvalues)
    top_GWAS_SNP <- gwas_data$snp[top_gwas_index]
    top_GWAS_pvalue <- gwas_data$pvalues[top_gwas_index]
    top_GWAS_position <- gwas_data$position[top_gwas_index]
    
    # Identify top eQTL SNP
    top_eqtl_index <- which.min(eqtl_data$pvalues)
    top_eQTL_SNP <- eqtl_data$snp[top_eqtl_index]
    top_eQTL_pvalue <- eqtl_data$pvalues[top_eqtl_index]
    top_eQTL_position <- eqtl_data$position[top_eqtl_index]
    
    # Calculate distances
    distance_top_SNPs <- abs(top_GWAS_position - top_eQTL_position)
    distance_top_GWAS_to_TSS <- abs(top_GWAS_position - TSS)
    distance_top_eQTL_to_TSS <- abs(top_eQTL_position - TSS)
    
    # Create a results data frame for the current file
    results <- data.frame(
      gennaam = ifelse(is.na(gennaam), "NA", gennaam),
      gensymbol = ifelse(is.na(gensymbol), "NA", gensymbol),
      PP3 = PP3,
      PP4 = PP4,
      top_GWAS_SNP = ifelse(is.na(top_GWAS_SNP), "NA", top_GWAS_SNP),
      top_GWAS_pvalue = ifelse(is.na(top_GWAS_pvalue), 1, top_GWAS_pvalue),
      top_eQTL_SNP = ifelse(is.na(top_eQTL_SNP), "NA", top_eQTL_SNP),
      top_eQTL_pvalue = ifelse(is.na(top_eQTL_pvalue), 1, top_eQTL_pvalue),
      distance_top_SNPs = ifelse(is.na(distance_top_SNPs), 0, distance_top_SNPs),
      distance_top_GWAS_to_TSS = ifelse(is.na(distance_top_GWAS_to_TSS), 0, distance_top_GWAS_to_TSS),
      distance_top_eQTL_to_TSS = ifelse(is.na(distance_top_eQTL_to_TSS), 0, distance_top_eQTL_to_TSS),
      stringsAsFactors = FALSE
    )
    
    # Generate a valid output filename
    output_file_name <- paste0(
      coloc_output_path, 
      "coloc_summary_", 
      make.names(basename(file)), 
      ".txt"
    )
    
    # Save results to a file
    write.table(results, file = output_file_name, sep = "\t", row.names = FALSE, quote = FALSE)
    print(paste("Summary results saved to:", output_file_name))
    
  }, error = function(e) {
    # Handle errors: Add the file name to the error list and print the error message
    print(paste("Error processing file:", file))
    print(e)
    error_files <<- c(error_files, file)
  })
}

# Save the list of error files
error_file_log <- paste0(coloc_output_path, "error_files_log.txt")
writeLines(as.character(error_files), con = error_file_log)
print(paste("Error files logged to:", error_file_log))
