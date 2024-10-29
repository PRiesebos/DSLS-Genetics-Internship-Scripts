library(coloc)


# Read data from tab-delimited text files, GWAS sum stats filtered on 5e-8
data1 <- read.table("/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/extra_scripts/coloc/gwas.tsv", header = TRUE, sep = "\t")
data2 <- read.table("/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/extra_scripts/coloc/sum.tsv", header = TRUE, sep = "\t")


trait1 <- list(
  beta = data1$beta,                       # Meta-analysis effect sizes (beta)
  varbeta = data1$varbeta,      # Variance of the beta (SE^2)
  snp = data1$rsid,                        # SNP IDs (Optional but useful for labeling)
  position = data1$base_pair_location,     # SNP positions
  type = "quant",                          # Specify "quant" for quantitative traits
  pvalues = data1$p_value,
  N = 398668,
  MAF = data1$effect_allele_frequency
)

trait2 <- list(
  beta = data2$MetaBeta,                   # Meta-analysis effect sizes (beta)
  varbeta = data2$varbeta,              # Variance of the beta (SE^2)
  snp = data2$SNP,                         # SNP IDs (Optional but useful for labeling)
  position = data2$SNPPos,                 # SNP positions
  type = "quant",                          # Specify "quant" for quantitative traits
  pvalues = data2$MetaP,
  N = data2$MetaPN,
  MAF = data2$SNPEffectAlleleFreq
)

check_dataset(trait1, suffix = "", req = c("type", "snp"), warn.minp = 1e-06)
check_dataset(trait2, suffix = "", req = c("type", "snp"), warn.minp = 1e-06)


# plot_dataset(trait1, main = "Dataset 1: Trait 1")
# check_dataset(trait1)

coloc_results <- coloc.abf(dataset1 = trait1, dataset2 = trait2)

# # View the results
print(coloc_results)

# # Plot the datasets
par(mfrow = c(2, 1))  # Set up the plotting area to have two plots, one above the other
plot_dataset(trait1, main = "Dataset 1: Trait 1")
plot_dataset(trait2, main = "Dataset 2: Trait 2")

# Define the output path based on input directory
output_path <- "/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/extra_scripts/coloc/coloc_results.tsv"

# Write coloc results to file
write.table(as.data.frame(coloc_results$summary), file = output_path, sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(coloc_results$results), file = gsub("coloc_results.tsv", "coloc_detailed_results.tsv", output_path), sep = "\t", quote = FALSE, row.names = TRUE)
