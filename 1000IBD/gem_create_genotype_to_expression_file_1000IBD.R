#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: gem_create_genotype_to_expression_file_1000IBD.R
# Function: create the genotype to expression file for the imputation pipeline, and merge the expression files for the MDL data 
############################################################################################################################

####################
# libraries        #
####################

# none

####################
# Functions        #
####################

#' for entries in the table missing the GSA_ID, see if the table contains another entry for that sample, where the GSA_ID was given, and use that GSA_ID
#' 
#' @param expression_locations list of paths where to search for expression data
#' @param prepends list of prepended strings placed before the sample names when searching for the expression data
#' @param appends subfolder and file containing the expression data
#' @param samples lists of samples to search expression data for
#' @param count_column which column in the STAR output to use as gene quantification
#' @param verbose print progress
#' @returns a dataframe containing the expression for all samples that were found
#' 
#' example:
#' expression_matrix_1 <- merge_expression_data(
#' expression_locations = c('/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/rna/alignment/output/batch1/star/', '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/rna/alignment/output/batch2/star'),
#' prepends = c('R', 'R'),
#' appends = c('_ReadsPerGene.out.tab.gz', '_ReadsPerGene.out.tab.gz'),
#' samples = list(rna_to_dna[['Novogene.ID.1']], rna_to_dna[['Novogene.ID.1']])
#' )
merge_expression_data <- function(expression_locations, prepends, appends, samples, count_column=2, verbose=T) {
  # create a list with all the samples
  expression_per_sample <- list()
  # keep track of all genes encountered
  all_genes <- c()
  # check each expression location
  for (i in 1 : length(expression_locations)) {
    # get the expression location
    expression_location <- expression_locations[i]
    # the prepend
    prepend <- prepends[[i]]
    # the append
    append <- appends[[i]]
    # the samples
    sample_batch <- samples[[i]]
    # check each sample
    for (sample in sample_batch) {
      # paste together the full path
      path_sample <- paste(expression_location, '/', prepend, sample, '/', prepend, sample, append, sep = '')
      # check if this file exists
      if (file.exists(path_sample)) {
        # print where we are
        if (verbose) {
          print(paste('reading sample', sample, 'at', path_sample))
        }
        # read the file
        expression_sample <- read.table(path_sample, header = F, sep = '\t', skip = 4)
        # sort by the gene
        expression_sample <- expression_sample[order(expression_sample[[1]]), ]
        # now set those as rownames
        rownames(expression_sample) <- expression_sample[[1]]
        # now only keep the counts
        expression_sample <- expression_sample[, c(count_column), drop = F]
        # set the sample as the sole column name
        colnames(expression_sample) <- c(sample)
        # add to the list
        expression_per_sample[[sample]] <- expression_sample
        # and update the list of genes
        new_genes <- setdiff(all_genes, rownames(expression_sample))
        if (length(new_genes) > 0) {
          all_genes <- c(all_genes, new_genes)
        }
      }
      else{
        if (verbose) {
          print(paste('skipping', sample))
        }
      }
    }
  }
  # go through each sample now
  for (sample in names (expression_per_sample)) {
    # get expression for sample
    expression_sample <- expression_per_sample[[sample]]
    # check which genes are missing
    missing_genes <- setdiff(all_genes, rownames(expression_sample))
    # if there are missing ones, we need to add them
    if (length(missing_genes) > 0) {
      # create a new dataframe
      missing_gene_expression <- data.frame(x = rep(0, times = length(missing_genes)))
      # with those genes as rownames
      rownames(missing_gene_expression) <- missing_genes
      # and the sample as column again
      colnames(missing_gene_expression) <- colnames(expression_sample)[1]
      # now add them together
      expression_sample <- rbind(expression_sample, missing_gene_expression)
      # and reorder it
      expression_sample <- expression_sample[order(rownames(expression_sample)), , drop = F]
      # and put that back in the list
      expression_per_sample[[sample]] <- expression_sample
    }
  }
  # merge all the expression
  expression_all <- do.call('cbind', expression_per_sample)
  return(expression_all)
}

####################
# Settings         #
####################

# none


####################
# Main Code        #
####################

# location of UMCG to RNA
umcg_to_rna_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/ongoing/metadata/gem_rna_sample_mapping.tsv'
# location of UMCG to DNA (GSA data)
umcg_to_dna_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/ongoing/metadata/gem_dna_sample_mapping.tsv'

# read the files
umcg_to_rna <- read.table(umcg_to_rna_loc, header = T, sep = '\t', stringsAsFactors = F)
umcg_to_dna <- read.table(umcg_to_dna_loc, header = T, sep = '\t', stringsAsFactors = F)

# merge the dna and rna, but only for the columns we care about
rna_to_dna <- merge(umcg_to_dna[, c('Research.ID', 'DNA.sample.ID')], umcg_to_rna[, c('Research.ID', 'Novogene.ID')], all = T, by = 'Research.ID')
# reformat the Novogene IDs, because the actual sample files are slightly differently named
rna_to_dna[['Novogene.ID']] <- gsub('-', '_', rna_to_dna[['Novogene.ID']])
rna_to_dna[['Novogene.ID.1']] <- gsub('Biopsy_', 'B', rna_to_dna[['Novogene.ID']]) # there are two possibilities
rna_to_dna[['Novogene.ID.2']] <- gsub('Biopsy', 'B', rna_to_dna[['Novogene.ID']]) # with or without the underscore after the B...

# the names in the genotype file are awfull. We'll fix them. Read the genotype file
genotype_fam_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/genotype/qced/ibd_PSIc_mix_mv-qc.fam.original'
genotype_fam <- read.table(genotype_fam_loc, sep = ' ')
# make the sample ID the one in the second column
genotype_fam[['V1']] <- genotype_fam[['V2']]
# save the result
genotype_fam_fixed_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/genotype/qced/ibd_PSIc_mix_mv-qc.fam'
write.table(genotype_fam, genotype_fam_fixed_loc, sep = ' ', quote = F, row.names = F, col.names = F)

# let's combine all the expression data
expression_matrix_1 <- merge_expression_data(
  expression_locations = c('/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/rna/alignment/output/batch1/star/', '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/rna/alignment/output/batch2/star'),
  prepends = c('R', 'R'),
  appends = c('_ReadsPerGene.out.tab.gz', '_ReadsPerGene.out.tab.gz'),
  samples = list(rna_to_dna[['Novogene.ID.1']], rna_to_dna[['Novogene.ID.1']])
)
expression_matrix_2 <- merge_expression_data(
  expression_locations = c('/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/rna/alignment/output/batch1/star/', '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/rna/alignment/output/batch2/star'),
  prepends = c('R', 'R'),
  appends = c('_ReadsPerGene.out.tab.gz', '_ReadsPerGene.out.tab.gz'),
  samples = list(rna_to_dna[['Novogene.ID.2']], rna_to_dna[['Novogene.ID.2']])
)
# merge them
expression_matrix <- merge(expression_matrix_1, expression_matrix_2, by = 0, all = T)
# set the rownames back
rownames(expression_matrix) <- expression_matrix[['Row.names']]
# and remove that column
expression_matrix[['Row.names']] <- NULL
# make NA zero
expression_matrix[is.na(expression_matrix)] <- 0
# add the gene name as an explicit column
expression_matrix <- cbind(data.frame(phenotype_id = rownames(expression_matrix)), expression_matrix)

# save the expression file
expression_loc <- gzfile('/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/processed/rna/alignment/output/merged/1000IBD_expression.tsv.gz')
write.table(expression_matrix, expression_loc, row.names = F, col.names = T, quote = F, sep = '\t')

# make a final column for the expression id
rna_to_dna[['expression_id']] <- NA
# set for novogene id 1 if that one is present in the expression data
rna_to_dna[!is.na(rna_to_dna[['Novogene.ID.1']]) & rna_to_dna[['Novogene.ID.1']] %in% colnames(expression_matrix), 'expression_id'] <- rna_to_dna[!is.na(rna_to_dna[['Novogene.ID.1']]) & rna_to_dna[['Novogene.ID.1']] %in% colnames(expression_matrix), 'Novogene.ID.1']
rna_to_dna[!is.na(rna_to_dna[['Novogene.ID.2']]) & rna_to_dna[['Novogene.ID.2']] %in% colnames(expression_matrix), 'expression_id'] <- rna_to_dna[!is.na(rna_to_dna[['Novogene.ID.2']]) & rna_to_dna[['Novogene.ID.2']] %in% colnames(expression_matrix), 'Novogene.ID.2']
# keep only entries with genotype and expression data
rna_to_dna_complete <- rna_to_dna[!is.na(rna_to_dna[['DNA.sample.ID']]) & !is.na(rna_to_dna[['expression_id']]), c('DNA.sample.ID', 'expression_id')]

# save the genotype to expression file
geno_to_exp_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/ongoing/metadata/gem_dna_rna_sample_mapping.tsv'
write.table(rna_to_dna_complete, geno_to_exp_loc, row.names = F, col.names = F, quote = F, sep = '\t')
