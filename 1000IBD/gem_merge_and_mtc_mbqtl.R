#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: gem_merge_and_mtc_mbqtl.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(qvalue)
library(stats)
library(optparse)


####################
# Functions        #
####################

#' add nominal p values cutoffs for each feature
#' 
#' @param input_table the QTL output table with the featuress
#' @param pval_column the nominal p value column to calculate the cutoff for
#' @param qval_column the column with the q values
#' @param fdr at which fdr to get the nominal p value cutoff for
#' @returns the input table, with the nominal p value cutoff added for each feature
#' 
add_nominal_cutoff <- function(input_table, pval_column='pval_beta', qval_column='qval', fdr=0.05) {
  # get the p values for effects that were significant after qvalue correction on the permuted p values
  lb_significant <- input_table[input_table[[qval_column]] <= fdr, pval_column]
  # order to get the smallest p first
  lb_significant <- lb_significant[order(lb_significant)]
  # get the p values for effects that were not significant after qvalue correction on the permuted p values
  ub_nonsignificant <- input_table[input_table[[qval_column]] > fdr, pval_column]
  # order to get the smallest p first
  ub_nonsignificant <- ub_nonsignificant[order(ub_nonsignificant)]
  
  # check if there were any significant effects
  if (length(lb_significant) > 0) {
    # get the biggest p values of the significant effects
    lb <- lb_significant[length(lb_significant)]
    # check if there were any not-significant effects
    if (length(ub_nonsignificant) > 0) {
      # get the smallest p value of the non-significant effects
      ub <- ub_nonsignificant[1]
      # the p value threshold
      pthreshold <- (lb + ub) / 2
    } else {
      # if there are no non-significant effects, the cutoff will be the biggest P value of the significant effects
      pthreshold <- lb
    }
    # fit beta distribution to set the threshold for each gene
    input_table$pval_nominal_threshold <- qbeta(pthreshold, input_table$beta_shape1, input_table$beta_shape2)
  } else {
    # if there were no significant effects, we will set the cutoff at zero
    input_table$pval_nominal_threshold <- 0
  }
  return(input_table)
}


#' add q value corrected value to QTL output
#' 
#' @param input_table the QTL output table with the featuress
#' @param pval_column the p value column to calculate q value for
#' @param qval_column the column with the q values to add
#' @param fdr at which fdr to get the nominal p value cutoff for
#' @returns the input table, with the q value and the nominal p value cutoff added for each feature
#' 
add_qval <- function(input_file, output_file, pval_column='pval_beta', qval_column='qval', fdr=0.05) {
  # read the input
  input_table <- read.table(input_file, header = T, sep = '\t', stringsAsFactors = F)
  # add the qval
  input_table[[qval_column]] <- qvalue(input_table[[pval_column]])$qvalues
  # add nominal significance cutoff
  input_table <- add_nominal_cutoff(input_table, pval_column, qval_column, fdr)
  # set the output location
  output_file_loc <- output_file
  # get the extention of the output file
  extension <- tail(strsplit(output_file, ".", fixed=T)[[1]], 1)
  # gzip the path if it ends with gz
  if (extension == 'gz') {
    output_file_loc <- gzfile(output_file)
  }
  # write the result
  write.table(input_table, output_file_loc, row.names = F, col.names = T, sep = '\t', quote = F)
}

merge_chromosome_outputs <- function(prepend, append, chroms=1:20, verbose=T) {
  # we'll put them all in a list first
  output_per_chrom <- list()
  # check each chromosome
  for (chrom in chroms) {
    # get the location of the file
    chrom_loc <- paste(prepend, chrom, append, sep = '')
    # print progress
    if (verbose) {
      print(paste('reading', chrom_loc))
    }
    # read the file
    output_chrom <- read.table(chrom_loc, header = T, sep = '\t')
    # put in the list
    output_per_chrom[[chrom]] <- output_chrom
  }
  # merge chromosomes
  all_chroms <- do.call('rbind', output_per_chrom)
  return(all_chroms)
}


get_mtc_top_effects <- function(prepend, append='-TopEffects.txt', chroms=1:20, pval_column='BetaAdjustedMetaP', verbose=T) {
  # get the top effects
  top_effects <- merge_chromosome_outputs(prepend, append, chroms)
  # print progress
  if (verbose) {
    print('calculating q values')
  }
  # add the qvalue
  top_effects[['feature_qval']] <- qvalue(top_effects[[pval_column]])$qvalues
  return(top_effects)
}


merge_and_mtc <- function(prepend, top_append='-TopEffects.txt', all_append='-AllEffects.txt.gz', chroms=1:20, pval_column='BetaAdjustedMetaP', verbose=T) {
  # get top effects with mtc
  mtc_top_effects <- get_mtc_top_effects(prepend, top_append, chroms, pval_column)
  # get the nominal effects
  all_effects <- merge_chromosome_outputs(prepend, all_append, chroms)
  # print progress
  if (verbose) {
    print('merging mtc top effects onto nominal effects')
  }
  # now merge on the feature
  top_effects_to_combine <- mtc_top_effects[match(all_effects[['Gene']], mtc_top_effects[['Gene']]), c('NrTestedSNPs','ProportionBetterPermPvals','BetaDistAlpha','BetaDistBeta','BetaAdjustedMetaP', 'feature_qval')]
  # now merge onto all the effects
  all_effects <- cbind(all_effects, top_effects_to_combine)
  return(all_effects)
}


write_merge_and_mtc <- function(prepend, top_append='-TopEffects.txt', all_append='-AllEffects.txt.gz', chroms=1:20, pval_column='BetaAdjustedMetaP', verbose=T) {
  # get the combined effects
  all_effects_merged <- merge_and_mtc(prepend, top_append, all_append, chroms, pval_column)
  # set an output file
  output_loc <- gzfile(paste(prepend, 'all_qval.tsv.gz', sep = ''))
  # print progress
  if (verbose) {
    print(paste('writing output to', output_loc))
  }
  # write the result
  write.table(all_effects_merged, output_loc, row.names = F, col.names = T, sep = '\t', quote = F)
}

# debug function
do_debug <- function() {
  write_merge_and_mtc(
    prepend='/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/ongoing/qtl/eqtl/output/bulk_', 
    top_append='-TopEffects.txt', 
    all_append='-AllEffects.txt.gz', 
    chroms=1:20, 
    pval_column='BetaAdjustedMetaP')
}

####################
# debug code      #
####################

#do_debug()


####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character",
              help="input base path", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# load parameters
input <- opt[['input']]

# do the run
write_merge_and_mtc(input)