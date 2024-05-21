#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: gem_create_qtl_g2e_1000IBD.R
# Function: create the gene-to-expression file required for mbQTL
############################################################################################################################

####################
# libraries        #
####################

#none


####################
# Functions        #
####################

#none


####################
# Options          #
####################

#none


####################
# Main Code        #
####################

# location of the g2e file used for the imputation
g2e_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/ongoing/metadata/gem_dna_rna_sample_mapping.tsv'
# read the file
g2e <- read.table(g2e_loc, header = F, sep = '\t')
# add the dataset column
g2e[['V3']] <- '1000IBD'
# set the output location
g2e2d_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/ongoing/qtl/eqtl/annotations/g2e2d_1000IBD.tsv'
# write the result
write.table(g2e, g2e2d_loc, sep = '\t', row.names = F, col.names = F, quote = F)
# now create a deduplicated version of the g2e2d
g2e2d_dedup <- g2e[!duplicated(g2e[['V1']]), ]
# and write that result
g2e2d_dedup_loc <- '/groups/umcg-weersma/tmp01/projects/gut_eqtl_meta_analysis/ongoing/qtl/eqtl/annotations/g2e2d_1000IBD_dedup.tsv'
write.table(g2e2d_dedup, g2e2d_dedup_loc, sep = '\t', row.names = F, col.names = F, quote = F)
