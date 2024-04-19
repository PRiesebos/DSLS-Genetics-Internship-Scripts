# gut biopsy scRNA-seq meta-analysis - 1000IBD data
This subfolder contains the scripts used for processing the 1000IBD data

## alignment
Alignment was performed using the RNA-seq pipeline here:
Files were renamed to be compatible using the script '*gem_rename_fastqcs.sh*'

## genotype preprocessing
Genotype imputation was performed using the eQTLgen imputation pipeline here: https://eqtlgen.github.io/eqtlgen-web-site/eQTLGen-p2-cookbook.html#2-genotype-imputation
A genotype to expression annotation file was created using the script '*gem_create_genotype_to_expression_file_1000IBD.R*' (this also merges the gene expression files)

DataQC was done using the template provided by eQTLgen, the modified version is named: '*submit_DataQc_pipeline_1000IBD.sh*'