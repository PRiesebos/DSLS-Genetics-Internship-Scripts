#!/usr/bin/env bash
###################################################################
#Script Name	  : gem_liftover_1000IBD_variants.sh
#Description	  : liftover the original 1000IBD genotypes to build 38
#Args           : 
#Author       	: Roy Oelen
# example        : 
###################################################################


# load libraries
ml HTSlib
ml BCFtools
ml picard

# rezip the vcf to have bgzip compression
gunzip -c european_GRCh37_maf0001.vcf.gz | bgzip > european_GRCh37_maf0001.bgz.vcf.gz

# create a sequence dictionary for the reference fasta
java -jar -Xms60g -Xmx60g \
    ${EBROOTPICARD}/picard.jar CreateSequenceDictionary \
      -R /groups/umcg-franke-scrna/tmp04/external_datasets/hg38/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
      -O /groups/umcg-franke-scrna/tmp04/external_datasets/hg38/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.dict

# subset each the vcf for each chromosome, as liftover of all chromosomes together was too memory-intensive
for i in {1..22}
    do
        bcftools view /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/paper_genotype/european_GRCh37_maf0001.bgz.vcf.gz --regions ${i} -o /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/paper_genotype/european_GRCh37_maf0001_chr${i}.vcf.gz -Oz
        tabix -p vcf /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/paper_genotype/european_GRCh37_maf0001_chr${i}.vcf.gz
    done

# now liftover each chromosome to b38
for i in {1..22}
    do
        java -jar -Xms60g -Xmx60g \
        ${EBROOTPICARD}/picard.jar LiftoverVcf \
            -I /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/paper_genotype/european_GRCh37_maf0001_chr${i}.vcf.gz \
            -O /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/paper_genotype/european_GRCh38_lifted_maf0001_chr${i}.vcf.gz \
            --CHAIN /groups/umcg-franke-scrna/tmp04/external_datasets/liftover/hg19ToHg38.over.chain.nochr.gz \
            --REJECT /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/paper_genotype/european_GRCh38_lifted_rejected_maf0001_chr${i}.vcf.gz \
            -R /groups/umcg-franke-scrna/tmp04/external_datasets/hg38/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
            --WARN_ON_MISSING_CONTIG true
    done