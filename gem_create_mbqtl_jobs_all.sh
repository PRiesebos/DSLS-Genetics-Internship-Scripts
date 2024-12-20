#!/usr/bin/env bash
###################################################################
#Script Name	  : gem_create_mbqtl_jobs_per_chrom_all.sh
#Description	  : create SBATCH jobs scripts to do tensorQTL mapping
#Args           :
#Author       	: Roy Oelen
#example        : ./gem_create_mbqtl_jobs_all.sh \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/genotype/werna_merged_filtered_chrs.vcf.gz \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/metadata/ \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/qtl/eqtl/output/ \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/qtl/eqtl/jobs/
###################################################################

# the parameters supplied
GENOTYPES_LOC=$1
FEATURES_LOC=$2
OUTPUT_LOC=$3
JOBS_LOC=$4

# the regex we use to get the files
FEATURES_REGEX='*.tsv.gz'

# and some system settings we usually don't have to change
RUNTIME='4:00:00'
CORES='16'
MEMORY_GB='32'

# and libraries
JAVA_LIB='Java/11-LTS'

# and the executable
MBQTL_LOC='/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/MbQTL-1.5.0-SNAPSHOT-jar-with-dependencies.jar'

# location of GTE
GTE_LOC='/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/qtl/eqtl/annotations/g2e2d_1000IBD_dedup.tsv'
ANNOTATIONS_LOC='/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/references/gencode_44_2023/annotation_file_build44_genes_no_version.tsv'

# Exp file override:
EXP_FILE='/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/combined_expression_matrix_protein_coding_filtered.txt.gz'

# settings related to the mapping
MAF_THRESHOLD='0.05'
SEED='7777'
N_PERM='100'

# the chromosomes
CHROMS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22')

# loop through files
files=(${FEATURES_LOC}"/"${FEATURES_REGEX})
# list all files that have the search
for f in ${files[@]}
do
        # ${f} is the features file full path, we also need just the name
        fbasename=$(basename "${f}" .tsv.gz)
        # make basename posix compliant (i.e. no spaces etc.)
        posix_basename=${fbasename//[^a-zA-Z0-9]/_}
        # now do each chromosome
        for chrom in ${CHROMS[@]}
        do
            # create base for the result files
            result_location=${OUTPUT_LOC}"/"${posix_basename}'_'${chrom}
            # create location for job file
            job_location=${JOBS_LOC}"/"${posix_basename}"_"${chrom}"_mbqtl_job_SBATCH.sh"
            # create header of sbatch file
            echo "#!/usr/bin/env bash
#SBATCH --job-name=map_"${posix_basename}"_"${chrom}"
#SBATCH --output="${JOBS_LOC}"/map_"${posix_basename}"_"${chrom}".out
#SBATCH --error="${JOBS_LOC}"/map_"${posix_basename}"_"${chrom}".err
#SBATCH --time="${RUNTIME}"
#SBATCH --cpus-per-task="${CORES}"
#SBATCH --mem="${MEMORY_GB}"GB
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

" > ${job_location}

# add the libraries
echo "ml "${JAVA_LIB}"

" >> ${job_location}

# add the permutation pass
echo "java -jar "${MBQTL_LOC}" \\
    --vcf "${GENOTYPES_LOC}" \\
    --exp "${EXP_FILE}" \\
    --annotation "${ANNOTATIONS_LOC}" \\
    --gte "${GTE_LOC}" \\
    --mode mbqtl \\
    --maf "${MAF_THRESHOLD}" \\
    --seed "${SEED}" \\
    --out "${result_location}" \\
    --perm "${N_PERM}" \\
    --chr "${chrom}" \\
    --mingenotypecount 2 \\
    --outputall
" >> ${job_location}

        done
    
done


# ./gem_create_mbqtl_jobs_all.sh \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/genotype/werna_merged_filtered_chrs.vcf.gz \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/rna/qc/output/ \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/qtl/eqtl/output_no_ver \
# /groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/qtl/eqtl/jobs

