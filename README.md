# DSLS Master's Internship Scripts - UMCG Genetics Department

This repository contains a collection of scripts developed to assist with various tasks during the DSLS Master's internship at the UMCG Genetics Department. These scripts streamline and automate tasks such as data preparation, pipeline job submission, and data inspection.

## Overview of Tasks

The scripts are used for:

- **Data Preparation**:  
  Formatting, cleaning, and organizing data to ensure readiness for analysis or integration into workflows.

- **Pipeline Job Submission**:  
  Automating the submission of computational jobs to pipelines, including setting configurations, managing dependencies, and tracking progress.

- **Data Inspection**:  
  Validating, exploring, and examining datasets to ensure data integrity and extract insights.

---

## Script Details

Below is a description of each script, including its purpose, expected input, and generated output.

### **General Python and Shell Scripts**
- **`CombineQCFiles.py`**  

  **_Purpose:_**  
  This script consolidates and integrates various quality metrics and statistics from alignment and genotyping files. The processed data is output as a tab-delimited file summarizing metrics for all samples, including RNA-seq metrics, STAR alignment logs, variant metrics, and missingness data.

  **_Input:_**  
  - **Command-line Arguments:**  
    1. `alignmentDir`: Directory containing alignment-related files and metrics.  
    2. `genotypingDir`: Directory containing genotyping-related files and metrics.  
    3. `outfile`: Path for the output file to write the consolidated metrics.  
  - **Alignment Files:**  
    - `multiple_metrics.alignment_summary_metrics.gz`: Alignment summary metrics.  
    - `multiple_metrics.insert_size_metrics.gz`: Insert size metrics.  
    - `*_rnaseqmetrics.gz`: RNA-seq metrics.  
    - `*_ReadsPerGene.out.tab.gz`: STAR read count files.  
    - `*_Log.final.out.gz`: STAR alignment summary logs.  
  - **Genotyping Files:**  
    - `output_stats_regions.tsv.gz`: Variant statistics file.  
    - `*_PRE_FILTER.imiss`: Missingness file before filtering.  
    - `filtered_output.vcf.gz-filtered.vcf.gz`: Filtered VCF file.  
    - `*_POST_FILTER.imiss`: Missingness file after filtering.  
  - Optional: `samplesWithPrediction_16_09_22_noOutliers.txt` for additional metadata like study and layout.

  **_Output:_**  
  - A tab-delimited file summarizing metrics for all samples, with samples as columns and metrics as rows.  
    - Example metrics include:  
      - `ALIGNMENT_METRICS_*` (e.g., alignment summary stats)  
      - `INSERT_METRICS_*` (e.g., median insert size)  
      - `RNASEQ_METRICS_*` (e.g., RNA-seq quality stats)  
      - `STAR_TAB_*` (e.g., STAR gene counts)  
      - `STAR_LOG_*` (e.g., STAR alignment logs)  
      - `TOTAL_VARIANTS`, `MISSINGNESS`, `VARIANTS_WITH_GENOTYPE_CALL`  
      - Filtered variants and missingness metrics (`*_AFTER_FILTERING`).  
      - Optional: `STUDY`, `LAYOUT`.  

- **`MergeMetrics.py`**  

  **_Purpose:_**  
  Consolidates multiple gzipped TSV files containing metrics into a single gzipped TSV file, aligning data by shared metric identifiers.

  **_Input:_**  
  - **Command-line Argument:**  
    - `-i` or `--input_path`: Path to the directory containing `_outfile.txt.gz` files from the CombineQCFiles.py script.  
  - **File Format:**  
    - Tab-delimited, with a `-` column for metric names and additional columns for sample values.

  **_Output:_**  
  - **Merged File:**  
    - `metric_files_merged.tsv.gz` saved in the same directory, containing all metrics as rows and sample names as columns.


- **`PCA_expression_data.py`**  

  **_Purpose:_**  
  Performs PCA (Principal Component Analysis) on the correlation matrix of the input data to reduce dimensions and visualize component relationships.

  **_Input:_**  
  - **Command-line Arguments:**  
    - `infile`: Path to the input tab-separated matrix file (samples as rows, features as columns).  
    - `outprefix`: Prefix for output files.  
    - `npcs`: Number of principal components to compute.
  - **File Format:**  
    - Tab-separated values with sample names as row indices.

  **_Output:_**  
  - **Principal Component Data:**  
    - `outprefix_PCs.txt`: Tab-separated file with computed principal components as columns and features as rows.  
  - **Scatter Plots:**  
    - `outprefix_PC1and2.png`: Scatter plot of PC1 vs. PC2.  
    - `outprefix_PC1and4.png`: Scatter plot of PC1 vs. PC4.  


- **`calc_lines_count.py`**  

  **_Purpose:_**  
  Counts the total number of lines across gzipped text files for chromosomes 1 through 21 (left out 22 on purpose).

  **_Input:_**  
  - **Files:**  
    - Files named `Colon_Sigmoid.v8.EUR.allpairs.chr{i}.txt.gz` for `i` ranging from 1 to 21.  
    - Files must be in gzipped text format (`.txt.gz`).

  **_Output:_**  
  - Total line count printed to the console.  
    - Example: `"Total number of lines across all files: {total_lines}"`


- **`combineSampleMapping.py`**  

  **_Purpose:_**  
  Modifies a file containing principal components (PCs) for combined samples by adding study names based on a mapping file.

  **_Input:_**  
  - **Folder:** The folder specified via the command line argument containing the following files:
    1. **`combined_samples_PCs.txt`**  
       - Tab-separated file where the first column contains sample names.
    2. **`study_samples_mapping.txt`**  
       - File containing study names followed by their corresponding sample names (space-separated).

  **_Output:_**  
  - A modified file named **`combined_samples_PCs_with_study.txt`** written to the input folder:
    - Includes all columns from `combined_samples_PCs.txt` with an additional **`Study`** column as the first column.
    - Any unmatched or malformed sample names are logged to the console.


- **`combine_expression_matrices.py`**  

  **_Purpose:_**  
  Combines multiple gzipped tab-delimited expression matrices, aligning them by a common gene identifier column. Handles missing values gracefully during alignment.

  **_Input:_**  
  - A list of gzipped expression matrix files:
    - Each file must have a column for gene identifiers (column named `"-"`) with optional version numbers that will be stripped (e.g., `ENSG00000123456.1` → `ENSG00000123456`).
    - Other columns represent sample expression values.

  **_Output:_**  
  - A single combined gzipped matrix:
    - Genes are merged using an outer join.
    - Missing values are not replaced by zeros in the final output.
    - Written to a specified gzipped file.


- **`create_metrics_matrix_small.py`**  

  **_Purpose:_**  
  Extracts and combines key alignment metrics from two sets of gzipped files located in `multiple_metrics` and `rna_seq_metrics` subdirectories. Generates a single tabular output summarizing `PCT_PF_READS_ALIGNED` and `PCT_CODING_BASES` for each sample.

  **_Input:_**  
  - **Directories:**
    - `multiple_metrics`: Contains `multiple_metrics.alignment_summary_metrics.gz` files for each sample.
    - `rna_seq_metrics`: Contains `*.rnaseqmetrics.gz` files for each sample.
  - **File Contents:**
    - `multiple_metrics.alignment_summary_metrics.gz`: Contains tab-delimited metrics with `PCT_PF_READS_ALIGNED` in the 7th column (PAIR or UNPAIRED rows).
    - `*.rnaseqmetrics.gz`: Contains metrics where `PCT_CODING_BASES` is the 17th value in the second row below a header.

  **_Output:_**  
  - **DataFrame (`alignment_metrics.txt.gz`):**
    - Columns: `Sample`, `PCT_PF_READS_ALIGNED`, `PCT_CODING_BASES`.
    - Rows: One per sample with extracted metrics.
    - Output file is saved as a gzipped tab-delimited text file.
 

- **`create_plink_hardy.sh`**  

  **_Purpose:_**  
  Automates Hardy-Weinberg equilibrium (HWE) analysis for all VCF files in a directory using PLINK 2. Results are saved in a specified output directory.

  **_Input:_**  
  - **VCF Files:**  
    Files matching the pattern `*filtered.vcf.gz` in the current directory.

  **_Output:_**  
  - **Output Directory:**  
    - Results are stored in the `../hardy` directory (created if it doesn't exist).
  - **Files:**  
    For each input VCF file, the script generates files prefixed with the VCF's base filename (without extension) in the output directory. These files contain Hardy-Weinberg equilibrium test results.


- **`create_plink_missing.sh`**  

  **_Purpose:_**  
  Automates the analysis of missing genotype data in VCF files using PLINK 2. Results are stored in a specified output directory.

  **_Input:_**  
  - **VCF Files:**  
    Files matching the pattern `*filtered.vcf.gz` in the current directory.

  **_Output:_**  
  - **Output Directory:**  
    - Results are stored in the `../missing` directory (created if it doesn't exist).
  - **Files:**  
    For each input VCF file, the script generates files prefixed with the VCF's base filename (without extension) in the output directory. These files contain missingness analysis results.


- **`filter_expression_data_on_samples.py`**  

  **_Purpose:_**  
  Filters an expression matrix file to include only columns corresponding to specified samples, writing the filtered result to a new gzipped file.

  **_Input:_**  
  - **Sample File:**  
    A text file (`file.txt`) containing a list of sample names to include.
  - **Expression File:**  
    A gzipped expression matrix (`file.gz`), where the first row contains sample names and subsequent rows contain expression values.

  **_Output:_**  
  - **Filtered Expression File:**  
    A new gzipped file (`file.txt.gz`) containing only the rows and columns corresponding to the specified samples. The output file includes a header with the sample names and corresponding expression data for those samples.


- **`fix_annot.py`**  

  **_Purpose:_**  
  Processes a whitespace-delimited input file, modifies specific columns, and writes the result to an output file.

  **_Input:_**  
  - **Input File:**  
    A whitespace-delimited text file. The script expects at least six columns. It processes the second column (chromosome) by removing the 'chr' prefix and limits the length of the third and fourth columns after the ':' to 10 characters.

  **_Output:_**  
  - **Output File:**  
    A new file with the same format as the input, where the second, third, and fourth columns are modified based on the script's rules. The header is preserved and written to the output file.


- **`fix_folder_structure.sh`**  

  **_Purpose:_**  
  This script creates a folder structure based on a given input directory and organizes various types of data into subfolders. It then copies or moves files from the input directory to the newly created structure.

  **_Input:_**  
  - **Input Location:**  
    The path to the project directory. This directory should contain subdirectories like `fastqc`, `gvcf`, `leafcutter`, `mark_duplicates`, `multiple_metrics`, `rmats`, `rna_seq_metrics`, `star`, and `genotypes`.

  **_Output:_**  
  - A new folder structure is created in the parent directory of the input location, containing subfolders for each sample and data type (e.g., `fastqc`, `gvcf`, `leafcutter`, etc.).  
  - Data from the input directory is copied or moved into these subfolders, depending on the file type and organization specified by the script.
 

- **`gem_create_mbqtl_jobs_all.sh`**  

  **_Purpose:_**  
  This script automates the creation of SBATCH job scripts to perform mbQTL mapping for multiple chromosomes using genotype and gene expression data. It generates jobs for each chromosome (1-22) based on feature files in the input directory and prepares them for submission on a SLURM cluster.

  **_Input:_**  
  - **Genotype Location:**  
    Path to the VCF file containing genotype data for the mapping.
  
  - **Features Location:**  
    Path to the directory containing feature files (e.g., gene expression files in `.tsv.gz` format) to be used for mapping. The script will look for files matching `*.tsv.gz`.

  - **Output Location:**  
    Path to the directory where the results from mbQTL mapping will be stored, with results for each chromosome saved in separate subdirectories.

  - **Jobs Location:**  
    Path to the directory where SBATCH job scripts will be written for each chromosome and feature file.

  **_Output:_**  
  - A set of SBATCH job scripts is created, one per chromosome and feature file, in the specified **Jobs Location**. Each script is configured to run mbQTL mapping on a specific chromosome.
  - Results for each chromosome will be stored in the **Output Location**, with files named based on the feature and chromosome (e.g., `feature_1` for chromosome 1).


- **`linkSampleToStudy.py`**  

  **_Purpose:_**  
  This script identifies and maps sample names from study directories within a root directory. It searches for files matching the pattern `ReadsPerGene.out.tab.gz`, and based on the file names, associates them with their corresponding study. The script outputs a mapping of study names and the associated sample names to a text file.

  **_Input:_**  
  - **Root Directory:**  
    The path to the root directory containing multiple study subdirectories. These subdirectories should start with `SRP` or `ERP`, indicating the study identifiers.
  
  **_Output:_**  
  - A text file (`study_samples_mapping.txt`) containing a mapping of study names and their associated sample names. Each line corresponds to a study, followed by the sample names (either starting with `SRR` or `ERR`) found within that study directory.
 

- **`move_expressionqc_output.sh`**  

  **_Purpose:_**  
  This script copies the `expressionqc_output` folder from multiple study directories into a final destination directory. It checks whether the source folder exists before copying it. If the folder is missing, the script prints a message indicating the absence of the directory.

  **_Input:_**  
  - **Source Folders:**  
    A list of study directories, including both `SRP` and `ERP` study codes (e.g., `SRP076426`, `SRP077046`, `ERP114636`, etc.). These directories should contain a subfolder named `expressionqc_output`.
  
  - **Base Source Path:**  
    The path to the root directory where the study folders are located.
  
  - **Destination Path:**  
    The target directory where the `expressionqc_output` folder for each study will be copied.

  **_Output:_**  
  - Copies the `expressionqc_output` folder from each study directory into the destination path. If a directory doesn't exist, an error message is printed for that specific study.


- **`move_filtered_samples_metrics.py`**  

  **_Purpose:_**  
  This script is designed to move specific samples (identified by sample names) from subdirectories (`star`, `multiple_metrics`, `rna_seq_metrics`) to a separate "removed_samples" directory. It checks the subdirectories for these sample names and moves the corresponding folders to a designated "removed_samples" folder while maintaining the original directory structure.

  **_Input:_**  
  - **Base Directory:**  
    The base directory path where the subdirectories (`star`, `multiple_metrics`, `rna_seq_metrics`) exist. This is supplied as a command-line argument.

  - **Samples to Remove:**  
    A list of sample names to be removed (e.g., `ERR3262403`, `ERR3262441`, `ERR3262409`). The script will check if any subfolder name contains one of these sample names and move them accordingly.

  **_Output:_**  
  - The specified samples (those in the `samples_to_remove` list) are moved from their respective subdirectories to new subdirectories under a `removed_samples` parent directory, within `star`, `multiple_metrics`, and `rna_seq_metrics` subdirectories.  
  - The script creates any missing directories if needed, including subdirectories for each category of removed sample (e.g., `star`, `multiple_metrics`, `rna_seq_metrics`).


- **`remove_filtered_samples_folders.py`**  

  **_Purpose:_**  
  This script removes unwanted sample folders from specified directories. It reads a list of sample identifiers from a file (`first_parts.txt`), and retains only the folders corresponding to these samples. Any folders that do not match the identifiers are deleted from the `star`, `multiple_metrics`, and `rna_seq_metrics` directories.

  **_Input:_**  
  - **Directories:**  
    Three directories are specified:  
    - `star_dir`: Path to the `star` directory
    - `multiple_metrics_dir`: Path to the `multiple_metrics` directory
    - `rna_seq_metrics_dir`: Path to the `rna_seq_metrics` directory

  - **File (`first_parts.txt`):**  
    This file contains a list of sample identifiers to be kept. Only folders with names matching any of these sample identifiers will be retained in the directories.

  **_Output:_**  
  - The script deletes any subfolders in the specified directories (`star`, `multiple_metrics`, `rna_seq_metrics`) that do not contain a sample identifier from the `first_parts.txt` file.  
  - It performs the deletions by using the `rm -rf` command to remove unwanted folders.


- **`run_imputation.py`**  

  **_Purpose:_**  
  This Python script generates and submits SLURM job scripts for imputation (eQTLgen phase 2 imputation pipeline, https://github.com/eQTLGen/eQTLGenImpute) tasks using `Nextflow` and `Singularity`. It creates an `sbatch` script for each study ID, which runs the imputation for each study. The script submits jobs sequentially, ensuring that each new job depends on the successful completion of the previous job.

  **_Input:_**  
  - **Study IDs:**  
    A list of study directories representing the studies to process.  
  - **Base QC Input Folder:**  
    Path to the base directory containing the study-specific data.
  - **Base Output Path:**  
    Path to the base directory where the output of the imputation will be stored.

  **_Output:_**  
  - SLURM job scripts are generated for each study, and these are submitted to the SLURM scheduler.  
  - Each job is submitted with a dependency on the previous job, ensuring that the jobs run sequentially.
  - The SLURM job ID for each submission is printed to the terminal, and the scripts are saved in the `sbatch_scripts/imputation` directory.
  - Imputed genotype data is stored to the output directory.

  **_Generated Files:_**  
  - A separate SLURM submission script (`submit_<study_id>.sh`) for each study ID is created and stored in the `sbatch_scripts/imputation` directory.
  - Log and error files are generated for each job run (e.g., `chr{chromosome}.log`, `chr{chromosome}.err`).


- **`run_mapping.py`**  

  **_Purpose:_**  
  This Python script generates and submits SLURM job scripts for processing genetic data for each chromosome (from 1 to 22) using the `mbQTL` tool (https://github.com/molgenis/systemsgenetics/tree/master/mbQTL). The script creates a separate SLURM job for each chromosome, which runs the tool to perform eQTL mapping. The jobs are submitted sequentially, ensuring that each job depends on the successful completion of the previous one.

  **_Input:_**  
  - **Chromosome Range:**  
    A range of chromosomes (1 through 22) to process, with each chromosome having its own SLURM job.
  - **`MbQTL` Tool:**  
    The tool (`MbQTL-1.5.0-SNAPSHOT-jar-with-dependencies.jar`) used to perform the analysis on the chromosome data.
  - Output path where mbQTL output will be stored.
  - **Input Files:**  
    - **VCF file:** Path to the VCF file for each chromosome (e.g., `chr{chromosome}.dose.filtered.vcf.gz`).
    - **Expression Data File:** A file containing expression data.
    - **Annotation File:** An annotation file used by `mbQTL` for mapping.
    - **Link File:** A file mapping the genotype identifiers to expression identifiers, and the dataset the samples are from.

  **_Output:_**  
  - SLURM job scripts are generated for each chromosome (`submit_chr{chromosome}.sh`) and stored in the `sbatch_scripts` directory.  
  - These scripts are submitted to the SLURM job scheduler with dependencies ensuring sequential execution based on successful completion of previous jobs.
  - The output for each job is written to log files in the `sbatch_scripts` directory (e.g., `chr{chromosome}.log`, `chr{chromosome}.err`).
  - Output from the mbQTL mapping software is stored to the output path.
  
  **_Generated Files:_**  
  - SLURM job submission scripts (`submit_chr{chromosome}.sh`) are created for each chromosome and stored in the `sbatch_scripts` directory.
 

- **`run_pubrna.py`**

  **_Purpose:_**  
  Automates the generation and submission of SLURM job scripts to run a Nextflow RNA-Seq pipeline (https://github.com/ogkourlias/pub-rna
) on multiple datasets. It ensures that each dataset is processed sequentially with job dependencies.

  **_Input:_**  
  - A list of input directories (e.g., `SRP076426`, `SRP077046`, etc.) corresponding to the RNA sequencing datasets to be processed.  
  - Path to the Nextflow pipeline script (`main.nf`).  
  - Reference files required for the pipeline, including genome reference (`GRCh38`), GTF annotations, and other necessary data files.  

  **_Output:_**  
  - Generated SLURM job scripts for each dataset, stored in the `sbatch_scripts/` directory. Each script is named `submit_{input_dir}.sh`.  
  - SLURM log files for each job:  
    - Standard output: `sbatch_scripts/{input_dir}.log`  
    - Standard error: `sbatch_scripts/{input_dir}.err`  
  - Sequential job submissions, with each job dependent on the completion of the previous one.
  - Aligned and genotyped files for each SRA ID study, with log files (STAR, rna seq metrics, input metrics, ...)  


- **`run_pubrnapipeline.sh`**  

  **_Purpose:_**  
  Automates the creation and execution of Nextflow pipeline jobs for multiple RNA-Seq datasets. It generates a custom Nextflow configuration file for each dataset and runs the RNA-Seq pipeline (`pub-rna`, https://github.com/ogkourlias/pub-rna
) on each dataset sequentially.

  **_Input:_**  
  - A list of input directories (`SRP063496`, `SRP068609`, `SRP076426`, etc.), each representing an RNA sequencing dataset to be processed.  
  - A template Nextflow configuration file (`nextflow.config`) containing the pipeline's configuration.  
  - The Nextflow pipeline script (`main.nf`) used to run the RNA-Seq analysis.  

  **_Output:_**  
  - For each dataset, a custom Nextflow configuration file (`temp_nextflow_{input_dir}.config`) is generated, with dataset-specific parameters like the output directory and input text file.  
  - The Nextflow job runs in the background for each dataset, creating corresponding output directories for each dataset's processed data.  
  - Optionally, the temporary configuration files can be removed after job execution.
  - Aligned and genotyped files for each SRA ID study, with log files (STAR, rna seq metrics, input metrics, ...)  


- **`run_qcpipeline.py`**  

  **_Purpose:_**  
  Automates the creation and submission of SLURM job scripts to run the QC analysis pipeline (`QCPipeline`, https://github.com/PRiesebos/QCPipeline) for multiple datasets using Nextflow. It ensures that each dataset is processed sequentially, with job dependencies between them.

  **_Input:_**  
  - A list of input directories (e.g., `ERP109626`, `ERP113396`, etc.), each corresponding to an RNA sequencing dataset to be processed.  
  - The Nextflow pipeline script (`main.nf`) used for the QC analysis.  
  - A template configuration file for Nextflow (`nextflow.config`) that includes various parameters required for the QC analysis.

  **_Output:_**  
  - Generated SLURM job scripts for each dataset, stored in the `sbatch_scripts/QC/` directory. Each script is named `submit_{input_dir}.sh`.  
  - SLURM log files for each job:  
    - Standard output: `sbatch_scripts/QC/{input_dir}.log`  
    - Standard error: `sbatch_scripts/QC/{input_dir}.err`  
  - Sequential job submissions, with each job dependent on the completion of the previous one.
  - genotype QCPipeline output files, PLINK2 format bed files for each study.  


- **`run_qcpipeline.sh`**  

  **_Purpose:_**  
  Automates the execution of the QC analysis pipeline (`QCPipeline`, https://github.com/PRiesebos/QCPipeline) using Nextflow for multiple datasets. This script runs the pipeline sequentially for each dataset, generating necessary configuration files and executing the pipeline with appropriate parameters.

  **_Input:_**  
  - A list of input directories (e.g., `SRP64952`, `SRP076426`, etc.), each corresponding to a dataset that needs to be processed.  
  - A template Nextflow configuration file (`nextflow.config`) that is modified for each dataset.  
  - The Nextflow pipeline script (`main.nf`) to run the QC analysis.

  **_Output:_**  
  - For each dataset, the script modifies the Nextflow configuration and submits the job to Nextflow for execution.  
  - Logs and outputs from the Nextflow job will be generated by the pipeline itself (i.e., standard output and error logs).  
  - A one-hour delay (`sleep 1h`) between each job to control execution timing.
  - genotype QCPipeline output files, PLINK2 format bed files for each study.  


- **`sample_study_mapping.py`**  

  **_Purpose:_**  
  Processes the alignment metrics data and associates each sample with its corresponding study based on the `study_samples_mapping.txt` file. This script reads the alignment metrics, adds the study information to the dataset, and saves the resulting dataframe to a compressed file.

  **_Input:_**  
  - `alignment_metrics.txt.gz`: A gzipped tab-delimited file containing alignment metrics (e.g., `PCT_CODING_BASES`, `PCT_PF_READS_ALIGNED`, etc.) for RNA-Seq samples.  
  - `study_samples_mapping.txt`: A text file mapping samples to their respective studies.

  **_Output:_**  
  - `alignment_metrics_with_study.txt.gz`: A gzipped tab-delimited file that contains the original alignment metrics, but with an additional column (`Study`) mapping each sample to its respective study.
 

- **`study_cleanup_remove.sh`**  

  **_Purpose:_**  
  Cleans up specific directories based on a list of identifiers. For each identifier, it attempts to remove the corresponding directory (e.g., `SRP077046-test`) if it exists.

  **_Input:_**  
  - A list of identifiers (e.g., `SRP077046`, `SRP096757`, etc.) that correspond to the directories to be cleaned.  
  - The script operates on directories named using the identifier appended with `-test` (e.g., `SRP077046-test`).

  **_Output:_**  
  - Directories corresponding to each identifier are removed (if they exist). No new files are created.


- **`submit_mbqtl.py`**  

  **_Purpose:_**  
  Automates the creation and submission of SLURM job scripts to run mbQTL (https://github.com/molgenis/systemsgenetics/tree/master/mbQTL) analysis across multiple chromosomes. Each chromosome job runs the mbQTL analysis on a specific chromosome, using various genomic and expression data files.

  **_Input:_**  
  - Genomic VCF file (`1000IBD_merged_grch38.vcf.gz`) containing genetic variants.  
  - Expression data file (`combined_expression_matrix_protein_coding_filtered_no_zeros.txt.gz`) containing expression levels for the genes of interest.  
  - Annotation file (`annotation_file_build44_genes_no_version.tsv`) for the genome version.  
  - A link file (`1000IBD_new_g2e_dedup.tsv`) linking genotype data to expression data.  
  - Various parameters like mode (`mbqtl`), number of permutations (100), minimum genotype count (2), and minor allele frequency (MAF) (0.05).
  - Output path where mbQTL output will be stored.

  **_Output:_**  
  - SLURM job submission scripts for each chromosome (1-22), saved in the `werna_og_geno_liftover_grch38` directory. Each script is named `submit_chr{chr_num}.sh`.
  - SLURM log files for each job:  
    - Standard output: `werna_og_geno_liftover_grch38-{chr_num}.log`  
    - Standard error: `werna_og_geno_liftover_grch38-{chr_num}.err`  
  - Each job is submitted to the SLURM scheduler for execution.
  - Output from the mbQTL mapping software is stored to the output path.


- **`submit_merge_vcf_files.sh`**  

  **_Purpose:_**  
  This script is used to submit a SLURM job to merge multiple VCF files into one combined VCF file using the `bcftools` tool. The merged VCF file is then saved as a compressed `.vcf.gz` file.

  **_Input:_**  
  - path to multiple VCF files, one path for each file

  **_Output:_**  
  - Merged VCF file (`combined_vcf_files.vcf.gz`) containing data from the input VCF files.  
  - SLURM log files:  
    - Standard output: `vcf_merge_%j.out`  
    - Standard error: `vcf_merge_%j.err`


- **`tabix_script.sh`**  

  **_Purpose:_**  
  This script automates the process of running `tabix` on all VCF files with a common prefix (`chr*`). It ensures that each file is indexed using `tabix` for efficient access.

  **_Input:_**  
  - VCF files with a filename starting with `chr` (e.g., `chr1.vcf`, `chr2.vcf`, etc.).

  **_Output:_**  
  - Indexed VCF files using `tabix`. Each file will have a corresponding index file (`.tbi`).  
  - Console output indicating any missing files (if a file starting with `chr` does not exist).


- **`tar_gzip_cleanup.sh`**  

  **_Purpose:_**  
  This script compresses the log directories associated with multiple identifiers into tarball files. For each identifier, it creates a tar.gz archive of the corresponding directory.

  **_Input:_**  
  - A list of identifiers (e.g., `SRP063496`, `SRP064952`, etc.).  
  - Corresponding directories named as `{identifier}-test` (e.g., `SRP063496-test`, `SRP064952-test`, etc.) containing log files.

  **_Output:_**  
  - A tarball for each identifier's directory (e.g., `SRP063496-logs.tar.gz`, `SRP064952-logs.tar.gz`, etc.). Each tarball contains the contents of the associated `{identifier}-test` directory.


---

### **coloc Folder**
- **`coloc.R`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`filter_eqtls_for_coloc.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`filter_gwas_for_eqtlsnps_for_coloc.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

---

### **eQTL_Scripts Folder**
- **`1-merge_top_effects.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`2-calculate_q_values.R`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`3-count_significant_effects.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

---

### **1000G_comparison_scripts Folder**
- **`step_1_filter_comparison.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`step_1_filter_comparison_test.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`step_2_freq_calc_plink.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`step_3_fix_annot.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

---

### **1000G_comparison_scripts (fixed_maf) Folder**
- **`step_1_filter_comparison_fixed_maf.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`step_2_freq_calc_plink_fixed_maf.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`step_3_fix_annot_fixed_maf.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

---

### **jupyter_notebooks Folder**
- **`1000IBD_paper_vs_own_data_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`1000IBD_sample_comparison.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`1000gComparison.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`1000gComparison_fixed_maf.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`MergedMetricsInspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`QcPlots.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`QcPlots_2.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`SRP113470_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`SRP113470_metrics_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`analyze_imputed.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`coloc_final_prep_and_output_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`coloc_prep.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`combined_count_matrix_pca.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`expression_files_comparison.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`filter_expression_data_on_protein_coding.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`filter_gencode_protein_coding.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`filter_mbQTL_output_protein_coding.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`final_mbqtl_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`gtex_eqtl_comparison_fastqtl_vs_mbqtl.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`gtex_mbqtl_preparations.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`mbQTL_output_protein_coding_filtered_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`mbqtl_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`mbqtl_pub_rna_output_inspection.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`merge_expression_data.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`minimac4_chunkSize_comparison.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`summarystats_top_vs_all_graphs.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`vcf_analysis.ipynb`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

---

### How to Use

1. Clone this repository to your local machine or server.
2. Navigate to the relevant folder containing the scripts you wish to use.
3. Make sure you have the required input files, and relevant tools/libraries.
4. Run the script according to the instructions provided in the script or within its corresponding notebook.

