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
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`MergeMetricFiles.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`PCA_expression_data.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`calc_lines_count.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`combineSampleMapping.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`combine_expression_matrices.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`create_metrics_matrix_small.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`create_plink_hardy.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`create_plink_missing.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`filter_expression_data_on_samples.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`fix_annot.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`fix_folder_structure.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`gem_create_mbqtl_jobs_all.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`linkSampleToStudy.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`move_expressionqc_output.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`move_filtered_samples_metrics.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`remove_filtered_samples_folders.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`run_imputation.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`run_mapping.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`run_pubrna.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`run_pubrnapipeline.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`run_qcpipeline.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`run_qcpipeline.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`sample_study_mapping.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`study_cleanup_remove.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`submit_mbqtl.py`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`submit_merge_vcf_files.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`tabix_script.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

- **`tar_gzip_cleanup.sh`**  
  _Purpose:_  
  _Input:_  
  _Output:_  

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

