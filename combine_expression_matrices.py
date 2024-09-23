# Author: Peter Riesebos
# Created: 19-09-2024
# Purpose: script used for combining various expression matrices

import pandas as pd
import gzip

def read_expression_file(filepath):
    """
    Reads a gzipped tab-delimited expression matrix file into a pandas DataFrame.
    The first column "-" is the gene identifier, and it strips version numbers.
    """
    with gzip.open(filepath, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
    
    # Rename the gene column and strip version numbers (e.g., ENSG00000123456.1 -> ENSG00000123456)
    df = df.rename(columns={"-": "gene"})
    df['gene'] = df['gene'].str.replace(r'\.\d+$', '', regex=True)
    
    return df

def combine_expression_matrices(files):
    """
    Combines multiple expression matrix files. It aligns on the gene column
    and fills missing values with 0.
    """
    combined_df = None

    for file in files:
        df = read_expression_file(file)

        if combined_df is None:
            combined_df = df
        else:
            combined_df = pd.merge(combined_df, df, on="gene", how="outer")

    # Fill missing values with 0 after merging
    combined_df = combined_df.fillna(0)

    return combined_df

def save_combined_matrix(df, output_file):
    """
    Saves the combined expression matrix to a gzipped text file.
    """
    with gzip.open(output_file, 'wt') as f:
        df.to_csv(f, sep="\t", index=False)

# List of gzipped expression matrix files to combine
files = ['/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/GTEx/tweaked_files/combined_corrected_exp.txt.gz', 
         '/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/Werna/rna/qc/output/9_covariate_correction/1000IBD_gene_counts-TMM.SampleSelection.ProbesWithZeroVarianceRemoved.Log2Transformed.forcenormal.covariatecorrected.txt.gz.CovariatesRemovedOLS.txt.gz', 
         '/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/pub_rna/final_files_pub_rna/merged_expression_data.txt.gz'
         ]

# Combine the expression matrices
combined_df = combine_expression_matrices(files)

# Save the result
output_file = '/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined_expression_matrix.txt.gz'
save_combined_matrix(combined_df, output_file)

print(f"Combined matrix saved to {output_file}")
