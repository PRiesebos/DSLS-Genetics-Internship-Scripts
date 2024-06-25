import pandas as pd
import gzip

# Read the alignment metrics file
alignment_metrics_file = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/alignment_metrics.txt.gz"
df = pd.read_csv(alignment_metrics_file, sep='\t')

# Read the study_samples_mapping.txt and create a dictionary mapping samples to studies
study_samples_mapping = {}
with open('/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/output/study_samples_mapping.txt', 'r') as f:
    for line in f:
        parts = line.strip().split()
        study = parts[0]
        samples = parts[1:]
        for sample in samples:
            study_samples_mapping[sample] = study

# Map the sample names to study names
df['Study'] = df['Sample'].map(study_samples_mapping)

# Reorder the columns
df = df[['Study', 'Sample', 'PCT_CODING_BASES', 'PCT_PF_READS_ALIGNED']]

# Save the final dataframe to a new file
final_output_file = "/scratch/hb-functionalgenomics/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/5-expressieQC/all_alignment_output/alignment_metrics_with_study.txt.gz"
with gzip.open(final_output_file, 'wt') as f:
    df.to_csv(f, index=False, sep='\t')

print("Data processing complete. The final dataframe has been saved to", final_output_file)
