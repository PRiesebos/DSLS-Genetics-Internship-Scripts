# author: Harm-Jan Westra, modified by: Peter Riesebos
# purpose: coloc prep step, filter eQTLs for coloc. By taking a 200kb window around the TSS we only keep variants 
# from the AllEffects files where we have overlapping genes and SNPPos, 
# taking the variants where the SNPPos is within the 200kb window of the TSS
# input: eQTL top effect file, eQTl AllEffect files
# output: filtered eQTL sum stats for each chromosome, eQTL snp ids for each chromosome

import gzip
import sys

chr = sys.argv[1]
fh = open('/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/mbqtl_output_combined_exp_fixed/merged_topeffects_final.txt','r')

fh.readline()
signgenes = set()
for line in fh:
	elems = line.strip().split("\t")
	signgenes.add(elems[0])
fh.close()

maxdist = 200000


for chr in [chr]:
	validsnps = set()

	file = f"/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/mbqtl_output_combined_exp_fixed/combined_chr{chr}-AllEffects.txt.gz"
	fileout = f"/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/mbqtl_output_combined_exp_fixed/combined_chr{chr}-AllEffects_200kb_significant.txt.gz"
	print(file)
	print(fileout)

	lctr = 0
	wctr = 0
	fh = gzip.open(file,'rt')
	fho = gzip.open(fileout,'wt')
	fho.write(fh.readline())
	for line in fh:
		elems = line.strip().split("\t")
		gene = elems[0]
		if gene in signgenes:
			snp = elems[5]
			snppos = int(elems[7])
			genepos = int(elems[2])
			if abs(snppos - genepos) <= maxdist:
				validsnps.add(snp)
				fho.write(line)
				wctr += 1

		lctr += 1
		if lctr % 10000 == 0:
			print(f"{lctr} lines, {wctr} written",end='\r')
	print(f"{lctr} lines, {wctr} written",end='\n')
	fho.close()
	fh.close()

	fileout = f"/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/mbqtl_output_combined_exp_fixed/combined_chr{chr}-AllEffects_200kb_significant-snps.txt.gz"
	fho = open(fileout,'wt')
	for snp in validsnps:
		fho.write(snp+"\n")
	fho.close()
