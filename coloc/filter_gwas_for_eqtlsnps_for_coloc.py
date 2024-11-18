# Author: Harm-Jan Westra, modified by: Peter Riesebos
# Purpose: filter GWAS sum stats on overlapping eQTL rsids, per chromosome
# input: filtered eQTL sum stat files (per chromosome), GWAS sum stats
# output: filtered GWAS variants per chromosome, tsv files
 
import gzip


fh = open('/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/combined/mbqtl_output_combined_exp_fixed/combined-AllEffects_200kb_significant-snps.txt','r')
snps = fh.readlines()
fh.close()

snpset = set()
for snp in snps:
	snpset.add(snp.strip())

fh = gzip.open('/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/GWAS/GCST90292538.h.tsv.gz','rt')


fhoarr=[None] * 24
header=fh.readline()
for i in range(1,23):
	fhoarr[i] = gzip.open(f"/groups/umcg-fg/tmp04/projects/gut-bulk/ongoing/2024-02-07-GutPublicRNASeq/datasets/GWAS/GCST90292538.h_filtered_eQTLSNPs_200kb-chr{i}.tsv.gz",'wt',3)
	fhoarr[i].write(header)

wctr = 0
lctr = 0

for line in fh:
	elems = line.strip().split("\t")
	rsid = elems[9]
	chr = 0
	try:
		chr = int(elems[0])
	except:
		chr = 0
	if chr > 0 and chr < 23:
		if rsid in snpset:
			fhoarr[chr].write(line)
			wctr += 1
	lctr += 1
	print(f"{lctr} / {wctr}",end='\r')
print(f"{lctr} / {wctr}",end='\n')

for i in range(1,23):
	fhoarr[i].close()
fh.close()

print(f"{wctr} written")
