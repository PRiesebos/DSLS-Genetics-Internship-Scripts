# Author: Harm-Jan Westra, Modified by: Peter Riesebos

import gzip

fh = open('path/to/sample/names/file.txt','r')
samples = set()
for line in fh:
    samples.add(line.strip())
fh.close()

fh = gzip.open('path/to/expression/file.gz','rt')
fho = gzip.open('path/to/output/expression/file.txt.gz','wt')
header  = fh.readline().strip().split("\t")
colsToInclude = []
outheader = 'Sample'
for i in range(0,len(header)):
    sample = header[i]
    if sample in samples:
        colsToInclude.append(i)
        outheader += "\t"+sample
fho.write(outheader+"\n")
for line in fh:
    elems = line.strip().split("\t")
    outln = elems[0]
    vals = []
    for col in colsToInclude:
        vals.append(elems[col])
    outln = outln +"\t"+ "\t".join(vals)
    fho.write(outln+"\n")
fh.close()
fho.close()