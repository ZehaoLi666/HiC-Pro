#!/usr/bin/env python3

import sys
from Bio import SeqIO

g = SeqIO.parse(open(sys.argv[1]),"fasta")
out = open((sys.argv[1]).replace("fasta","chrom.sizes"),"w")
for i in g:
	out.write(i.id + "\t" + str(len(i)) + "\n")
out.close()