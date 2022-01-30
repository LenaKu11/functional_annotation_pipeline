#!/usr/bin/env python3
"""

python3 identify_genes_internal_stop_codons.py file.fasta genes.txt
"""

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
gene_file = sys.argv[2]  # Output file, one gene name per line

with open(fasta_file, "r") as fa, open(gene_file, "w") as li:
    for record in SeqIO.parse(fa, "fasta"):
        if record.seq.count('.') != 0:
            genename = record.id
            genename = genename[:-2]
            li.write(genename + '\n')

fa.close()
li.close()
