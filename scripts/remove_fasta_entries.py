#!/usr/bin/env python3
"""

python3 remove_fasta_entries.py file.fasta genes.txt out.fasta
"""

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file
gene_file = sys.argv[2]  # File listing fasta record ids, one gene per line (without ">")
filtered_fasta = sys.argv[3] # Output fasta file, records from list removed

with open(fasta_file, "r") as fa, open(gene_file, "r") as li, open(filtered_fasta, "w") as fi:
    notwanted = [r.strip() for r in li]
    for record in SeqIO.parse(fa, "fasta"):
        if record.id in notwanted:
            continue
        else:
            fi.write(record.format("fasta"))

fa.close()
li.close()
fi.close()
