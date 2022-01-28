from Bio import SeqIO
import sys

max_len = 0
max_description = ""

for record in SeqIO.parse(sys.argv[1], "fasta"):
    if len(record) > max_len:
        max_len = len(record)
        max_description = record.description

print(max_description)
print(max_len)
