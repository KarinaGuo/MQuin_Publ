import sys
from Bio import AlignIO

filename = sys.argv[1]


alignment = AlignIO.read(open(filename), "fasta")
print("Alignment length %i" % alignment.get_alignment_length())
for record in alignment:
    if (record.seq.count('N') / alignment.get_alignment_length() ) < 0.4:
       print(">" + record.id + "\n" + record.seq + "\n")

