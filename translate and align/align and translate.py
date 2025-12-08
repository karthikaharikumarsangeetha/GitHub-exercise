# importing required packages from biopython

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from scipy.stats import linregress
from collections import defaultdict
import os

# File containing all FASTA files
file_path = "C:/Users/pj25643/Documents/proteins/GitHub-exercise/translate and align/translate and align"
# combined DNA sequences of the samples
input_file = r"C:/Users/pj25643/Documents/proteins/GitHub-exercise/translate and align/combined_sequences.fasta"
print(f"reading combined DNA sequences from: {input_file}")
# Output File for translated sequences
output_file = r"C:/Users/pj25643/Documents/proteins/GitHub-exercise/translate and align/translated_sequence.fasta"
print(f" Translated proteins will be saved to: {output_file}")

# combine fragements per sample
sample_sequences = defaultdict(str)

#
# Parse the input and combine the fragments
for record in SeqIO.parse(input_file, "fasta"):
    sample_name = record.id.split('part')[0]  # for extracting sample names
    # concatenate the fragments
    sample_sequences[sample_name] += str(record.seq)

print(f"Total unique samples found: {len(sample_sequences)}")

# storing translated AA sequences from ORF
aa_sequence = []


# loop through each sequence grouped by samples

for idx, (sample_name, dna_seq) in enumerate(sample_sequences.items(), start=1):
    print(f"/nProcessing sample {idx}/{len(sample_sequences)}:{sample_name}")

# convert string DNA to seq object
    seq_obj = Seq(dna_seq)
    best_orf = ""
    best_frame = None
# Translate the reading frames
    for frame in range(3):
        protein = seq_obj[frame:].translate(to_stop=False)
# split the protein sequence to orf
        fragments = str(protein).split("*")
        longest_fragment = max(fragments, key=len)
        print(
            f" Frame {frame}: longest fragment length = {len(longest_fragment)}")

        if len(longest_fragment) > len(best_orf):
            best_orf = longest_fragment
            best_frame = frame

# print the frame and ORF
    print(f" selected frame: {best_frame}, ORF length: {len(best_orf)}")

    aa_sequence.append(
        SeqRecord(Seq(best_orf), id=sample_name,
                  description=f"Longest ORF in frame{best_frame}")
    )

print(f"\nsaving{len(aa_sequence)} translated protein sequences...")

# saving translated sequences to a FASTA File
SeqIO.write(aa_sequence, output_file, "fasta")
print(f"\nALL sequences sved to:\n{output_file}\n")
