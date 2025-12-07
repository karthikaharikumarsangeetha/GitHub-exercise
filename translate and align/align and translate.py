# importing required packages from biopython

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from scipy.stats import linregress
from collections import defaultdict
import os

#
file_path = "C:/Users/pj25643/Documents/proteins/GitHub-exercise/translate and align/translate and align"
# folder containing the FASTA files
input_file = r"C:/Users/pj25643/Documents/proteins/GitHub-exercise/translate and align/combined_sequences.fasta"
print(f"reading combined DNA sequences from: {input_file}")
# Output File for translated sequences
output_file = r"C:/Users/pj25643/Documents/proteins/GitHub-exercise/translate and align/translated_sequence.fasta"
print(f" Translated proteins will be saved to: {output_file}")

# combine fragements per sample
sample_sequences = defaultdict(str)

for record in SeqIO.parse(input_file, "fasta"):
    sample_name = record.id.split('part')[0]
    sample_sequences[sample_name] += str(record.seq)

print(f"Total unique samples found: {len(sample_sequences)}")
