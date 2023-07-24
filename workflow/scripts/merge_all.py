import os
from Bio import SeqIO
import pandas as pd

# Directory path containing the FASTA files
fasta_directory = "results/alignment/"
samples_small = pd.read_csv(str(snakemake.input["sample"]), index_col="sample", sep='\t')

sample_small = [s for s in list(samples_small.index)]


# Dictionary to store the merged sequences
merged_sequences = {}
# Iterate through the FASTA files
for file_name in os.listdir(fasta_directory):
    file_path = os.path.join(fasta_directory, file_name)
    if os.path.isfile(file_path) and file_name.endswith(".fasta"):
        # Read sequences from the current FASTA file
        sample_all=[]
        for record in SeqIO.parse(file_path, "fasta"):
                        # Extract the sample name from the header
            sample = record.id
            sequence = str(record.seq)
            sample_all.append(sample)
            # Check if the sample already exists in the merged_sequences dictionary
            if sample in merged_sequences:
             #   # Concatenate the sequences
                merged_sequences[sample] += str(record.seq)
            else:
            #    # Add the sample and sequence to the dictionary
                merged_sequences[sample] = str(record.seq)


with open(str(snakemake.output), "w") as output_file:
    for sample_name, sequence in merged_sequences.items():
        output_file.write(f">{sample_name}\n{sequence}\n")