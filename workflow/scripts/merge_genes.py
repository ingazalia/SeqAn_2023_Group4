from Bio import SeqIO
from os import listdir, makedirs
from os.path import isfile, join

fasta_files = [f for f in listdir("results/cutgeneseqs") if isfile(join("results/cutgeneseqs", f))]

merged_sequences = {}

# Iterate through the FASTA files
for fasta_file in fasta_files:
    # Read sequences from the current FASTA file
    for record in SeqIO.parse("results/cutgeneseqs/"+fasta_file, "fasta"):
        gene_id = record.id
        sequence = str(record.seq)
        if sequence.count('N') / len(sequence) < float(snakemake.params["threshold"]):
            gene_id = gene_id.replace(')', '_')
            gene_id = gene_id.replace('(', '_')
            # Check if the gene ID already exists in the merged_sequences dictionary
            if gene_id in merged_sequences:
                # Concatenate the sequences
                merged_sequences[gene_id] += ">" + fasta_file + "\n" + sequence + "\n"
            else:
                # Add the gene ID and sequence to the dictionary
                merged_sequences[gene_id] = ">" + fasta_file + "\n" + sequence + "\n"

gene_names = []

# Iterate over the keys and save them to the list
for gene_name in merged_sequences.keys():
    count = merged_sequences[gene_name].count(">")
    if count == len(fasta_files):
      gene_names.append(gene_name)

makedirs("results/genes", exist_ok=True)

# Iterate over the gene names and write the sequences to separate FASTA files

for gene_name in gene_names:
    with open("results/genes/" + gene_name + ".fasta", "w") as output_file:
        for gene_id, sequence in merged_sequences.items():
            if gene_id == gene_name:
                output_file.write(f"{sequence}\n")

