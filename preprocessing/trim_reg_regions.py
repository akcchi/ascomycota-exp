# Script truncates sequences FASTA files to a specified length,
# by removing either upstream or downstream sequence

# IMPORTANT: define input and output directories

from Bio import SeqIO # Biopython
import os

# Trim (if necessary) sequences in FASTA files to a target length, 
# by removing UPSTREAM sequence
# i.e. use for promoter (1000 bp) and 5' UTR (300 bp)
def remove_upstream(input, output,target):
    target_length = int(target)

    # Read and parse FASTA file
    working_sequence = SeqIO.parse(open(input), "fasta") # filepath, format = "fasta"

    with open(output, "w") as out_file:
        for fasta in working_sequence:
            header = ">" + fasta.id
            sequence = str(fasta.seq)
            
            # Remove upstream if longer than target; 
            # by keeping last (target length) nucleotides
            if len(sequence) > target_length:
                sequence = sequence[-target_length:]

            # Write to output file
            out_file.write(header + "\n")
            out_file.write(sequence + "\n")

# Trim (if necessary) sequences in FASTA files to a target length, 
# by removing DOWNSTREAM sequence
# i.e. use for 3' UTR (350 bp) and terminator (500 bp)
def remove_downstream(input, output, target):
    target_length = int(target)

    # Read and parse FASTA file
    working_sequence = SeqIO.parse(open(input), "fasta") # filepath, format = "fasta"

    with open(output, "w") as out_file:
        for fasta in working_sequence:
            header = ">" + fasta.id
            sequence = str(fasta.seq)
            
            # Remove downstream if longer than target; 
            # by keeping first (target length) nucleotides
            if len(sequence) > target_length:
                sequence = sequence[:target_length]

            # Write to output file
            out_file.write(header + "\n")
            out_file.write(sequence + "\n")


# Trim promoter sequences to 1000 bp
working_directory = "./project_data/inputs/padded/PROMOTER_single_pad1000"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/inputs/padded_cut/PROMOTER_final"
    input_file = os.path.join(working_directory, filename)
    output_file = os.path.join(output_directory, filename + ".cut")

    remove_upstream(input_file, output_file, 1000)
    print("Trimmed promoters to 1000 bp: " + filename)

# Trim 5' UTR sequences to 300 bp
working_directory = "./project_data/inputs/padded/UTR5_single_pad300"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/inputs/padded_cut/UTR5_final"
    input_file = os.path.join(working_directory, filename)
    output_file = os.path.join(output_directory, filename + ".cut")

    remove_upstream(input_file, output_file, 300)
    print("Trimmed 5' UTR sequences to 300 bp: " + filename)

# Trim 3' UTR sequences to 350 bp
working_directory = "./project_data/inputs/padded/UTR3_single_pad350"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/inputs/padded_cut/UTR3_final"
    input_file = os.path.join(working_directory, filename)
    output_file = os.path.join(output_directory, filename + ".cut")

    remove_downstream(input_file, output_file, 350)
    print("Trimmed 3' UTR sequences to 350 bp: " + filename)

# Trim terminator sequences to 500 bp
working_directory = "./project_data/inputs/padded/TERMINATOR_single_pad500"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/inputs/padded_cut/TERMINATOR_final"
    input_file = os.path.join(working_directory, filename)
    output_file = os.path.join(output_directory, filename + ".cut")

    remove_downstream(input_file, output_file, 500)
    print("Trimmed terminators to 500 bp: " + filename)
