# Script pads promoter/5'UTR/3'UTR/terminator sequences in FASTA files 
# to their repsective target lengths (1000/300/350/500 bp respectively)

# IMPORTANT: define input and output directories

import os

# Read input FASTA file sequences, pad if necessary to specified target length, 
# output new FASTA file
def pad_sequence(in_file, out_file, target_length):
    with open(in_file, "r") as file_in, open(out_file, "w") as file_out:
        for line in file_in:
            line = line.strip()
            if line.startswith(">"):
                file_out.write(line + "\n")  # Write header line as is
            else:
                sequence = line
                if len(sequence) <= int(target_length):
                    # Zero pad front of sequence to meet target length
                    padded_sequence = "0" * (int(target_length) - len(sequence)) + sequence
                    # Write padded sequence
                    file_out.write(padded_sequence + "\n")
                else:
                    # Keep sequence as is, if already longer than target length
                    file_out.write(sequence + "\n")


# Directory containing promoter single line FASTA files
working_directory = "./project_data/PROMOTER_single"
for filename in os.listdir(working_directory):
    # Output directory, edit this
    output_directory = "./project_data/inputs/padded/PROMOTER_single_pad1000"
    
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".pad1000")

    # Target length: 1000 bp
    pad_sequence(input, output, 1000)
    print("Promoter file processed: " + filename)

# Directory containing 5' UTR single line FASTA files
working_directory = "./project_data/UTR5_single"
for filename in os.listdir(working_directory):
    # Output directory, edit this
    output_directory = "./project_data/inputs/padded/UTR5_single_pad300"
    
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".pad300")

    # Target length: 300 bp
    pad_sequence(input, output, 300)
    print("5' UTR file processed: " + filename)

# Directory containing 3' UTR single line FASTA files
working_directory = "./project_data/UTR3_single"
for filename in os.listdir(working_directory):
    # Output directory, edit this
    output_directory = "./project_data/inputs/padded/UTR3_single_pad350"
    
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".pad350")

    # Target length: 350 bp
    pad_sequence(input, output, 350)
    print("3' UTR file processed: " + filename)

# Directory containing terminator single line FASTA files
working_directory = "./project_data/TERMINATOR_single"
for filename in os.listdir(working_directory):
    # Output directory, edit this
    output_directory = "./project_data/inputs/padded/TERMINATOR_single_pad500"
    
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".pad500")

    # Target length: 500 bp
    pad_sequence(input, output, 500)
    print("Terminator file processed: " + filename)
