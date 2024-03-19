# Script converts all multi-line sequences in FASTA files
# to single-line sequences
# Outputs new FASTA files

# IMPORTANT: define input and output directories

import bioinfokit as bioinfokit
from bioinfokit.analys import Fasta
import os
import shutil

# Directory containing promoter sequences
working_directory = "./project_data/PROMOTER"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/PROMOTER_single"
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".single")

    Fasta.multi_to_single_line(file=input)
    
    shutil.move("output.fasta", output)

    print("(Promoter) File processed: " + filename)

# Directory containing terminator sequences
working_directory = "./project_data/TERMINATOR"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/TERMINATOR_single"
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".single")

    Fasta.multi_to_single_line(file=input)
    
    shutil.move("output.fasta", output)

    print("(Terminator) File processed: " + filename)

# Directory containing 5' UTR sequences
working_directory = "./project_data/UTR5"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/UTR5_single"
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".single")

    Fasta.multi_to_single_line(file=input)
    
    shutil.move("output.fasta", output)

    print("(5' UTR) File processed: " + filename)

# Directory containing 3' UTR sequences
working_directory = "./project_data/UTR3"
for filename in os.listdir(working_directory):
    output_directory = "./project_data/UTR3_single"
    input = os.path.join(working_directory, filename)
    output = os.path.join(output_directory, filename + ".single")

    Fasta.multi_to_single_line(file=input)
    
    shutil.move("output.fasta", output)

    print("(3' UTR) File processed: " + filename)
