# Script replaces all FASTA header gene names with transcript names
# As AGAT tool only generates FASTAs with gene names for intergenic regions 
# (promoters and terminators)
# Uses previously made reference files

import os
import re

# Replace gene names in these FASTA files with transcript names
# using previously made reference files
# Example line:
# KAF2020338      BU24DRAFT_416053
def gene_to_transcript(reference_file, input_fasta, output_fasta):
    # Read reference file and assign contents to dictionary
    transcripts_dict = {}
    with open(reference_file, "r") as file:
        for line in file:
            gene_transcript = line.strip().split("\t")
            if len(gene_transcript) == 2:
                transcript, gene = gene_transcript 
                # gene ID = key; transcript name = value
                transcripts_dict[gene] = transcript

    # Read whole FASTA to memory (it's amateur hour)
    with open(input_fasta, "r") as input_file:
        contents = input_file.read()

    for gene_id, transcript_id in transcripts_dict.items():
        # Replace gene name with transcript name
        contents = re.sub(re.escape(gene_id), transcript_id, contents)
            # arguments: replace all of "genes" with "transcript" in "contents"
            # re.escape so "gene" is treated as literal string
        
        # Replace "gene" in fasta header with "transcript"
        contents = re.sub(re.escape("gene"), "transcript", contents)

    # Write to output file
    with open(output_fasta, "w") as output_file:
        output_file.write(contents)

# Get PROMOTER fastas with transcript names
reference_directory = "./project_data/inputs/temp/PROMOTER_FINAL-transcript-gene"
fasta_directory = "./project_data/inputs/padded/PROMOTER_single_pad1000"
output_directory = "./project_data/inputs/padded_cut/PROMOTER_final"
for filename in os.listdir(reference_directory):
    ref = os.path.join(reference_directory, filename)
        # Example: Yarrowia_lipolytica_gca_003367965.YarliYB420_FINAL-transcript-gene.txt
    
    in_fasta = os.path.join(fasta_directory, re.sub(re.escape("_FINAL-transcript-gene.txt"), ".dna.toplevel.fa.PROMOTER.single.pad1000.cut", filename))
        # Example: Yarrowia_lipolytica_gca_003367965.YarliYB420.dna.toplevel.fa.PROMOTER.single.pad1000.cut

    out_fasta = os.path.join(output_directory, re.sub(re.escape("_FINAL-transcript-gene.txt"), ".dna.toplevel.fa.PROMOTER.final", filename))
        # Example: Yarrowia_lipolytica_gca_003367965.YarliYB420.dna.toplevel.fa.PROMOTER.final

    gene_to_transcript(ref, in_fasta, out_fasta)
    print("Promoter file processed: " + out_fasta)

# Get TERMINATOR fastas with transcript names
reference_directory = "./project_data/inputs/temp/TERMINATOR_FINAL-transcript-gene"
fasta_directory = "./project_data/inputs/padded/TERMINATOR_single_pad500"
output_directory = "./project_data/inputs/padded_cut/TERMINATOR_final"
for filename in os.listdir(reference_directory):
    ref = os.path.join(reference_directory, filename)
        # Example: Yarrowia_lipolytica_gca_003367965.YarliYB420_FINAL-transcript-gene.txt
    
    in_fasta = os.path.join(fasta_directory, re.sub(re.escape("_FINAL-transcript-gene.txt"), ".dna.toplevel.fa.TERMINATOR.single.pad500.cut", filename))
        # Example: Yarrowia_lipolytica_gca_003367965.YarliYB420.dna.toplevel.fa.TERMINATOR.single.pad500.cut

    out_fasta = os.path.join(output_directory, re.sub(re.escape("_FINAL-transcript-gene.txt"), ".dna.toplevel.fa.TERMINATOR.final", filename))
        # Example: Yarrowia_lipolytica_gca_003367965.YarliYB420.dna.toplevel.fa.TERMINATOR.final

    gene_to_transcript(ref, in_fasta, out_fasta)
    print("Terminator file processed: " + out_fasta)
