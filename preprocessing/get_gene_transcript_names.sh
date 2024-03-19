#!/bin/bash

# Script finds the corresponding transcript name for every given gene name,
# for each genome; by searching through the relevant GFF file for that genome
# Writes the matching transcript names and gene names to output file

# For each genome, extract list of gene names the promoter sequences belong to
for FILE in ~/project_data/inputs/padded_cut/PROMOTER_final/*; do
    out_file=`basename $FILE| sed "s/$/_genes.txt/"`
    grep ">" ${FILE} | sed "s/>gene://" > ~/project_data/inputs/temp/PROMOTER_gene-names/${out_file}
done

# For each genome, extract list of gene names the terminator sequences belong to
for FILE in ~/project_data/inputs/padded_cut/TERMINATOR_final/*; do
    out_file=`basename $FILE| sed "s/$/_genes.txt/"`
    grep ">" ${FILE} | sed "s/>gene://" > ~/project_data/inputs/temp/TERMINATOR_gene-names/${out_file}
done

# Now for each genome, search corresponding GFF file for transcript names
# that match the list of gene names
# Output the matching transcript name and gene name to file

# Example line in output file:
# KAF2020338      BU24DRAFT_416053

# Do for promoter sequences
for FILE in ~/project_data/inputs/temp/PROMOTER_gene-names/*; do
    gff_file=`basename $FILE| sed "s/dna.toplevel.fa.PROMOTER.single.pad1000.cut_genes.txt/56.gff3/"`
    out_file=`basename $FILE| sed "s/dna.toplevel.fa.PROMOTER.single.pad1000.cut_genes.txt/transcript-gene.txt/"`
    
    # grep '-f' flag to search GFF file using the list of gene names
    grep -f ${FILE} ~/project_data/gff/${gff_file} | grep "transcript" | cut -f 9 | cut -d ";" -f 1,2 \
    | sed "s/ID=transcript://" | sed "s/;Parent=gene:/\t/" | grep -v "ID" \
    > ~/project_data/inputs/temp/PROMOTER_transcript-gene/${out_file}
done

# Now repeat for terminator sequences
for FILE in ~/project_data/inputs/temp/TERMINATOR_gene-names/*; do
    gff_file=`basename $FILE| sed "s/dna.toplevel.fa.TERMINATOR.single.pad500.cut_genes.txt/56.gff3/"`
    out_file=`basename $FILE| sed "s/dna.toplevel.fa.TERMINATOR.single.pad500.cut_genes.txt/transcript-gene.txt/"`
    grep -f ${FILE} ~/project_data/gff/${gff_file} | grep "transcript" | cut -f 9 | cut -d ";" -f 1,2 \
    | sed "s/ID=transcript://" | sed "s/;Parent=gene:/\t/" | grep -v "ID" \
    > ~/project_data/inputs/temp/TERMINATOR_transcript-gene/${out_file}
done