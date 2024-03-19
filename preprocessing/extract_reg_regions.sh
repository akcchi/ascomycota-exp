#!/bin/bash

# Script uses the AGAT toolkit to extract features from a (reference genome) FASTA file,
# using the annotations from a corresponding GFF file

# For each genome, four regulatory regions will be extracted from each gene:
# 	PROMOTER (1000 bp upstream of 5' UTR)
# 	5' UTR
# 	3' UTR
# 	TERMINATOR (500 bp downstream of 3' UTR)
# Each feature has its own directory, where FASTA files will be output

# Assumes two existing directories "gff" and "fasta_dna":
# 	"gff" containing .gff3 files; example name
# 		Aaosphaeria_arxii_cbs_175_79_gca_010015735.Aaoar1.56.gff3
# 	"fasta_dna" containing .fa files; example name
# 		Aaosphaeria_arxii_cbs_175_79_gca_010015735.Aaoar1.dna.toplevel.fa
# NB Each species should have one .gff and .fa file
# NB Files for the same species should have the same name apart from the suffix

# Acknowledgements: script uses code snippets shared by A. Zelezniak

feature="PROMOTER"
mkdir -p $feature
for FILE in ~/project_data/gff/*; do
	dna_file=`basename $FILE| sed "s/.56.gff3$/.dna.toplevel.fa/"`
	agat_sp_extract_sequences.pl -g $FILE -f ~/project_data/fasta_dna/$dna_file -t gene --eo --up 1000 -o $feature/${dna_file}.${feature}
done

feature="UTR5"
mkdir -p $feature
for FILE in gff/*; do
	dna_file=`basename $FILE| sed "s/.56.gff3$/.dna.toplevel.fa/"`
	agat_sp_extract_sequences.pl -g $FILE -f fasta_dna/$dna_file -t five_prime_UTR --full -o $feature/${dna_file}.${feature}
done

feature="UTR3"
mkdir -p $feature
for FILE in gff/*; do
	dna_file=`basename $FILE| sed "s/.56.gff3$/.dna.toplevel.fa/"`
	agat_sp_extract_sequences.pl -g $FILE -f fasta_dna/$dna_file -t three_prime_UTR --full -o $feature/${dna_file}.${feature}
done

feature="TERMINATOR"
mkdir -p $feature
for FILE in ~/project_data/gff/*; do
	dna_file=`basename $FILE| sed "s/.56.gff3$/.dna.toplevel.fa/"`
	agat_sp_extract_sequences.pl -g $FILE -f ~/project_data/fasta_dna/$dna_file -t gene --eo --down 500 -o $feature/${dna_file}.${feature}
done
