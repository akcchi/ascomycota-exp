#!/bin/bash

# Script compares FASTA headers across promoter, 5'UTR, 3'UTR, and terminator files

# For each species, return list of transcript names with a complete set 
# of all four regulatory regions

for FILE in ~/project_data/inputs/padded_cut/PROMOTER_final/*; do
    species=`basename $FILE| sed "s/.dna.toplevel.fa.PROMOTER.final//"`
    outfile=`basename $FILE| sed "s/.dna.toplevel.fa.PROMOTER.final/.transcripts_completeset.txt/"`
    
    # Compare FASTA headers across the four files
    comm -12 <(grep ">" ${FILE} | sed 's/>transcript://' | sort) <(grep ">" ~/project_data/inputs/padded_cut/UTR5_final/${species}* | sed 's/>transcript://' | sort) \
    | comm -12 - <(grep ">" ~/project_data/inputs/padded_cut/UTR3_final/${species}* | sed 's/>transcript://' | sort) \
    | comm -12 - <(grep ">" ~/project_data/inputs/padded_cut/TERMINATOR_final/${species}* | sed 's/>transcript://' | sort) \
    > ~/project_data/inputs/temp/transcripts_completeset_index/${outfile}
    
    echo "Generated: ${outfile}"
done