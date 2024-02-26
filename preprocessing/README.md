Mixture of Bash and Python scripts to preprocess nucleotide sequence data.

## Dependencies
TODO: add conda environment file

## Usage and order
It is assumed that FASTA files and corresponding GFF files have already been retrieved for each genome.<br>
Rename file/directory names in scripts as appropriate.

### extract_reg_regions.sh
For every gene in each genome, extract the four following regulatory regions: 'promoter region' (1000 bp upstream of 5' untranslated region), 5' untranslated region, 3' untranslated region, and 'terminator region' (500 bp downstream of 3' untranslated region). For each genome, outputs four FASTA files (one for each regulatory region). This script uses the [AGAT toolkit](https://github.com/NBISweden/AGAT).

### fasta_single_line.py

### get_utr_len_gc.sh

### pad_reg_regions.py

### trim_reg_regions.py

### get_gene_transcript_names.sh

### gene_to_transcript_names.py

### get_reg_set.sh

### get_gfp_variables.py

### assemble_npz.py



