Mixture of Bash and Python scripts to preprocess the nucleotide sequences of fungi (Ascomycota) genomes into input datasets for use with the _Saccharomyces cerevisiae_ expression neural network model published by [Zrimec et al.](https://doi.org/10.1038/s41467-020-19921-4)

## Dependencies
See conda [environment.yml](/preprocessing/environment.yml) file.

## Usage and order
It is assumed that FASTA files and corresponding GFF files have already been retrieved for each genome.<br>
Rename file/directory names in scripts as appropriate.

### extract_reg_regions.sh
Where possible for every gene in each genome, extract the four following regulatory regions: 'promoter region' (1000 bp upstream of 5' untranslated region), 5' untranslated region, 3' untranslated region, and 'terminator region' (500 bp downstream of 3' untranslated region). For each genome, outputs four FASTA files (one for each regulatory region). This script uses the [AGAT toolkit](https://github.com/NBISweden/AGAT). 

### fasta_single_line.py
For all FASTA files, convert all multi-line sequences into single-line sequences. This script uses the [bioinfokit toolkit](https://github.com/reneshbedre/bioinfokit). 

### get_utr_len_gc.sh
For each genome, calculate the lengths and GC content of the 5' and 3' untranslated regions. For each genome, outputs four plaintext files (two variables, for both untranslated regions). This script uses the [bioawk extension](https://github.com/lh3/bioawk). 

### pad_reg_regions.py
Zero-pad (where necessary) sequences to the specified length for the type of regulatory region: 1000 bp (promoter region), 300 bp (5' untranslated region), 350 bp (3' untranslated region), or 500 bp (terminator region). 

### trim_reg_regions.py
Truncate (where necessary) sequences to the specified length for the type of regulatory region: 1000 bp (promoter region), 300 bp (5' untranslated region), 350 bp (3' untranslated region), or 500 bp (terminator region). Upstream sequence is removed from the promoter and 5' untranslated regions, while downstream sequence is removed from the terminator and 3' untranslated regions. 

### get_gene_transcript_names.sh
For each genome, find the corresponding transcript name for every given gene name. This is needed as the FASTA files containing promoter or terminator regions only contain the gene name in headers, not the transcript name (see [extract_reg_regions.sh](#extract_reg_regionssh) outputs). For each genome, outputs one plaintext file for each of the promoter and terminator regions.

### gene_to_transcript_names.py
For all of the promoter and terminator region FASTA files, use the [previously generated files](#get_gene_transcript_namessh) to replace all gene names in FASTA headers with their corresponding transcript names.

### get_reg_set.sh
For each genome, get a list of transcript names with a complete set of the four regulatory regions (not all transcripts have these regions annotated in the GFF file): promoter region, 5' untranslated region, 3' untranslated region, and terminator region.

### get_gfp_variables.py
Calculate the following variables for the coding sequence (CDS) of enhanced GFP (sequence supplied via a FASTA file): length, GC content at codon positions 1/2/3, and codon usage frequencies. This script uses some code snippets by Zrimec et al.; see code comments for specific acknowledgments. 

### assemble_npz.py
For each genome, generate a .npz file containing an input dataset in the required format/dimensions for the deep learning model. Each dataset includes: one-hot encoded 2150 bp regulatory sequence, as well as an additional 72 variables (as outlined in [Zrimec et al.](https://doi.org/10.1038/s41467-020-19921-4); see Methods section). This script uses some code snippets by Zrimec et al.; see code comments for specific acknowledgments.
