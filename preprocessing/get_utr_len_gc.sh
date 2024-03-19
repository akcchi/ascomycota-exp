#!/bin/bash

# Script uses the bioawk extension to calculate lengths and GC content of  
# 5' and 3' untranslated regions

# Each species gets four output files: 
# a list of lengths and a list of GC content, for both 5' and 3' untranslated regions

# Example output (length):
# 	transcript:KAF2020338   314
# 	...
# Example output (GC):
# 	transcript:KAF2020338   0,496815
# 	...

# Calculate lengths
for FILE in UTR5/*; do
	out_file=`basename ${FILE}| sed "s/$/_len/"`
	bioawk -c fastx '{ print $name, length($seq) }' ${FILE} > ~/project_data/inputs/lengths/UTR5_len/${out_file}
done

for FILE in UTR3/*; do
	out_file=`basename ${FILE}| sed "s/$/_len/"`
	bioawk -c fastx '{ print $name, length($seq) }' ${FILE} > ~/project_data/inputs/lengths/UTR3_len/${out_file}
done

# Calculate GC content
for FILE in UTR5/*; do
	out_file=`basename ${FILE}| sed "s/$/_GC/"`
	bioawk -c fastx '{ print $name, gc($seq) }' ${FILE} > ~/project_data/inputs/GC/UTR5_GC/${out_file}
done

for FILE in UTR3/*; do
	out_file=`basename ${FILE}| sed "s/$/_GC/"`
	bioawk -c fastx '{ print $name, gc($seq) }' ${FILE} > ~/project_data/inputs/GC/UTR3_GC/${out_file}
done
