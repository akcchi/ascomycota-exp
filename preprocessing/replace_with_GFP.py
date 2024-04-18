# Replace the CDS variables for the yeast/Sacce1 test dataset with that of EGFP:
## Length of CDS
## GC content at codon positions 1, 2, 3
## Codon usage frequencies

# Yeast test dataset previously extracted from Zenodo repository published by Jan Zrimec
# (https://doi.org/10.5281/zenodo.3905252)
# filename: scerevisiae.rsd1.lmbda_22.npz; arr_1 (onehot) and arr_3 (72 variables)

import numpy as np


length_EGFP = int(720) # Previously calculated, hardcoded, sorry
EGFP_GC_file = "./project_data/inputs/EGFP/EGFP_codon_GC.npy"
EGFP_codon_file = "./project_data/inputs/EGFP/EGFP_codon_usage.npy"

# Load files
EGFP_GC_arr = np.load(EGFP_GC_file)
EGFP_codon_arr = np.load(EGFP_codon_file)


yeast_npz_filepath = "./Saccharomyces_cerevisiae_R64-1-1.Sacce1.npz"
yeast_og = np.load(yeast_npz_filepath, allow_pickle=True)

print(yeast_og.files)

og_hot = yeast_og["hot"] # no need to modify this
og_var = yeast_og["var"] # modify this so CDS variables are that of EGFP


for i in range(len(og_var)):
    # Replace CDS length
    (og_var[i])[1] = length_EGFP

    # Replace GC content for codon positions 1, 2, 3
    for index in range(5, 8):
        replacement = round(float(EGFP_GC_arr[index - 5])*1000)
        (og_var[i])[index] = replacement

    # Replace codon usage frequencies
    for index in range(8, 72):
        replacement = EGFP_codon_arr[index - 8]
        (og_var[i])[index] = replacement


outfile = "./Saccharomyces_cerevisiae_R64-1-1.Sacce1_MODIFIED.npz"
np.savez(outfile, og_hot=og_hot, og_var=og_var)
print("Done!")
