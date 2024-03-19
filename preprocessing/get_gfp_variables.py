# Script calculates the following variables for the EGFP coding sequence:
## Length
## GC content at codon positions 1, 2, and 3
## Codon usage frequencies

# Acknowledgements: script uses some functions shared and authored by Zrimec et al.
# See comments for specifics

from Bio import SeqIO
import numpy as np
import os

# Function shared and authored by Zrimec et al.
# Calculate GC content at codon positions 1/2/3
# Takes DNA sequence as input "file"
def count_codon_GC(file):
    out1=list()
    out2=list()
    out3=list()
    for i in range(0,len(file)-2,3):
        out1.append(int(file[i] in ['G','C']))
        out2.append(int(file[i+1] in ['G','C']))
        out3.append(int(file[i+2] in ['G','C']))
    return 3*sum(out1)/len(file),3*sum(out2)/len(file),3*sum(out3)/len(file)

# Function shared and authored by Zrimec et al.
# Calculate codon usage frequencies
# Takes DNA sequence as input "file"
def count_codons(file):
    '''codon frequency counter'''
    CodonsDict = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}
    
    # make the codon dictionary local
    codon_count = CodonsDict.copy()
    # iterate over sequence and count all the codons in the string.
    # make sure the sequence is upper case
    if str(file).islower():
        dna_sequence = str(file).upper()
    else:
        dna_sequence = str(file)
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            #raise TypeError("illegal codon %s" % (codon))
            print("illegal codon %s" % (codon))
    # return values in dict with sorted keys alphabetically
    out=list()
    for key,value in sorted(codon_count.items()):
        out.append(value)
    
    return np.asarray(out)


egfp_sequence = SeqIO.parse("EGFP.fasta", "fasta")
    # Create SeqRecord object from EGFP.fasta, specify "fasta" format
outputs = "./project_data/inputs/EGFP"

for sequence in egfp_sequence:
    # Length of EGFP CDS
    length = len(sequence.seq)

    # GC content at codon positions 1/2/3
    codon_GC = count_codon_GC(sequence.seq)
    codon_GC_arr = np.array(codon_GC)
    # Save to .npy file
    codon_GC_filepath = os.path.join(outputs, "EGFP_codon_GC.npy")
    np.save(codon_GC_filepath, codon_GC_arr)

    # Codon usage frequencies
    codon_usage_file = os.path.join(outputs, "EGFP_codon_usage.npy")
    # Save to .npy file
    np.save(codon_usage_file, count_codons(sequence.seq))

    # Write length and codon GC to another file
    with open(os.path.join(outputs, "EGFP_variables.txt" ), "w") as output_file:
        output_file.write(str(length) + "\n" + str(codon_GC))

    print("Done")