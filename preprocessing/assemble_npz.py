# Script produces input feature arrays for the neural network model by Zrimec et al.

# RUN LAST: Requires previously generated files

# Acknowledgements: script uses some functions shared and authored by Zrimec et al.
# See comments for specifics

import re
import numpy as np
from keras.utils import to_categorical
from Bio import SeqIO
import glob

# Read text file containing list of transcript names, and
# Save each line into list "transcripts_ID_list"
def get_transcripts_ID_list(ID_file):
    with open (ID_file, "r") as infile:
        contents = infile.readlines()
        transcripts_ID_list = [line.strip() for line in contents]
    
    return transcripts_ID_list

# Read text files containing list of 
# transcript names + length of UTR5 or UTR3 (bioawk output file)
# Example line: transcript:KAF2020338	314
# Save into dictionary {transcriptID: length}
def get_transcripts_length_dict(length_file):
    length_dict = {}
    with open(length_file, "r") as infile:
        for line in infile:
            ID_length = line.strip().split("\t")            
            if len(ID_length) == 2:
                ID, length = ID_length 
                ID = re.sub(re.escape("transcript:"), "", ID)
                    # remove the "transcript:" in front of ID
                length_dict[ID] = length    
    return(length_dict)

# Read text files containing list of 
# transcript names + GC content (decimal) of UTR5 or UTR3 (bioawk output file)
# Example line: transcript:KAF2020362   0,434629
# Save into dictionary {transcriptID: GC content}
def get_transcripts_GC_dict(GC_file):
    GC_dict = {}
    with open(GC_file, "r") as infile:
        for line in infile:
            ID_GC = line.strip().split("\t")
            if len(ID_GC) == 2:
                ID, GC = ID_GC
                ID = re.sub(re.escape("transcript:"), "", ID)
                    # remove the "transcript:" in front of ID
                GC = GC.replace(",", ".")
                    # replace comma with decimal point
                GC_dict[ID] = GC
    return(GC_dict)

# Create 2150 bp sequence from separate prom/UTR5/UTR3/term fasta files
def create_2150_NT_seq(prom_file, UTR5_file, UTR3_file, term_file, ID_list):
    prom_full = SeqIO.index(prom_file, "fasta")
    UTR5_full = SeqIO.index(UTR5_file, "fasta")
    UTR3_full = SeqIO.index(UTR3_file, "fasta")
    term_full = SeqIO.index(term_file, "fasta")

    regulatory_seq_list = [] # store 2150bp NT sequence

    for entry in ID_list:
        regulatory_seq_working = str()
        temp_identifier = "transcript:" + entry
        regulatory_seq_working += str((prom_full[temp_identifier]).seq)
        regulatory_seq_working += str((UTR5_full[temp_identifier]).seq)
        regulatory_seq_working += str((UTR3_full[temp_identifier]).seq)
        regulatory_seq_working += str((term_full[temp_identifier]).seq)

        regulatory_seq_list.append(regulatory_seq_working.upper())
    
    return regulatory_seq_list

# Function shared and authored by Zrimec et al.
# One hot encodes nucleotide sequence, returns numpy array
# Needs uppercase input
def NT_to_onehot(data):
    alphabet = 'ACGT'
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
    integer_encoded = []
    for char in data:
        if char in alphabet:
            integer_encoded.append(char_to_int[char])
        else:
            integer_encoded.append(0)
    onehot = np.array(to_categorical(integer_encoded,num_classes=4),dtype=np.int8)
    onehot[[char not in alphabet for char in data]] = [0,0,0,0]
    return onehot

# Create array of 2150 bp one hot encoded sequences
def get_x_hot_arr(regulatory_nt_list):
    # Create list of 2150 bp one hot encoded sequences
    regulatory_hot_list = []
    for sequence in regulatory_nt_list:
        regulatory_hot_list.append(NT_to_onehot(sequence.upper()))
    # Save list of 2150 bp one hot encoded sequences to numpy array
    x_hot = np.asarray(regulatory_hot_list, dtype=np.int8)
    return x_hot

# Create array containing the 72 additional variables:
## 8 mRNA stability variables: 
###     lengths of UTR5, CDS, UTR3 (total = 3)
###     GC of UTR5, UTR3, CDS codon positions 1, 2, 3  (total = 5)
## 64 codon usage counts
def get_x_var_arr(ID_list, length_UTR5_dict, length_EGFP, length_UTR3_dict, 
                  GC_UTR5_dict, GC_UTR3_dict, GC_EGFP_list, EGFP_codon_list):
    add_var_list = []
    for entry in ID_list:
        add_var_working = []
        
        # Append length values
        add_var_working.append(length_UTR5_dict[entry])
        add_var_working.append(length_EGFP)
        add_var_working.append(length_UTR3_dict[entry])

        # Append GC values
        # Multiply GC content by 1000 and round to nearest integer
        add_var_working.append(round(float(GC_UTR5_dict[entry])*1000))
        add_var_working.append(round(float(GC_UTR3_dict[entry])*1000))
        for decimal in GC_EGFP_list:
            add_var_working.append(round(float(decimal)*1000))

        # Extend with list of codon usage counts
        add_var_working.extend(EGFP_codon_list)

        # Append working list for this entry to main list
        add_var_list.append(add_var_working)

    x_var = np.asarray(add_var_list, dtype=np.int16)
    return x_var


# Input filepaths
# CHANGE THESE AS NECESSARY
working_directory = "./project_data/inputs"

transcripts_ID_directory = working_directory + "/temp/transcripts_completeset_index"

transcripts_length_UTR5_directory = working_directory + "/lengths/UTR5_len"
transcripts_length_UTR3_directory = working_directory + "/lengths/UTR3_len"

transcripts_GC_UTR5_directory = working_directory + "/GC/UTR5_GC"
transcripts_GC_UTR3_directory = working_directory + "/GC/UTR3_GC"

EGFP_GC_file = working_directory + "/EGFP/EGFP_codon_GC.npy"
EGFP_codon_file = working_directory + "/EGFP/EGFP_codon_usage.npy"

prom_directory = working_directory + "/padded_cut/PROMOTER_final_V2"
UTR5_directory = working_directory + "/padded_cut/UTR5_final"
UTR3_directory = working_directory + "/padded_cut/UTR3_final"
term_directory = working_directory + "/padded_cut/TERMINATOR_final_V2"

# Output directory
npz_output_directory = working_directory + "/npz_files"

# Make list containing all species names
# from txt file containing each species on new line
species_list = get_transcripts_ID_list(working_directory + "/species_list.txt")


# !!BEGIN ASSEMBLING DATASETS!!

for species in species_list:
    print("Now processing species: " + species)
    
    # Define input and output filepaths
    transcripts_ID_file = glob.glob(str(transcripts_ID_directory + "/" + species + "*"))
    
    transcripts_length_UTR5_file = glob.glob(str(transcripts_length_UTR5_directory + "/" + species + "*"))
    transcripts_length_UTR3_file = glob.glob(str(transcripts_length_UTR3_directory + "/" + species + "*"))
    
    transcripts_GC_UTR5_file = glob.glob(str(transcripts_GC_UTR5_directory + "/" + species + "*"))
    transcripts_GC_UTR3_file = glob.glob(str(transcripts_GC_UTR3_directory + "/" + species + "*"))
    
    prom_file = glob.glob(str(prom_directory + "/" + species + "*"))
    UTR5_file = glob.glob(str(UTR5_directory + "/" + species + "*"))
    UTR3_file = glob.glob(str(UTR3_directory + "/" + species + "*"))
    term_file = glob.glob(str(term_directory + "/" + species + "*"))

    npz_out_file = npz_output_directory + "/" + species + ".npz"


    # MAKE x_hot array
    print("Creating one hot encoded sequences...")
    # Load in list of transcripts
    ID_list = get_transcripts_ID_list(transcripts_ID_file[0]) 
        # glob output is technically a list so use [0]
    # Create list of 2150 bp NT sequences
    regulatory_nt_list = create_2150_NT_seq(prom_file[0], UTR5_file[0], UTR3_file[0], term_file[0], ID_list)
    # Now make the array
    x_hot = get_x_hot_arr(regulatory_nt_list) # np.save(x_hot_output_file, x_hot)
    print("Done")


    # MAKE x_var array; 72 aforementioned variables
    print("Processing additional 72 variables...")
    # Get lengths
    length_UTR5_dict = get_transcripts_length_dict(transcripts_length_UTR5_file[0])
    length_EGFP = int(720)
    length_UTR3_dict = get_transcripts_length_dict(transcripts_length_UTR3_file[0])
    # Get GC
    GC_UTR5_dict = get_transcripts_GC_dict(transcripts_GC_UTR5_file[0])
    GC_EGFP_list = (np.load(EGFP_GC_file)).tolist()
    GC_UTR3_dict = get_transcripts_GC_dict(transcripts_GC_UTR3_file[0])
    # Get EGFP codon usage counts
    EGFP_codon_list = (np.load(EGFP_codon_file)).tolist()
    # Now make the array
    x_var = get_x_var_arr(ID_list, length_UTR5_dict, length_EGFP, length_UTR3_dict, 
                    GC_UTR5_dict, GC_UTR3_dict, GC_EGFP_list, EGFP_codon_list)
    print("Done")


    # MAKE x_name array
    print("Saving names to array...")
    x_name = np.asarray(ID_list, dtype=object)
    print("Done")


    # MAKE output .npz file
    print("Creating output .npz file...")
    np.savez(npz_out_file, x_hot=x_hot, x_var=x_var, x_name=x_name)
    print("Done")

    print("Completed assembling dataset for species: " + species)

print("All done!")