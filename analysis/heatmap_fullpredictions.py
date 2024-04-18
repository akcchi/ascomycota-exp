# For the >300,000 gene sequences
# Create heatmap comparing the means of the following gene regulatory 
# sequence variables to that of yeast using standard score:
## Expression level
## Length of 5' UTR
## Length of 3' UTR
## GC content of 5' UTR
## GC content of 3'UTR

import numpy as np
from statistics import mean
from statistics import stdev
import pandas as pd
from statsmodels.stats.weightstats import ztest as ztest
import seaborn as sns
import matplotlib.pyplot as plt
import glob

# Take in 72-variable feature array of one species
# Return one list containing four elements, means of:
#   5UTR length
#   3UTR length
#   5UTR GC
#   3UTR GC
def extract_means(var_arr):
    len_5UTR = []
    len_3UTR = []
    gc_5UTR = []
    gc_3UTR = []
    for arr in var_arr:
        len_5UTR.append(arr[0])
        len_3UTR.append(arr[2])
        gc_5UTR.append(arr[3])
        gc_3UTR.append(arr[4])
    means_species = [mean(len_5UTR), mean(len_3UTR), mean(gc_5UTR), mean(gc_3UTR)]
    return(means_species)

# Take in 72-variable feature array of one species
# Return 4 lists:
#   Lengths of 5UTR
#   Lengths of 3UTR
#   GC of 5UTR
#   GC of 3UTR
def extract_values(var_arr):
    len_5UTR = []
    len_3UTR = []
    gc_5UTR = []
    gc_3UTR = []
    for arr in var_arr:
        len_5UTR.append(arr[0])
        len_3UTR.append(arr[2])
        gc_5UTR.append(arr[3])
        gc_3UTR.append(arr[4])
    # means_species = [mean(len_5UTR), mean(len_3UTR), mean(gc_5UTR), mean(gc_3UTR)]
    return(len_5UTR, len_3UTR, gc_5UTR, gc_3UTR)

# Take in 72-variable feature array of one species
# Return one list containing four elements
# Standard scores of: means of four variables, compared to that of yeast
#   Lengths of 5UTR
#   Lengths of 3UTR
#   GC of 5UTR
#   GC of 3UTR 
# Standard score = (mean - mean(yeast)) / std(yeast)
def standardise_per_species(predictions_list, features_arr):
    len_5UTR = []
    len_3UTR = []
    gc_5UTR = []
    gc_3UTR = []

    # Extract the relevant feature from each gene
    # Append to its respective list
    for arr in features_arr:
        len_5UTR.append(arr[0])
        len_3UTR.append(arr[2])
        gc_5UTR.append(arr[3])
        gc_3UTR.append(arr[4])
    

    # Calculate standard score (normalised to yeast values)
    stdscore_exp = stdscore(mean(predictions_list), mean(yeast_exp), stdev(yeast_exp))
    stdscore_len_5UTR = stdscore(mean(len_5UTR), mean(yeast_len_5UTR), stdev(yeast_len_5UTR))
    stdscore_len_3UTR = stdscore(mean(len_3UTR), mean(yeast_len_3UTR), stdev(yeast_len_3UTR))
    stdscore_gc_5UTR = stdscore(mean(gc_5UTR), mean(yeast_gc_5UTR), stdev(yeast_gc_5UTR))
    stdscore_gc_3UTR = stdscore(mean(gc_3UTR), mean(yeast_gc_3UTR), stdev(yeast_gc_3UTR))

    stdscore_register = [stdscore_exp, stdscore_len_5UTR, stdscore_len_3UTR, 
                         stdscore_gc_5UTR, stdscore_gc_3UTR]

    return(stdscore_register)

# Function for standard score
# Standard score = (mean(working species) - mean(yeast)) / std(yeast)
def stdscore(mean, yeast_mean, yeast_stdev):
    standard_score = (mean - yeast_mean) / yeast_stdev
    return(standard_score)

# Filepaths
test_npz = "./inputs_testing/Aaosphaeria_arxii_cbs_175_79_gca_010015735.Aaoar1.npz"
yeast_npz = "./inputs_testing/Saccharomyces_cerevisiae_R64-1-1.Sacce1_MODIFIED.npz"

fpath_predictions = "./predictions/"
fpath_features = "./npz_files/"
fpath_predictions_yeast = "./predictions_YEAST/Saccharomyces_cerevisiae_R64-1-1.Sacce1_MODIFIED.npy"

output_directory = "./figures/"


# Get list of filepaths to PREDICTIONS
fpath_pred_list = glob.glob(fpath_predictions + "*")
fpath_pred_list.sort()

# Get list of filepaths to FEATURE ARRAYS (npz files)
fpath_features_list = glob.glob(fpath_features + "*")
fpath_features_list.sort()

# Make list of species names, same order as filepaths
species_names = list()
for file in fpath_pred_list:
    working_name = file.replace("./predictions/", "")
    working_name = working_name.replace("./predictions_YEAST/", "")
    working_name = working_name.replace("_MODIFIED", "")
    # Split names by periods, get abbreviated name
    working_name = working_name.split(".")
    species_names.append(working_name[1])
del working_name


# Load yeast variables
yeast = np.load(yeast_npz, allow_pickle=True)
yeast_var = yeast["og_var"]
yeast_len_5UTR, yeast_len_3UTR, yeast_gc_5UTR, yeast_gc_3UTR = extract_values(yeast_var)
del yeast
# Load yeast predicted exp levels
data_load = np.load(fpath_predictions_yeast, allow_pickle=True)
yeast_exp = [array[0] for array in data_load]
del data_load



# print(fpath_pred_list)


# For each species:
#   Load in predictions (one list)
#   Load in features (arrays of 72 features, x_var in .npz file)
#
#   Now do STANDARD SCORE compared to yeast, for
#       Predicted expression
#       Length of 5UTR
#       Length of 3UTR
#       GC of 5UTR
#       GC of 3UTR
stdscore_master_list = []
for i, file in enumerate(fpath_pred_list):
    # Get predictions
    pred_load = np.load(file, allow_pickle=True)
    pred_dataset = [array[0] for array in pred_load]
    del pred_load
    
    # Get features
    features_load = np.load(fpath_features_list[i], allow_pickle=True)
    features_dataset = features_load["x_var"]
    del features_load

    # Now get standard scores: this returns a list of 5 elements
    stdscore_register = standardise_per_species(pred_dataset, features_dataset)
    # Append to master list
    stdscore_master_list.append(stdscore_register)

print("standardising done...")




# Column names for dataframe
column_names = ["Predicted exp. level", "Length 5'UTR", 
                "Length 3'UTR", "GC content 5'UTR", "GC content 3'UTR"]

# Convert zscore master list to dataframe
dframe = pd.DataFrame(stdscore_master_list, index=species_names, columns=column_names)

# Sort dataframe by expression levels (highest at top/descending)
dframe_sorted = dframe.sort_values(by=["Predicted exp. level"], ascending=False)
print("dataframe done...")

# # Get index values (to plot names)
# index = list(dframe_sorted.index.values)

# print(dframe_sorted)


fig, ax = plt.subplots(figsize=(10, 14))
sns.set_theme(style="white", rc={"figure.dpi": 300})

sns.heatmap(dframe_sorted, vmin=-3, vmax=3, center=0, cmap="vlag", 
            yticklabels=True, linewidths=0.5)

plt.xlabel("Gene regulatory sequence variables compared to Sacce1")
ax.xaxis.labelpad = 10

plt.savefig(output_directory + "heatmap_2.1.svg", format = "svg", bbox_inches="tight")
print("Figure done!")
