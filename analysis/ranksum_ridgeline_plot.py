# For the >300,000 predicted expression levels for 52 fungi species + yeast
# Conduct 2-tailed Wilcoxon rank-sum test (each species compared to yeast)
# And display in ridgeline plot (ordered by median)
# Also add vertical line indicating median expression level for yeast dataset

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
from scipy.stats import ranksums


fpath_predictions = "./predictions/"
fpath_predictions_yeast = "./predictions_YEAST/Saccharomyces_cerevisiae_R64-1-1.Sacce1_MODIFIED.npy"
output_directory = "./figures/"

fpath_list = glob.glob(fpath_predictions + "*")
fpath_list.sort()
fpath_list.append(fpath_predictions_yeast)

# Make list of species names
species_names = list()
for file in fpath_list:
    working_name = file.replace("./predictions/", "")
    working_name = working_name.replace("./predictions_YEAST/", "")
    working_name = working_name.replace("_MODIFIED", "")
    # Split names by periods, get abbreviated name
    working_name = working_name.split(".")
    species_names.append(working_name[1])
del working_name

series_list = []
# for i, file in enumerate(fpath_list):
for file in fpath_list:
    data_load = np.load(file, allow_pickle=True)
    data_temp = [array[0] for array in data_load]
    series = pd.Series(data_temp)
    series_list.append(series)
del data_load
del data_temp
del series

# Create a dataframe containing 2 columns:
## fpath_list index
## median expression lvl of that species
median_list = []
counter = []
for i, species in enumerate(series_list):
    median_list.append(species.median())
    counter.append(i)
median_df = pd.DataFrame({"fpath_list_counter": counter, "median": median_list})

# Sort dataframe by median value, highest at top (descending)
median_df_sorted = median_df.sort_values(by=["median"], ascending=False)

# Get the ordered fpath_list_counter, use this to plot
fpath_list_counter = median_df_sorted["fpath_list_counter"].tolist()

plot_series_list = []
for i in fpath_list_counter:
    # print(fpath_list[i])
    data_load = np.load(fpath_list[i], allow_pickle=True)
    data_temp = [array[0] for array in data_load]
    series = pd.Series(data_temp)
    plot_series_list.append(series)
del data_load
del data_temp
del series

# Load in yeast predictions for ranksum tests
yeast_load = np.load(fpath_list[52], allow_pickle=True)
yeast_predictions = [array[0] for array in yeast_load]
del yeast_load

# plot_species_names = []
# for i in fpath_list_counter:
#     plot_species_names.append(species_names[i])

plot_species_names = []
# plot_species_sizes = []
# Do ranksums, append stat. sig. to species name
for i in fpath_list_counter:
    # This looks insane but needed for sample size labels on right hand side
    buffer = " " * 2
    buffer_temp = buffer + (" " * (11 - len(species_names[i])))

    if i == 52: # yeast test dataset
        append_temp = "    " + species_names[i] + buffer_temp + "n=" + str(len(yeast_predictions))
        plot_species_names.append(append_temp)
        break
    
    data_load = np.load(fpath_list[i], allow_pickle=True)
    data_temp = [array[0] for array in data_load]
    del data_load
    
    w, p = ranksums(data_temp, yeast_predictions, alternative="two-sided")
    # print(p)
    if p <= 0.001:
        append_temp = "*** "+ species_names[i] + buffer_temp + "n=" + str(len(data_temp))
    elif p <= 0.01:
        append_temp = " ** "+ species_names[i] + buffer_temp + "n=" + str(len(data_temp))
    elif p <= 0.05:
        append_temp = "  * "+ species_names[i] + buffer_temp + "n=" + str(len(data_temp))
    else:
        append_temp = "    "+ species_names[i] + buffer_temp + "n=" + str(len(data_temp))

    plot_species_names.append(append_temp)


dframe = pd.concat(plot_series_list, axis=1)

dframe.columns = plot_species_names

# print(dframe)

dframe_melted = dframe.melt()
    # get:
    # # variable, value
    # # Aaor1, decimal
    # # Aaor1, decmial
    # # Sacce1, decimal
    # etc.

# name columns to species and prediction
dframe_melted.columns = ["species", "prediction"]

# print(dframe_melted)

dframe_melted.dropna(axis=0, inplace=True) # remove rows with NaN

dframe_melted["prediction"] = pd.to_numeric(dframe_melted["prediction"], downcast="float")
    # convert predictions to float

print("dframe done")

sns.set_theme(style="white", 
              rc={"figure.dpi": 300, 
                  "axes.facecolor": (0, 0, 0, 0)})

# g = sns.FacetGrid(dframe_melted, row="species")
g = sns.FacetGrid(dframe_melted, row="species", aspect=40, height=0.25) # make plots wider
    # Each subplot on new row, each row is 
    # the same species as specified in the "species" column

g.map_dataframe(sns.kdeplot, x="prediction", color="gainsboro", fill=True, alpha=0.5)
g.map_dataframe(sns.kdeplot, x="prediction", color="black", zorder=100) # Trace black outline

g.figure.subplots_adjust(hspace=-0.8) # Reduce height space (overlap)

g.set_titles("") # Remove titles

g.set(yticks=[], # Hide y axis ticks
      xlabel="Predicted gene expression level (TPM)") # Label bottom x-axis

g.despine(left=True) # Hide y-axis line

for i, ax in enumerate(g.axes.flat):
    # Label each distribution with species name
    if i == 52: # Red for yeast
        ax.set_ylabel(plot_species_names[i], rotation=0, loc="bottom", color="red", fontname="Liberation Mono")
    else:
        ax.set_ylabel(plot_species_names[i], rotation=0, loc="bottom", fontname="Liberation Mono")

    ax.yaxis.set_label_coords(-0.3, 0.0) # Move y-labels to the left

    ax.spines["bottom"].set_color("lightgrey")

    ax.axvline(x=median_list[52], color="red", linestyle="solid")


# g.axvline(x=median_list[52], ymax=1, color="red", linestyle="dotted")

plt.savefig(output_directory + "ordered_ridge_9.svg", format = "svg", bbox_inches="tight")
print("Done")