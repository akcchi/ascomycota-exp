# For the 70,225 prom/term combinations expression levels
# Conduct 2-tailed Wilcoxon rank-sum test (combinations compared to yeast)
# Display as violin plot
# Also display the top 5 expressing genes in yeast as vertical lines

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
from scipy.stats import ranksums

fpath_predictions = "./predictions/"
fpath_predictions_yeast = "./predictions_YEAST/Saccharomyces_cerevisiae_R64-1-1.Sacce1_MODIFIED.npy"

fpath_predictions_combo = "./predictions_COMBO/swapped_PREDICTIONS.npy"

output_directory = "./figures/"

# Load combo predictions
combo_predictions_load = np.load(fpath_predictions_combo, allow_pickle=True)
data_temp = [array[0] for array in combo_predictions_load]
combo_predictions = pd.Series(data_temp)
del combo_predictions_load

# Load in yeast predictions
yeast_load = np.load(fpath_predictions_yeast, allow_pickle=True)
yeast_temp = [array[0] for array in yeast_load]
yeast_predictions = pd.Series(yeast_temp)
del yeast_load

# Save top 5 expressing yeast genes
yeast_top = (yeast_predictions.nlargest(n=5)).tolist()

# Conduct Wilcoxon ranksum test
w, p_val = ranksums(data_temp, yeast_temp, alternative="two-sided")
del data_temp
del yeast_temp

plot_series_list = [yeast_predictions, combo_predictions]

dframe = pd.concat(plot_series_list, axis=1)
dframe.columns = ["Sacce1 (n=425)", "Combinations (n=70225)"]
# print(dframe)

dframe_melted = dframe.melt()
dframe_melted.columns = ["source", "prediction"]
dframe_melted.dropna(axis=0, inplace=True) # remove rows with NaN
# print(dframe_melted)

print("predictions loaded")

# Format figure style
params = {"figure.dpi": 300, 
          "axes.edgecolor": "black", 
          "axes.spines.top": False, "axes.spines.right": False,
          "xtick.bottom": True, "ytick.left": True}

sns.set_theme(style="whitegrid", 
              rc=params)

# Violinplot
g = sns.violinplot(data=dframe_melted, x="prediction",
                   y="source", inner=None,
                   color="gainsboro", linecolor="black", linewidth=0.9)

# Boxplot inside violinplot
sns.boxplot(data=dframe_melted, x="prediction", y="source",
             color="gainsboro", linecolor="black",
               boxprops={"zorder": 2}, ax=g, 
             width=0.15, showfliers=False, showcaps=False)

# Format plot labels and (if any) with statistical significance
if p_val <= 0.001:
    g.set(yticklabels=["Sacce1\n(n=" + str(len(yeast_predictions)) + ")", 
                   "*** Combinations\n(n=" + str(len(combo_predictions)) + ")"])
elif p_val <= 0.01:
    g.set(yticklabels=["Sacce1\n(n=" + str(len(yeast_predictions)) + ")", 
                   "** Combinations\n(n=" + str(len(combo_predictions)) + ")"])
elif p_val <= 0.05:
    g.set(yticklabels=["Sacce1\n(n=" + str(len(yeast_predictions)) + ")", 
                   "* Combinations\n(n=" + str(len(combo_predictions)) + ")"])
else:
    g.set(yticklabels=["Sacce1\n(n=" + str(len(yeast_predictions)) + ")", 
                   "Combinations\n(n=" + str(len(combo_predictions)) + ")"])


# Add max expressing yeast levels
for pred in yeast_top:
    plt.axvline(x=pred, color="red", linestyle="--", linewidth=0.8)

# Axes labels
plt.xlabel("Predicted gene expression level (TPM)")
plt.ylabel("Gene origin")

plt.savefig(output_directory + "combinations_plot_v5.1.svg", format = "svg", bbox_inches="tight")
print("graph done")