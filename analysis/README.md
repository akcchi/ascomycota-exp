Code to reproduce the data analysis and visualisations, following gene expression level predictions (see ["predictions"](/predictions)).<br>
Uses the [SciPy](https://scipy.org/) and [seaborn](https://seaborn.pydata.org/) libraries.

# Dependencies
See conda [environment.yml](./environment.yml) file.


## ranksum_ridgeline_plot.py
Conduct two-tailed Wilcoxon rank-sum test for the predicted gene expression levels from each of 52 fungi species, compared to those of _Saccharomyces cerevisiae_; and generate ridgeline plot ordered by median.

## heatmap_fullpredictions.py
For the >300,000 gene sequences, create heatmap comparing the means of the following gene regulatory sequence variables to that of _S. cerevisiae_ using standard score:
+ Predicted expression level
+ Length of 5' UTR
+ Length of 3' UTR
+ GC content of 5' UTR
+ GC content of 3'UTR

## generate_combinations.ipynb
Get the predicted top 5 expressing genes from each species (including _S. cerevisiae_), and create "promoter/terminator" combinations, where 
+ "promoter" is simplified to promoter + 5'UTR regions.
+ "terminator" is simplified to 3'UTR + terminator regions.

## violinplot_combinations.py
Conduct two-tailed Wilcoxon rank-sum test for the predicted expression levels of the promoter/terminator combinations, compared to those of _Saccharomyces cerevisiae_; and visualise via violin plot.

## heatmap_combinations.py
For the promoter/terminator combinations, create heatmap comparing the means of the following gene regulatory sequence variables to that of _S. cerevisiae_ using standard score:
+ Predicted expression level
+ Length of 5' UTR
+ Length of 3' UTR
+ GC content of 5' UTR
+ GC content of 3'UTR
