Code repository for research project (January to March 2024) conducted as part of the BSc Biomedical Science programme at KCL.

## Predicting and Optimising Gene Expression in _Saccharomyces cerevisiae_
This repository was intended as supplementary information for the project write-up.

[preprocessing](/preprocessing): contains scripts used to extract and preprocess nucleotide sequences (from the [Ensembl Fungi](https://fungi.ensembl.org/) database) into input datasets compatible with the _Saccharomyces cerevisiae_ gene expression neural network model trained by Zrimec et al., as described in their 2020 paper (DOI: [10.1038/s41467-020-19921-4](https://doi.org/10.1038/s41467-020-19921-4)).

[predictions](/predictions): code to generate predictions with the aforementioned pre-trained neural network model, using the prepared input datasets.

[analysis](/analysis): code to reproduce the data analysis and visualisations as reported in the write-up. Some of the figures can be browsed at [figs](/figs); see below for one of them.

---
<p align="center">
  <img src=https://github.com/akcchi/ascomycota-exp/blob/main/figs/ridgeline_all.png alt="ridgeline plot comparing predicted gene expression levels" width="600">
</p>
<p>
  Figure: Predicted gene expression levels of foreign Ascomycota genes in <i>Saccharomyces cerevisiae</i> generated using the pre-trained neural network model. Ridgeline plot showing prediction distributions of <i>n</i> number of genes grouped by genomic origin and ordered by median. Predicted gene expression levels as mRNA abundance in transcripts per million (TPM). Red line indicates median predicted expression level for native <i>S. cerevisiae</i> (Sacce1) gene sample. Two-tailed Wilcoxon rank-sum test was used to compare each prediction distribution to that of Sacce1. Asterisks (*) denote statistical significance as follows: * p < 0.05, ** p < 0.01, *** p < 0.001.
</p>
