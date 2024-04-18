Code to generate predicted expression levels for: 
+ previously assembled dataset from 52 fungi genomes (see [here](../preprocessing/assemble_npz.py))
+ _Saccharomyces cerevisiae_ test dataset, extracted from the [Zenodo repository](https://doi.org/10.5281/zenodo.3905252) by Jan Zrimec
+ 70,225 promoter/terminator combinations (see hereTBA)

Uses the _Saccharomyces cerevisiae_ expression neural network model published by [Zrimec et al.](https://doi.org/10.1038/s41467-020-19921-4);
move the contents of the [Zenodo repository](https://doi.org/10.5281/zenodo.3905252) by J. Zrimec into a directory named "data".

# Dependencies
See conda [environment.yml](/predictions/environment.yml) file.

# Acknowledgements
Code copied and modified (MIT License) from [https://github.com/JanZrimec/DeepExpression/blob/master/scripts/Chapter_3_1.ipynb](https://github.com/JanZrimec/DeepExpression/blob/master/scripts/Chapter_3_1.ipynb)<br>
Original author: Jan Zrimec ([https://github.com/JanZrimec](https://github.com/JanZrimec))

Code modified to take the different input datasets (as described above).

