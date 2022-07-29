# meta-visualization
A spectral method for assessing and combining multiple data visualizations.

For a collection of visualizations of a dataset, the method provides a quantitative measure -- eigenscore -- of relative performance of each visualization in a pointwise manner, which induces a ranking of these visualizations in terms of their concordance to the underlying true structures of the data. On the other hand, the method automatically combines strengths and ameliorates weakness of multiple visualizations, leading to a meta-visualization which is  provably better than all the candidate visualizations. 

# Content

The R script `main_function.R` contains the main functions for producing multiple data visualizations, their eigenscores, and the meta-visualization.

The directory `Data` contains several example datasets for demonstrating our methods.

The directory `R Codes` includes R scripts for analyzing the example datasets.

The directory `Python Code` includes Python implementation of the method.

# System Requirements

The meta-visualization package requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method are tested under R version 4.1.1, and require the R packages: `rARPACK`,`MASS`,`lle`,`dimRed`,`uwot`,`cluster`,`Rtsne`,`phateR`.


# Get Started

The main functions for meta-visualization are contained in `main_function.R`. The function `candidate.visual()` helps to produce diverse candidate visualizations based on our choice of dimension reduction methods. The users can also generate candidate visualizations on their own. The function `ensemble.v.local()` takes the candidates visualizations as inputs and returns their eigenscores and a final meta-visualization.

To apply our method to an example dataset, follow the three steps below.

1. Download the dataset `AllBooks_baseline_DTM_Labelled.csv` from the directory `Data`. 
2. Load the R script `main_function.R`.
3. Run the R script `ReligionText.R` in the directory `R Codes`.



For further questions and inquiries, please contact Rong Ma (rongm@stanford.edu).
