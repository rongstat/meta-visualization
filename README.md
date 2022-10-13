# Spectral Method for Assessing and Combining Multiple Data Visualizations


We present an effcient spectral method for assessing and combining multiple visualizations of a given dataset produced by diverse algorithms. The proposed method provides a quantitative measure -- the visualization eigenscore -- of the relative performance of the visualizations for preserving the structure around each data point. Then it leverages the eigenscores to obtain a consensus visualization, which has much improved quality over the individual visualizations in capturing the underlying true data structure. Our approach is flexible and works as a wrapper around any visualizations.

The method is based on the paper:

Ma, R., Sun, E., and Zou, J. (2022) A Spectral Method for Assessing and Combining Multiple Data Visualizations. 

# Content

The R script `main_fun.R` contains the main functions for producing multiple data visualizations, their eigenscores, and the meta-visualization.

The directory `Data` contains several example datasets for demonstrating our methods.

The directory `R Codes` includes R scripts for analyzing the example datasets, and for reproducing the data analysis in the manuscript.

The directory `Python Code` includes Python implementation of the method.

# System Requirements

The meta-visualization package requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method are tested under R version 4.1.1, and require the R packages: `rARPACK`,`MASS`,`lle`,`dimRed`,`uwot`,`cluster`,`Rtsne`,`phateR`.


# Get Started

The main functions for meta-visualization are contained in `main_fun.R`. The function `candidate.visual()` helps to produce diverse candidate visualizations based on our choice of dimension reduction methods. The users can also generate candidate visualizations on their own. The function `ensemble.viz()` takes the candidates visualizations as inputs and returns their eigenscores and a final meta-visualization.

To apply our method to an example dataset, follow the three steps below.

1. Download the dataset `AllBooks_baseline_DTM_Labelled.csv` from the directory `Data`. 
2. Load the R script `main_function.R` in the directory `R Codes`.
3. Run the R script `ReligionText.R` in the directory `R Codes`.

Note that the main difference between the function `ensemble.viz()` in `main_fun.R` and the function `ensemble.v.local` in `main_function.R` is whether the pairwise distance matrices for the candidate visualizations are returned. For general applications, we recommend `ensemble.viz()` in `main_fun.R` as it requires less memory, especially for large datasets.

For further questions and inquiries, please contact Rong Ma (rongm@stanford.edu).
