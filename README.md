# meta-visualization
A spectral method for assessing and combining multiple data visualizations.

For a collection of visualizations of a dataset, the method provides a quantitative measure -- eigenscore -- of relative performance of each visualization in a pointwise manner, which induces a ranking of these visualizations in terms of their concordance to the underlying true structures of the data. On the other hand, the method automatically combines strengths and ameliorates weakness of multiple visualizations, leading to a meta-visualization which is  provably better than all the candidate visualizations. 

# Use Guide

The R script `main_function.R` contains the main functions for producing multiple data visualizations and the meta-visualization.

To apply our method to an example dataset of religious and biblical text, download the dataset `AllBooks_baseline_DTM_Labelled.csv` and simply run the R codes in `ReligionText.R`.



