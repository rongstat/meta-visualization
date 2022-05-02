# meta-visualization
A spectral method for assessing and combining multiple data visualizations.

For a collection of visualizations of a dataset, the method provides a quantitative measure -- eigenscore -- of relative performance of each visualization in a pointwise manner, which induces a ranking of these visualizations in terms of their concordance to the underlying true structures of the data. On the other hand, the method automatically combines strengths and ameliorates weakness of multiple visualizations, leading to a meta-visualization which is  provably better than all the candidate visualizations. 

# Content

The R script `main_function.R` contains the main functions for producing multiple data visualizations and the meta-visualization.

The directory `Data` contains several example datasets for demonstrating our methods.

The directory `R Codes` includes R scripts for analyzing the example datasets.

# Get Started

To apply our method to an example dataset, follow the three step below.

1. Download the dataset `AllBooks_baseline_DTM_Labelled.csv` from the directory `Data`. 
2. Run the R script `main_function.R`.
3. Run the R script `ReligionText.R` in the directory `R Codes`.



For further questions and inquiries, please contact Rong Ma (rongm@stanford.edu).
