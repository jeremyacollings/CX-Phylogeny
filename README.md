# CX-Phylogeny

Analysis of phylogenetic and trait influence on coexistence. 

Currently, the code is broken up into the following scripts: 

- __1_Data.Prep.R:__ Brings in raw trait data, fitness data, and phylogenetic distances and prepares dataframes for analysis. 
- __2_Fit.Pop.Models.R:__ Brings in fitness data and fits Ricker population model to estimate intrinsic growth rates and competition coefficients. 
- __2.5_Plot.Pop.Params.R:__ Plots estimates of intrinsic growth rates and competition coefficients.
- __3_Fit.Path.Analysis.R:__ Brings in trait and phylogenetic data as well as Ricker model estimates, calculates coexistence metrics (niche differences, demographic ratios, competitive ratios, fitness differences, and coexistence probabilities) and trait distances, and fits path analyses.