.R files, capture history, and environmental covariate data for the Steelhead PDAT Project (under review)

##Last Updated 2025-09-17 by EMG

#This code accompanies the model used in the following publication:

#Quantifying component mortality estimates of out-migrating juvenile steelhead (Oncorhynchus mykiss) 
#Authors: E. M. Greenheck1, C. Michel2, B. Lehman2, L. Takata2, N. Demetras2, T. R. Nelson1
#1 George Mason University, Department of Environmental Science and Policy
#2 University of California Santa Cruz in affiliation with NOAA-NMFS-SWFSC
#[in progress]

#Kéry and Schaub 2012 (Bayesian population analysis using WinBUGS: a hierarchical perspective) was used to build these models

#GroupedModel_GitHub.R outlines how to estimate component mortality within groups (release group of fish)
#CovariateModel_GitHub.R outlines how to incorporate environmental and biological covariates when estimating component mortality

#Obs_matrix_all_new_20250711.csv is the capture history matrix for the project (need for GroupedModel_GitHub.R)
#The observation matrix also includes the grouping and individual-level covariates (e.g., Fork Length, Weight)
#turb_it.csv, discharge_it.csv, and temp_it.csv are environmental covariates (for each individual i at each time t) needed for CovariateModel_GitHub.R inaddition to the capture history matrix

