# Workflow for the STOPPP model with real data

## GENERAL INFORMATION
Workflow to demonstrate an implementation of the STOPPP model with previously published serosurveillance data. Each subdirectory contains the same files specific to each species. Data and R code in each file will re-create the analyses and figures presented in Clancey et al. "Using serosurveys to optimize surveillance for zoonotic pathogens." 

## DATA & CODE FILE OVERVIEW
The repository is split into two subdirectoires: 
1. **E_helvum_cam** - This directory contains information to follow a workflow for *E. helvum* with the following files: <br>
   a. E_helvum_data.csv <br>
   b. Data_Packages_E_helvum.R <br>
   c. Functions_E_helvum.R <br>
   d. Fit_data_STAN_lac_Ehelvum.R <br>
   e. Fit_data_STAN_EBOV_E_helvum.R <br>
   f. Plot_E_helvum_Fit.R
2. **H_monstrosus_cam** - This directory contains information to follow a workflow for *H. monstrosus* with the following files: <br>
   a. H_monstrosus_data.csv <br>
   b. Data_Packages_H.monstrosus.R <br>
   c. Functions_H.monstrosus.R <br>
   d. Fit_data_STAN_Lac_H.monstrosus.R <br>
   e. Fit_data_STAN_EBOV_H.monstrosus.R <br>
   f. Plot_H_monstrosus_Fit.R

## Instructions to run the R code
Download files a-f for each species and save them in separate folders. Files d and e will source the information in files a, b, and c and fit the data using rSTAN. File f will print the results from interpolation and model fitting and plot the results. 

 
