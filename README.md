# Serosurveys to Optimize Peak Pathogen Predictions (STOPPP) Model

## GENERAL INFORMATION
Repository for the Mathematica notebook, R code and previoulsy published data used in analysesfor the manuscript:

Erin Clancey $^{1,\ast}$, Scott L. Nuismer $^2$, and Stephanie N. Seifert $^1$, Using serosurveys to optimize surveillance for zoonotic pathogens

1. Current Address: Paul G. Allen School for Global Health, Washington State University, Pullman, WA 99164 USA;
2. Department of Biological Sciences, University of Idaho, Moscow, ID 83844 USA:

$\ast$ Corresponding author; e-mail: erin.clancey@wsu.edu

Zoonotic pathogens pose significant risk to human health, with spillover into human populations contributing to chronic disease and epidemics. Despite the widely recognized burden of zoonotic spillover, our ability to identify which animal populations serve as primary reservoirs remains incomplete. This challenge is compounded when prevalence in reservoir populations reaches detectable levels only at specific times of year. In these cases, statistical models designed to predict the timing of peak prevalence could guide field sampling for active infections or predict when spillover risk is likely to be greatest. We develop a general mathematical model that leverages routinely collected serosurveillance data to optimize sampling for elusive pathogens. Using simulated data, we show that our methodology reliably identifies times when pathogen prevalence is expected to peak. Then, we demonstrate an implementation of our method using previously published serosurvey data in straw-colored fruit bats (\textit{Eidolon helvum}). The generality and simplicity of our methodology make it broadly applicable to a wide range of putative reservoir species where seasonal patterns of birth lead to cyclic, but potentially short-lived, pulses of pathogen prevalence.

All data used in this study was previously published and is available from Pleydell, D. R. J. (2023). R, C and NIMBLE code for Yaounde-Ebola-E.helvum modelling paper (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8193102

Funding was provided by the Centre for Research in Emerging Infectious Diseases - East and Central Africa (CREID-ECA) grant number U01AI151799 in support of EC and SNS, PIPP Phase I: Predicting Emergence in Multidisciplinary Pandemic Tipping-points (PREEMPT) from the U.S. National Science Foundation (NSF) grant number 2200140 in support of EC, Verena (viralemergence.org) from NSF including NSF BII 2021909 and NSF BII 2213854 and the US National Institute of Allergy and Infectious Disease/National Institutes of Health (NIAID/NIH) grant number U01AI151799 in support of SNS, NIH 2R01GM122079-05A1 awarded to SNL and NSF DEB 2314616 awarded to SNL. Funding was also provided by the cooperative agreement CDC-RFA-FT-23-0069 from the CDC’s Center for Forecasting and Outbreak Analytics in support of EC and SNS. Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the Centers for Disease Control and Prevention.

## DATA & CODE FILE OVERVIEW
The repository is split into three subdirectoires: 
1. **Mathematica** - This directory contains the Mathematica Notebook.
2. **Simulations** - This directory contains R code used to create and analyze the simulated data. 
3. **Real Data** - This directory contains E. helvum data and R code used to analyze the data. 
