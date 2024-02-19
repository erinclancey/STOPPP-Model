# STOPPP-Model

## GENERAL INFORMATION
Repository for the Mathematica notebook, R code and previoulsy publisheed data used in analysesfor the manuscript:

Erin Clancey $^{1,\ast}$, Scott L. Nuismer $^2$, and Stephanie N. Seifert $^1$, Using serosurveys to optimize surveillance for zoonotic pathogens

1. Current Address: Paul G. Allen School for Global Health, Washington State University, Pullman, WA 99164 USA;
2. 2Department of Biological Sciences, University of Idaho, Moscow, ID 83844 USA:

$\ast$ Corresponding author; e-mail: erin.clancey@wsu.edu

Zoonotic pathogens pose a significant risk to human health, with spillover into human populations contributing to chronic disease, sporadic epidemics, and occasional pandemics. Despite the widely recognized burden of zoonotic spillover, our ability to identify which animal populations serve as primary reservoirs for these pathogens remains incomplete. This challenge is compounded when prevalence reaches detectable levels only at specific times of year. In these cases, statistical models designed to predict the timing of peak prevalence could guide field sampling for active infections.  Here we develop a general model that leverages routinely collected serosurveillance data to optimize sampling for elusive pathogens. Using simulated data sets we show that our methodology reliably identifies times when pathogen prevalence is expected to peak. We then apply our method to two putative \textit{Ebolavirus} reservoirs, straw-colored fruit bats (\textit{Eidolon helvum}) and hammer-headed bats (\textit{Hypsignathus monstrosus}) to predict when these species should be sampled to maximize the probability of detecting active infections. In addition to guiding future sampling of these species, our method yields predictions for the times of year that are most likely to produce future spillover events. The generality and simplicity of our methodology make it broadly applicable to a wide range of putative reservoir species where seasonal patterns of birth lead to predictable, but potentially short-lived, pulses of pathogen prevalence.

All data used in this study was previously published and is available from Pleydell, D. R. J. (2023). R, C and NIMBLE code for Yaounde-Ebola-E.helvum modelling paper (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8193102

Funding was provided by the Centre for Research in Emerging Infectious Diseases—East and Central Africa (CREID-ECA) grant number U01AI151799 (SNS), Verena (viralemergence.org) from the U.S. National Science Foundation (NSF), including NSF BII 2021909 and NSF BII 2213854 and the US National Institute of Allergy and Infectious Disease/National Institutes of Health (NIAID/NIH), grant number U01AI151799 (SNS), NIH 2R01GM122079-05A1 (SLN) and NSF DEB 2314616 (SLN).

## DATA & CODE FILE OVERVIEW
The repository is split into four subdirectoires: 
1. **E_Helvum_Cam** - This directory contains the data and R code for the *E. helvum* population.
2. **H_Monstrosus_Cam** - This directory contains the data and R code for the *H. monstrosus* population.
3. **Mathematica** - This directory contains the Mathematica Notebook.
4. **Simulations** - This directory contains R code used to create and analyze the simulated data. 
