# disinfection_sensitivity_dist

# What is this page?
This page contains the R script to estimate disinfectant sensitivity distributions using a Bayesian framework.

# Objective
Disinfection efficacy against viruses has been extensively studied. By integrating reported experimental data on inactivation rate constants for various viruses, we propose a statistical framework to analyze variation in viral sensitivity to disinfectants through a Bayesian approach. This R code constructs parametric disinfectant sensitivity distributions for three disinfectants: ultraviolet, ozone, and free chlorine.

# Files
virus_ssd_github_v3.R
An example R code for analysis and visualization. This code requires an input data set.

Summary_k_UV.xlsx, Summary_k_Ozone.xlsx, Summary_k_FreeChlorine.xlsx
These files include aggregated values of inactivation rate constants for ultraviolet (UV), ozone, and free chlorine, respectively. The original data were collected from the papers below:
Rockey, N. C.; Henderson, J. B.; Chin, K.; Raskin, L.; Wigginton, K. R. Predictive Modeling of Virus Inactivation by UV. Environ. Sci. Technol. 2021, 55 (5), 3322–3332.
Morrison, C. M.; Hogard, S.; Pearce, R.; Gerrity, D.; von Gunten, U.; Wert, E. C. Ozone Disinfection of Waterborne Pathogens and Their Surrogates: A Critical Review. Water Res. 2022, 214, 118206.
Chaplin, M.; Henderson, J.; Leung, K.; Szczuka, A.; Hansen, B.; Rockey, N. C.; Wigginton, K. Linear Mixed Model of Virus Disinfection by Chlorine to Harmonize Data Collected Across Broad Environmental Conditions. bioRxiv, 2023, 2023.09.12.557160. https://doi.org/10.1101/2023.09.12.557160.


# Publication
In preparation
