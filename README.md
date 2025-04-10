=================================================================================

SETUP INSTRUCTIONS FOR ZENODO AND GITHUB REPOSITORIES (https://doi.org/10.5281/zenodo.14020702):

The following provides an introduction and guide to the Zenodo and github repositories, for the pre-print manuscriupt presenting a parameterisation framework for Arctic sea ice emulation (SASIEv.1). The following scripts were run using the progammng interface MATLAB R2024b, with no additional add-ons or libraries.

This code is intended to showcase a parameterisation framework that compliments the in-text parameterisations and analysis provided, rather than a off the shelf, runnable tool. We intend for this setup to be a framework that can be reimplemented from the in-text parameterisations provided, and the functions and analysis provided in this repository. The scripts we provide here are therefore intended to be a reference that can be used to further understand the intext described process, it is not intended to be a ‘single click’ runnable tool. 

Data: The scripts in this repository is setup to run using only the datasets described in the adjoining paper. This version does not support flexible data types. It can only be run using the CMIP6 models outlined in the paper over the 1850-2100 time period, alongside the observational datasets referenced, the MAGICC global-mean surface temperature ensemble and the RCMIP CO2 emission datasets.
Observational sea ice products can be found at: https://www.cen.uni-hamburg.de/en/icdc/data/cryosphere/uhh-sea-ice-area-product.html
CMIP6 data can be found at: https://cmip6.science.unimelb.edu.au/search and https://aims2.llnl.gov/search/cmip6/

The framework requires each script to be run separately in a specific order (defined in [1], [2] and [3]), using the data described under the ‘Data’ section above. Our setup does not currently support other datasets than those described in Section 2 of the adjoining paper, however the setup does accomodate other CMIP6 models and their ensembles than those currently used to run the setup. The following paragraph indicates the order in which each script can be run to reimplement the method in the paper:

[1] All files ending in ‘parameterisation.m’, are the parameterisations that make up the emulator and are referred to in the manuscript. These scripts must be run first, in any order.
- CMIP6_Arctic_Amplification_Parameterisation.m (Section 2.3)
- AMST_parameterisation.m (Section 2.4)
- SIA_max_Parameterisation.m (Section 2.5.1)
- SIA_CMIP6_Parameterisation.m (Section 2.5.1)

[2] To calibrate the above parameterisations with the CMIP6 data, the following scripts must be run in the following order:
- Arctic_Seasonal_Temperature_Calibration.m
- SIA_max_Calibration.m
- SIA_Calibration.m

While this repository intends to provide a series of scripts to re-implement the method described in the Chilcott and Meinshausen, (2025) manuscript, (the scripts in [2] therefore intend to show the calibration of the first ensemble member of each CMIP6 model used), it is possible to use these scripts to calibrate to other ensemble members of the same CMIP6 model. 

[3] The calibration parameters generated from running the scripts in [2] are then used in the following scripts to constrain the CMIP6 calibrations to observations:
- Arctic_Amplification_Observational_Constraint.m (Section 2.3.1)
- Observationally_Constrained_Emulator_bias_corrections.m (Section 2.4.1 and Section 2.5.2)

[5] We evaluated our model performance to understand if the parameterisations provided could project the non-linearity of Arctic sea ice loss outside of the calibration period using the following scripts:
- Assessing_calibrations_to_2300.m (Section 3.1)

[4] We applied the emulator's parameterisation framework to questions in the sea ice discourse to further understand future sea ice loss (calulating both the probability of an ice-free Arctic Ocean and the remaining carbon budget to prevent seasonally ice-free conditions) using the following scripts:
- Probability_Calculations.m (Section 3.3)
- Carbon_Budget_Calculations.m (Section 3.3)

[6] Appendix.m: This file mainly generates the figures and tables in the supplementary material.

You can then run the both the CMIP6 parameterisation framwork of the emulator, and the observationally constrained framework.

=================================================================================
