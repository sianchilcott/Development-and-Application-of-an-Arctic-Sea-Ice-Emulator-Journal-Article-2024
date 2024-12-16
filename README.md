=================================================================================

SETUP INSTRUCTIONS FOR ZENODO AND GITHUB REPOSITORIES (https://doi.org/10.5281/zenodo.14020702):

The following provides an introduction and guide to the Zenodo and github repositories, for the pre-print manuscriupt presenting an Arctic Sea Ice Emulator. 

[1] CMIP6 SIA data between 1850 and 2100, alongside MAGICC global mean

[2]All files ending in ‘parameterisation.m’, are the parameterisations that make up the emulator and are referred to in the manuscript. These scripts must be run first, in any order. 

[3] The following scripts must be run in the following order:
- Arctic_Seasonal_Temperature_Calibration.m
- Arctic_Seasonal_Temperature_Calibration.m
- SIA_max_Calibration.m
- SIA_Calibration.m

These scripts firstly calibrate the parameterisations to CMIP6 data, the calibration parameters are then used in the following scripts to constrain the CMIP6 calibrations to observations:
- Arctic_Amplification_Observational_Constraint.m
- Observationally_Constrained_Emulator_bias_corrections.m

[4] To apply the emulator to questions in the sea ice discourse e.g calulating both the probability of an ice-free Arctic Ocean and the remaining carbon budget to prevent seasonally ice-free conditions use:
- Probability_Calculations.m
- Carbon_Budget_Calculations.m

[5] To test the calibration parameters in extended runs use:
- Assessing_calibrations_to_2300.m

[5] Appendix.m: This file mainly generates the figures and tables in the appendix.

You can then run the both the CMIP6 portion of the emulator, and the observationally constrained emulator.

=================================================================================
