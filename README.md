#Code supporting spatial and fitness modelling analyses of the paper "The spread of sexually transmissible drug-resistant shigellosis"

##Data
You can find data used in the analysis in the "Data" folder

##Run analysis
The code used to generate Figure 1 and 2 are available in "Figure1" and "Figure2" files. To run the modelling analysis, you must run the code in the following order:

###For relative groups model: 
- MCMCDataPrep_RelGroups.R
- MCMC_RelGroups.R

###For all other models:
- ModelDataPrep_AZ.R
- MCMC_IndependentAZ.R/MCMC_Average.R/MCMC_Sharedf0.R
