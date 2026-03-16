#Code supporting spatial and fitness modelling analyses of the paper "The spread of sexually transmissible drug-resistant shigellosis"

##Data
You can find data used in the analysis in the "Data" folder

##Run analysis
The code used to generate Figure 1 and 2 are available in "Figure1" and "Figure2" files. 
To run the modelling analysis, you must run the code in the following order:

###For relative groups model: 
- MCMCDataPrep_RelGroups.R
- MCMC_RelGroups.R

###For all other models:
- ModelDataPrep_AZ.R
- MCMC_IndependentAZ.R/MCMC_Average.R/MCMC_Sharedf0.R

###Supplemental Analyses 
1. For relative groups model: 
  - MCMCDataPrep_RelGroups.R (Uncomment subsampling commands)
  - MCMC_Sharedf0_Supp.R 
2. Analyses on CRO and CIP resistance
  - Prepare data: ModelDataPrep_Supp_CIP_CRO.R
  - For non-switch analyses (Average, Independent Models), read above data frames into MCMC_IndependentAZ.R/MCMC_Average.R)
  - For switch analyses: MCMC_Sharedf0_Supp.R
