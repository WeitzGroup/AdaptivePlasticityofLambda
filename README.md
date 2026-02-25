# _Code for_: The Adaptive Plasticity of Phage Lambda
_Authors: Sylvain Gandon, Tapan Goel, Joshua S. Weitz, Sebastian Lion_

Code has been archived on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18711589.svg)](https://doi.org/10.5281/zenodo.18711589)

## General description:
This repository contains all the code and data needed to generate Figures S2,S3 and S5 in the paper "The Adapative Plasticity of Phage Lambda". The code is written in MATLAB 2024a and MATLAB 2025b. 
The code associated with the entire mansucript can be found on the Zenodo archive.

## File descriptions:

1. **FigureS2.m** : Generates figure S2.
2. **FigureS3.m** : Generates figure S3.
3. **FigureS5.m** : Generates figure S5 using data from Data/ComparisonData.xlsx.

4. **Utils/ODE_SELV_MOI2.m** : Contains the ODE equations for the SELV model with coinfection for a maximal multiplicity of infection = 2   and any number of viral strains.
5. **Utils/MultiSpeciesTrajectories.m** : Generates population trajectory of a polymorphic population where one of the viral types is       the ESS and the others are sampled randomly. Also provides the mean and standard deviations of the trait values as functions of time. This code generates the data needed for Figure S3.
6. **Utils/PairInvasionDynamics.m** : Simulates pairwise invasions to a given resident strategy by user specified mutant strategies. Also compuates the initial invasion growth rate for each mutant strategy. This code generates data needed for Figure S2. Uses the MATLAB parallel computing toolbox.
7. **Utils/MutantGrowthRate.m** : Takes the resident-mutant timeseries and calculates the mutant invasion growth rate by fitting an exponential to the peaks of the mutant population fraction. 

8. **Data/ComparisonData.xlsx** : Contains data used to generate Figure S5.
9. **Data/FigureData.mat** : Contains data used to generate Figure S2. Currently not in use. The data for the figure is generated insitu using the PairInvasionDynamics function. 





