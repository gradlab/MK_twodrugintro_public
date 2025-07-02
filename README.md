This is the code used for the analysis and figures in the manuscript "Comparing Strategies to Introduce Two New Antibiotics for Gonorrhea: A Modeling Study" by Kline et al, currently in submision.

__Code__: \
The `code/` folder contains the code used to run the analyses shown in the manuscript. Files that end in "`.R`" are the R scripts that were run and files ending in "`.sh`" are the bash scripts used to run the R scripts that correspond in name. The "`logs/`" subfolder contains the log files from the bash scripts. The scripts correspond to the following analyses: 
* Helper files: \
      `3drug_stochastic_models.R`: Contains the Odin models used in these analyses. \
      `3drug_parameters_6_12.R`: Contains the baseline set of parameters and initial conditions. \
      `3drug_stochastic_functions.R`: Contains the helped functions coded to run the analyses.
* Files running baseline parameter simulations: \
      `baselineparms_1000runs_ea.R`: Code used to run odin model `equal_allocation_stochastic_6_7_brief` from `3drug_stochastic_models.R` using run function `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` for 1000 simulations.
      `baselineparms_1000runs_seq.R`: Code used to run odin model `sequential_stochastic_6_7_brief` from `3drug_stochastic_models.R` using run function `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` for 1000 simulations.
