This is the code used for the analysis and figures in the manuscript "Comparing Strategies to Introduce Two New Antibiotics for Gonorrhea: A Modeling Study" by Kline et al, currently in submision and available as a preprint at https://www.medrxiv.org/content/10.1101/2025.07.01.25330638v1. 

__Code__: \
The `code/` folder contains the code used to run the analyses shown in the manuscript. Files that end in "`.R`" are the R scripts that were run and files ending in "`.sh`" are the bash scripts used to run the R scripts that correspond in name. The "`logs/`" subfolder contains the log files from the bash scripts. The scripts correspond to the following analyses: 
* Helper files: 
     - `3drug_stochastic_models.R`: Contains the Odin models used in these analyses. 
     - `3drug_parameters_6_12.R`: Contains the baseline set of parameters and initial conditions. 
     - `3drug_stochastic_functions.R`: Contains the helped functions coded to run the analyses. <br><br>
* Files running baseline parameter simulations: 
     - `baselineparms_1000runs_ea.R`: Code used to run odin model `equal_allocation_stochastic_6_7_brief` from `3drug_stochastic_models.R` using run function `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` for 1000 simulations.
     - `baselineparms_1000runs_seq.R`: Code used to run odin model `sequential_stochastic_6_7_brief` from `3drug_stochastic_models.R` using run function `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` for 1000 simulations. <br><br>
     - `threshold_analysis.R` uses  `equal_allocation_stochastic_6_7_brief` and `sequential_stochastic_6_7_brief` from `3drug_stochastic_models.R` and `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` to run baseline parameter models with resistance prevalence thresholds of 1% or 10% instead of default 5% for 1000 simulations for each strategy.  <br><br>
* Files for processing baseline parmeter simulation results into summarized outputs and making figures from them: 
     - `baselineparms_1000runs_output.R`: takes the output generated from `baselineparms_1000runs_ea.R` and `baselineparms_1000runs_seq.R`, which was too large to include in the github repo, and summarizes it into files written to `output/baseline_parms_analysis/`:
          - `ea_baselineparms_avglost.csv`
          - `ea_baselineparms_lostall3.csv`
          - `ea_baselineparms_numslost.csv`
          - `ea_baselineparms_propdf.csv`
          - `seq_baselineparms_avglost.csv`
          - `seq_baselineparms_lostall3.csv`
          - `seq_baselineparms_numslost.csv`
          - `seq_baselineparms_propdf.csv`
          - `combined_baselineparms_numslost.csv`
          - `combined_baselineparms_numslost_filtered.csv`: becomes Supplementary Table 3, and the other outputs are used for generating main and supplementary figures, as described.
    - `threshold_analysis_output_plots.R `: takes output generated from `threshold_analysis.R`, which were too large to include in the github repo, and creates summarized outputs and plots from them. The output files from this script are written to `output/baseline_parms_analysis/` and include: \
              - `seq_1perc_locksumdf.csv`
      

  



