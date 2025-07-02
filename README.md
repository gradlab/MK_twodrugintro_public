This is the code used for the analysis and figures in the manuscript "Comparing Strategies to Introduce Two New Antibiotics for Gonorrhea: A Modeling Study" by Kline et al, currently in submision and available as a preprint at https://www.medrxiv.org/content/10.1101/2025.07.01.25330638v1. 

Code was run on the FASRC Cannon High Performance Computing cluster at Harvard University using R 4.4.1.

__Code__: \
The `code/` folder contains the code used to run the analyses shown in the manuscript. Files that end in "`.R`" are the R scripts that were run and files ending in "`.sh`" are the bash scripts used to run the R scripts that correspond in name. The "`logs/`" subfolder contains the log files from the bash scripts. The scripts correspond to the following analyses: 
* Calibration files:
  - `calibration_runscript.R`: Runs a drug A only model to run calibration, which determines values for baseline parameter set found in  `3drug_parameters_6_12.R`. This script writes to output files to `output/calibration/`:
    - `6_23_cluster_calibration_expparms.csv`: which exponentiates the dataframe of calibrated parameters, giving interpretable calibrated parameters for rate parameters 
    - `6_23_cluster_calibration_ilogitparms.csv`: which takes the inverse logit of the dataframe of calibrated parameters, giving interpretable calibrated parameters for probability parameters <br><br>
* Helper files:
  - `3drug_stochastic_models.R`: Contains the Odin models used in these analyses. 
  - `3drug_parameters_6_12.R`: Contains the baseline set of parameters and initial conditions.
  - `3drug_stochastic_functions.R`: Contains the helped functions coded to run the analyses. <br><br>
* Running determinstic models:
  - `deterministic_runscript.R`: runs the deterministic versions of the equal allocation (`equal_allocation_simplified_discrete_5_28`) and sequential (`sequential_simplified_discrete_5_28`) from `3drug_stochastic_models.R` using function `run_stochastic_savebrief_6_3` from `3drug_stochastic_functions.R` and produces output and figures.
  - Outputs are written to `output/deterministic_output/` and include:
    - `equal_allocation_prevs_deterministic.csv`
    - `equal_allocation_compartments_deterministic.csv`
    - `sequential_prevs_deterministic.csv`
    - `sequential_compartments_deterministic.csv`
  - Figures are written to `figures/deterministic_figures/` and together form Supplementary Figure 2.  <br><br>
* Files running stochastic baseline parameter simulations: 
  - `baselineparms_1000runs_ea.R`: Code used to run odin model `equal_allocation_stochastic_6_7_brief` from `3drug_stochastic_models.R` using run function `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` for 1000 simulations.
  - `baselineparms_1000runs_seq.R`: Code used to run odin model `sequential_stochastic_6_7_brief` from `3drug_stochastic_models.R` using run function `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` for 1000 simulations. <br><br>
  - `threshold_analysis.R` uses  `equal_allocation_stochastic_6_7_brief` and `sequential_stochastic_6_7_brief` from `3drug_stochastic_models.R` and `run_stochastic_brief_6_27` from `3drug_stochastic_functions.R` to run baseline parameter models with resistance prevalence thresholds of 1% or 10% instead of default 5% for 1000 simulations for each strategy.  <br><br>
* Files for processing baseline parameter simulation results into summarized outputs and making figures from them:  
  - `baselineparms_1000runs_output.R`: takes the output generated from `baselineparms_1000runs_ea.R` and `baselineparms_1000runs_seq.R`, which was too large to include in the GitHub repo, and summarizes it into files written to `output/baseline_parms_analysis/`:
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
  - `baselineparms_1000runs_plots.R`: takes in the outputs generated in `baselineparms_1000runs_output.R` and generates figures. Figures generated are written to `figures/baseline_parms_figures/` and include:
    - `baseline_prevplot_50y.tiff`: becomes Figure 1 panel A
    - `ea_propres_heatmap.tif` and `seq_propres_heatmap.tif`: become Figure 1 panel B
    - `baseline_avglost_50y.tiff`: becomes Supplementary Figure 3 panel A
    - `ea_drugsremaining_heatmap.tiff` and `seq_drugsremaining_heatmap.tiff`: become Supplementary Figure 3 panel B
    - `baseline_all3lost_50y.tiff`: becomes Supplementary Figure 4 panel A
    - `alldrugslost_heatmap_ea_yrs.tiff` and `alldrugslost_heatmap_seq_yrs.tiff`: become Supplementary Figure 4 panel B
  - `threshold_analysis_output_plots.R`: takes output generated from `threshold_analysis.R`, which was too large to include in the GitHub repo, and creates summarized outputs and plots from them. The output files from this script are written to `output/baseline_parms_analysis/` and include:
    - `seq_1perc_locksumdf.csv`
    - `seq_5perc_locksumdf.csv`
    - `seq_10perc_locksumdf.csv`
    - `ea_1perc_locksumdf.csv`
    - `ea_5perc_locksumdf.csv`
    - `ea_10perc_locksumdf.csv`
  - The plots that are generated from `threshold_analysis_output_plots.R` include the components of Figure 2, written to `/figures/baseline_parms_figures/`:
    - `threshold_compare_plot.tiff` which becomes panel A
    - `threshold_heatmap_ea.tiff` and `threshold_heatmap_seq.tiff` which become panel B.
* Sensitivity Analyses:
  - `make_LHS_draws.R`: Draws two sets of 1000 sets of parameter combinations. The first is of Group 1 parameters, to which we did not expect our results to be sensitive. The dataframe of these parameter combinations is written to `output/6_7_unsensitive_LHS_df.csv`. The second is of Group 2 parameters, which we expected that our results could be sensitive to. The dataframe of these parameter combinations is written to `output/6_11_sensitive_LHS_df_fixed.csv`
  - Sensitivity analysis of Group 1 parameters (termed informally "unsensitivity analysis"):
    - `unsensitivity_analysis_ea.R`: Uses function `run_unsens_analysis` from `3drug_stochastic_functions.R` to run equal allocation stochastic model `equal_allocation_stochastic_6_7_brief` from `3drug_stochastic_models.R` for each set of parameters from `output/6_7_unsensitive_LHS_df.csv`. For each parameter combination under the equal allocation strategy, it runs 100 simulations and saves only datapoints from the times specified. 
    - `unsensitivity_analysis_seq.R`: Uses function `run_unsens_analysis` from `3drug_stochastic_functions.R` to run equal allocation stochastic model `sequential_stochastic_6_7_brief` from `3drug_stochastic_models.R` for each set of parameters from `output/6_7_unsensitive_LHS_df.csv`. For each parameter combination under the sequential strategy, it runs 100 simulations and saves only datapoints from the times specified.
    - `unsensitivity_analysis_output_plots.R`:  Reads in results from `unsensitivity_analysis_ea.R` and `unsensitivity_analysis_seq.R`, which are too large to store in github, summarizes the results, and makes plots.
        - The outputs are written to `output/unsensitivity_analysis_outputs/` and include:
            - `seqfulldf.csv`: The proportion of the 100 sequential simulations for a given LHS parameter combination from `output/6_7_unsensitive_LHS_df.csv` that had already hit the resistance prevalence threshold at that time point. 
            - `eafulldf.csv`: The proportion of the 100 equal allocation simulations for a given LHS parameter combination from `output/6_7_unsensitive_LHS_df.csv` that had already hit the resistance prevalence threshold at that time point. 
            -  `diffdfall.csv`: The proportion of the 100 equal allocation simulations for a given LHS parameter combination from `output/6_7_unsensitive_LHS_df.csv` that had already hit the resistance prevalence threshold at that time point minus the proportion of the 100 sequential simulations for a given LHS parameter combination from `output/6_7_unsensitive_LHS_df.csv` that had already hit the resistance prevalence threshold at that time point. Used to generate Supplementary Figure 5. 
        -  The figures are written to `figures/2025-06-30_unsensitivity_figures` and include `diffmap_2yr.tiff`, `diffmap_5yr.tiff`, and `diffmap_10yr.tiff` which form Supplementary Figure 5.
    -  `seqsens_runscript.R`: Uses function `run_sens_analysis_6_13` from `3drug_stochastic_functions.R` to run equal allocation stochastic model `sequential_stochastic_6_7_brief` from `3drug_stochastic_models.R` for each set of parameters from `output/6_11_sensitive_LHS_df_fixed.csv`. For each parameter combination under the sequential strategy, it runs 100 simulations and saves only datapoints from the times specified. 
    -  `ea_sens_runscript.R`: Uses function `run_sens_analysis_6_13` from `3drug_stochastic_functions.R` to run equal allocation stochastic model `equal_allocation_stochastic_6_7_brief` from `3drug_stochastic_models.R` for each set of parameters from `output/6_11_sensitive_LHS_df_fixed.csv`. For each parameter combination under the equal allocation strategy, it runs 100 simulations and saves only datapoints from the times specified.
    -  `sens_anal_plots.R`: takes in outputs from `seqsens_runscript.R` and `ea_sens_runscript.R`, which were too large to upload to github, and creates plots and outputs which include:
        - `output/6_23_edgecase_sens_analy_all.csv` which includes the list of 22 edge cases from the 1000 sets of parameters. This becomes Supplementary Table 4.
        - `figures/sens_analysis_figures/drugslost_scatter_5yrs.tiff` which becomes Figure 3 panel A
        - `figures/sens_analysis_figures/drugslost_scatter_10yrs.tiff` which becomes Figure 3 panel B
        - `figures/sens_analysis_figures/drugslost_scatter_20yrs.tiff` which becomes Figure 3 panel C
        - `figures/sens_analysis_figures/drugB_prevcompare_scatter_10yrs_all.tiff` which becomes Figure 3 panel D
        - `figures/sens_analysis_figures/drugC_prevcompare_scatter_10yrs_all.tiff` which becomes Figure 3 panel E
        -   <br><br>
        

      

  



