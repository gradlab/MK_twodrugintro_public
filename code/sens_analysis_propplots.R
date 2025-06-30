#load packages
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(odin)
library(grid)

#load any necessary helper functions
source("3drug_stochastic_models.R")
source("3drug_stochastic_functions.R")
source("3drug_parameters_6_12.R")

#load in list of parameter combinations to run:
parm_combs_df <- read_csv("../output/6_23_edgecase_sens_analy_all.csv")
parm_list <- parm_combs_df |> pull(parameter_comb) |> unlist()
sensitive_lhs_filled <- read_csv("../output/6_11_sensitive_LHS_df_fixed.csv")
for(c in parm_list){
  #run EA model by getting seeds
  EA_seed <- read_csv(paste0("../output/2025-06-14ea_sens_fixed_613/ea_sensseed_list_", c, ".csv")) |> unlist()
  ea_prop_df <- run_stochastic_brief_6_11(model = equal_allocation_stochastic_6_7_brief,
                                          parms = generate_parms_sens_6_13(sensitive_lhs_filled[c,], threshold = 0.05),
                                          seeds = EA_seed,
                                          runs = 100,
                                          threshold = 0.05,
                                          dt = seq(0,365*50, 1)) |> mutate(strategy = "equal allocation")
  seq_seed <- read_csv(paste0("../output/2025-06-14seq_sens_fixed_613/seq_sensseed_list_", c, ".csv")) |> unlist()
  seq_prop_df <- run_stochastic_brief_6_11(model = sequential_stochastic_6_7_brief,
                                           parms = generate_parms_sens_6_13(sensitive_lhs_filled[c,], threshold = 0.05),
                                           seeds = seq_seed,
                                           runs = 100,
                                           threshold = 0.05,
                                           dt = seq(0,365*50, 1)) |> mutate(strategy = "seq")

  comb_df <- rbind(ea_prop_df, seq_prop_df) |> pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above))
  write_csv(comb_df, paste0("../output/sens_anal_outputs/", c, "comb_prop_df.csv"))
  ggsave(paste0("../figures/sens_analysis_figures/prop_plots/", c, "propplot.tiff"),make_prop_plot_years_solid(comb_df, threshold = 0.05) + ggtitle(paste("Combination:", c)), width = 6, height = 3, dpi = 900)

}
