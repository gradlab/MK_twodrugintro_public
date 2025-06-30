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

baseline_1perc_parms <- baseline_parms_6_12
baseline_10perc_parms <- baseline_parms_6_12
baseline_1perc_parms$threshold <- 0.01
baseline_10perc_parms$threshold <- 0.10
seq_1perc_ <- run_stochastic_brief_6_27(model = sequential_stochastic_6_7_brief, parms = baseline_1perc_parms, runs = 1000, dt = seq(0,365*50, 1), newfoldername = "seq_1perc_threshold")
seq_10perc <- run_stochastic_brief_6_27(model = sequential_stochastic_6_7_brief, parms = baseline_10perc_parms, runs = 1000, dt = seq(0,365*50, 1), newfoldername = "seq_10erc_threshold")
ea_1perc <- run_stochastic_brief_6_27(model = equal_allocation_stochastic_6_7_brief, parms = baseline_1perc_parms, runs = 1000, dt = seq(0,365*50, 1), newfoldername = "ea_1perc_threshold")
ea_10perc_ <- run_stochastic_brief_6_27(model = equal_allocation_stochastic_6_7_brief, parms = baseline_10perc_parms, runs = 1000, dt = seq(0,365*50, 1), newfoldername = "ea1_10perc_threshold")



