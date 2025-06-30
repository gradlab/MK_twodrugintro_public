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

seq <- run_stochastic_brief_6_27(model  = sequential_stochastic_6_7_brief, 
                                 parms = baseline_parms_6_12, runs = 1000, dt = seq(0,365*50, 1), newfoldername = "baselineparms_1000runs_seq")
