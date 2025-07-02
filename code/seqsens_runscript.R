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

sensitive_lhs_filled <- read_csv("../output/6_11_sensitive_LHS_df_fixed.csv")

run_sens_analysis_6_13(model = sequential_stochastic_6_7_brief, strategy = "sequential", threshold = 0.05,
                  paramdf = sensitive_lhs_filled, runs = 100, years = 28, times = c(730, 1000, 1825, 3650, 5000, 7300, 10000),
                  newfoldername = "seq_sens_fixed_613", prefix = "seq_sens")
