library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(odin)
library(grid)

#load any necessary helper functions
source("3drug_stochastic_models.R")
source("3drug_stochastic_functions.R")

unsensitive_lhs_filled <- read_csv("../output/6_7_unsensitive_LHS_df.csv")
run_unsens_analysis(model = sequential_stochastic_6_7_brief, strategy = "sequential", threshold = 0.05,
                       paramdf = unsensitive_lhs_filled, runs = 100, years = 28, times = c(730, 1000, 1825, 3650, 5000, 7300, 10000),
                       newfoldername = "seq_unsens", prefix = "seq_unsens")
