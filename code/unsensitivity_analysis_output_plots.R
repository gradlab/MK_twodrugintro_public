#load necessary packages
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(odin)
library(grid)

source("3drug_stochastic_functions.R")

#load in output files from unsensitivity analysis
seqA <- read_csv("../output/2025-06-29seq_unsens/seq_unsenssequentialprevAdf.csv")
seqB <- read_csv("../output/2025-06-29seq_unsens/seq_unsenssequentialprevBdf.csv")
seqC <- read_csv("../output/2025-06-29seq_unsens/seq_unsenssequentialprevCdf.csv")

eaA <- read_csv("../output/2025-06-29ea_unsens/ea_unsensequal allocationprevAdf.csv")
eaB <- read_csv("../output/2025-06-29ea_unsens/ea_unsensequal allocationprevBdf.csv")
eaC <- read_csv("../output/2025-06-29ea_unsens/ea_unsensequal allocationprevCdf.csv")

make_plots <- make_heat_maps_from_fxn_output_yrs(seqA, seqB, seqC, eaA, eaB, eaC, foldername = "_unsensitivity_figures")

#save outputs
seqfulldf <- make_plots[[1]]
eafulldf <- make_plots[[2]]
difdfall <- make_plots[[3]]

write_csv(seqfulldf ,"../output/unsensitivity_analysis_outputs/seqfulldf.csv")
write_csv(eafulldf ,"../output/unsensitivity_analysis_outputs/eafulldf.csv")
write_csv(difdfall ,"../output/unsensitivity_analysis_outputs/diffdfall.csv")

