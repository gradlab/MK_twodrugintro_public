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

seq1000_baselineparms_prop_df <- read_csv("../output/baseline_parms_analysis/seq_baselineparms_propdf.csv")
seq1000_baselineparms_numlost_df <- read_csv("../output/baseline_parms_analysis/seq_baselineparms_avglost.csv")
seq1000_baselineparms_all3lost_df <- read_csv("../output/baseline_parms_analysis/seq_baselineparms_lostall3.csv")

ea1000_baselineparms_prop_df <- read_csv("../output/baseline_parms_analysis/ea_baselineparms_propdf.csv")
ea1000_baselineparms_numlost_df <- read_csv("../output/baseline_parms_analysis/ea_baselineparms_avglost.csv")
ea1000_baselineparms_all3lost_df <- read_csv("../output/baseline_parms_analysis/ea_baselineparms_lostall3.csv")

propplot_1000sims_linetype_50yr <- rbind(seq1000_baselineparms_prop_df, ea1000_baselineparms_prop_df) |>
  filter(Years <= 50) |>
  pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
  mutate(full_name = paste(strategy, name)) |>
  ggplot(aes(x = Years, y = value, color = full_name, linetype = full_name)) + geom_line() +
  theme_classic() +
  scale_color_manual(values = c("equal allocation prop_runs_prevA_above" = '#d73027',
                                "sequential prop_runs_prevA_above" = '#041c5a',
                                "equal allocation prop_runs_prevB_above" = '#fc8d59',
                                "sequential prop_runs_prevB_above" = '#2171b5',
                                "equal allocation prop_runs_prevC_above" = '#fdb863',
                                "sequential prop_runs_prevC_above" = '#6baed6'),
                     labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A",
                                "Sequential, Drug B", "Sequential, Drug C"),
                     name = "Strategy and Drug Resistance Prevalence") +
  scale_linetype_manual(values = c("equal allocation prop_runs_prevA_above" = 1,
                                   "sequential prop_runs_prevA_above" = 2,
                                   "equal allocation prop_runs_prevB_above" = 1,
                                   "sequential prop_runs_prevB_above" = 2,
                                   "equal allocation prop_runs_prevC_above" = 1,
                                   "sequential prop_runs_prevC_above" = 2),
                        labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A",
                                   "Sequential, Drug B", "Sequential, Drug C"),
                        name = "Strategy and Drug Resistance Prevalence") +
  ylab("Proportion of simulations")

ggsave("../figures/baseline_parms_figures/baseline_prevplot_lt_50y.tiff", propplot_1000sims_linetype_50yr, width = 8, height = 4,
       dpi = 900)

propplot_1000sims_50yr <- rbind(seq1000_baselineparms_prop_df, ea1000_baselineparms_prop_df) |>
  filter(Years <= 50) |>
  pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
  mutate(full_name = paste(strategy, name)) |>
  ggplot(aes(x = Years, y = value, color = full_name, linetype = full_name)) + geom_line() +
  theme_classic() +
  scale_color_manual(values = c("equal allocation prop_runs_prevA_above" = '#d73027',
                                "sequential prop_runs_prevA_above" = '#041c5a',
                                "equal allocation prop_runs_prevB_above" = '#fc8d59',
                                "sequential prop_runs_prevB_above" = '#2171b5',
                                "equal allocation prop_runs_prevC_above" = '#fdb863',
                                "sequential prop_runs_prevC_above" = '#6baed6'),
                     labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A",
                                "Sequential, Drug B", "Sequential, Drug C"),
                     name = "Strategy and Drug Resistance Prevalence") +
  scale_linetype_manual(values = c("equal allocation prop_runs_prevA_above" = 1,
                                   "sequential prop_runs_prevA_above" = 1,
                                   "equal allocation prop_runs_prevB_above" = 1,
                                   "sequential prop_runs_prevB_above" = 1,
                                   "equal allocation prop_runs_prevC_above" = 1,
                                   "sequential prop_runs_prevC_above" = 1),
                        labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A",
                                   "Sequential, Drug B", "Sequential, Drug C"),
                        name = "Strategy and Drug Resistance Prevalence") +
  ylab("Proportion of simulations")

ggsave("../figures/baseline_parms_figures/baseline_prevplot_50y.tiff", propplot_1000sims_50yr, width = 8, height = 4,
       dpi = 900)

#make heatmaps
ea_heatmap1000sims_years <- ea1000_baselineparms_prop_df |>
  filter(step %in% c(730, 1825, 3650, 10950)) |>
  pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
  ggplot(aes(x = name, y = as.factor(Years), fill = value)) + geom_tile(color = "white") +
  geom_text(aes(label = value), color = "black", size = 3) +
  scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Proportion of runs") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
  ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("Years") + xlab("")

ggsave("../figures/baseline_parms_figures/ea_propres_heatmap.tiff", ea_heatmap1000sims_years, width = 3.5, height = 5,
       dpi = 900)

seq_heatmap1000sims_years <- seq1000_baselineparms_prop_df |>
  filter(step %in% c(730, 1825, 3650, 10950)) |>
  pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
  ggplot(aes(x = name, y = as.factor(Years), fill = value)) + geom_tile(color = "white") +
  geom_text(aes(label = value), color = "black", size = 3) +
  scale_fill_gradient(low = "#deebf7", high = '#08519c', limits = c(0, 1), name = "Proportion of runs") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
  ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("Years") + xlab("")

ggsave("../figures/baseline_parms_figures/seq_propres_heatmap.tiff", seq_heatmap1000sims_years, width = 3.5, height = 5,
       dpi = 900)


#avg drugs lost plots
avg_drugs_lost_plot_years_linetype_50yrs <- rbind(seq1000_baselineparms_numlost_df, ea1000_baselineparms_numlost_df) |>
  filter(Years <= 50) |>
  ggplot(aes(x = Years, y = avg_drugs_lost, color = strategy, linetype = strategy)) +
  geom_line() + theme_classic() + ylim(0,3) + scale_color_manual(name = "Strategy", values = c("sequential" = '#08519c',
                                                                                               "equal allocation" = '#d73027'),
                                                                 labels = c("sequential" = "Sequential",
                                                                            "equal allocation" = "Equal allocation")) +
  scale_linetype_manual(name = "Strategy", values = c("sequential" = 2, "equal allocation" = 1),
                        labels = c("sequential" = "Sequential",  "equal allocation" = "Equal allocation")) +
  ylab("Average Drugs Lost")

ggsave("../figures/baseline_parms_figures/baseline_avglost_lt_50y.tiff", avg_drugs_lost_plot_years_linetype_50yrs, width = 6, height = 4,
       dpi = 900)

avg_drugs_lost_plot_years_50yrs <- rbind(seq1000_baselineparms_numlost_df, ea1000_baselineparms_numlost_df) |>
  filter(Years <= 50) |>
  ggplot(aes(x = Years, y = avg_drugs_lost, color = strategy, linetype = strategy)) +
  geom_line() + theme_classic() + ylim(0,3) + scale_color_manual(name = "Strategy", values = c("sequential" = '#08519c',
                                                                                               "equal allocation" = '#d73027'),
                                                                 labels = c("sequential" = "Sequential",
                                                                            "equal allocation" = "Equal allocation")) +
  scale_linetype_manual(name = "Strategy", values = c("sequential" = 1, "equal allocation" = 1),
                        labels = c("sequential" = "Sequential",  "equal allocation" = "Equal allocation")) +
  ylab("Average Drugs Lost")

ggsave("../figures/baseline_parms_figures/baseline_avglost_50y.tiff", avg_drugs_lost_plot_years_50yrs, width = 6, height = 4,
       dpi = 900)


#make avg drugs lost heatmaps
ea_drugs_remaining_plot_yrs <- ea1000_baselineparms_numlost_df  |>
  filter(Years %in% c(2, 5, 10, 30)) |>
  ggplot(aes(x = strategy, y = as.factor(Years), fill = avg_drugs_lost)) + geom_tile(color = "white") + 
  geom_text(aes(label = avg_drugs_lost), color = "black", size = 3) +
  scale_fill_gradient(low = '#fde0dd', high = "#d73027", limits = c(0, 3), name = "Average Drugs Lost") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = c("Equal Allocation")) + ylab("Years") + xlab("")

ggsave("../figures/baseline_parms_figures/ea_drugsremaining_heatmap.tiff", ea_drugs_remaining_plot_yrs, width = 3.5, height = 4, dpi = 900)

seq_drugs_remaining_plot_yrs <- seq1000_baselineparms_numlost_df |>
  filter(Years %in% c(2, 5, 10, 30)) |>
  ggplot(aes(x = strategy, y = as.factor(Years), fill =avg_drugs_lost)) + geom_tile(color = "white") + 
  geom_text(aes(label = avg_drugs_lost), color = "black", size = 3) +
  scale_fill_gradient(low = '#deebf7', high = "#08519c", limits = c(0, 3), name = "Average Drugs Lost") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(labels = c("Sequential")) + ylab("Years") + xlab("")


ggsave("../figures/baseline_parms_figures/seq_drugsremaining_heatmap.tiff", seq_drugs_remaining_plot_yrs, width = 3.5, height = 4, dpi = 900)

#make all 3 remaining plots

alldrugslost_plot_linetype_50y <- rbind(seq1000_baselineparms_all3lost_df, ea1000_baselineparms_all3lost_df) |>
  filter(Years <= 50) |>
  ggplot(aes(x = Years, y = avg_lost3, color = strategy, linetype = strategy)) +
  geom_line() + theme_classic() + ylim(0,1) + scale_color_manual(values = c("sequential" = '#08519c',
                                                                            "equal allocation" = '#d73027')) +
  ylab("Proportion of runs with all 3 drugs lost")

ggsave("../figures/baseline_parms_figures/baseline_all3lost_lt_50y.tiff", alldrugslost_plot_linetype_50y, width = 6, height = 4,
       dpi = 900)


alldrugslost_plot_50y <- rbind(seq1000_baselineparms_all3lost_df, ea1000_baselineparms_all3lost_df) |>
  filter(Years <= 50) |>
  ggplot(aes(x = Years, y = avg_lost3, color = strategy)) +
  geom_line() + theme_classic() + ylim(0,1) + scale_color_manual(values = c("sequential" = '#08519c',
                                                                            "equal allocation" = '#d73027')) +
  ylab("Proportion of runs with all 3 drugs lost")

ggsave("../figures/baseline_parms_figures/baseline_all3lost_50y.tiff", alldrugslost_plot_50y, width = 6, height = 4,
       dpi = 900)

avglost3_ea_yrs <- ea1000_baselineparms_all3lost_df |>
  filter(Years %in% c(10, 30, 50)) |>
  ggplot(aes(x = strategy, y = as.factor(Years), fill = avg_lost3)) + geom_tile(color = "white") +
  geom_text(aes(label = avg_lost3), color = "black", size = 3) +
  scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Average Runs") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("Equal Allocation")) + ylab("Years") + xlab("")

avglost3_seq_yrs <- seq1000_baselineparms_all3lost_df |>
  filter(Years %in% c(10, 30, 50)) |>
  ggplot(aes(x = strategy, y = as.factor(Years), fill = avg_lost3)) + geom_tile(color = "white") +
  geom_text(aes(label = avg_lost3), color = "black", size = 3) +
  scale_fill_gradient(low = "#deebf7", high = '#08519c', limits = c(0, 1), name = "Average Runs") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c("Sequential")) + ylab("Years") + xlab("")

ggsave("../figures/baseline_parms_figures/alldrugslost_heatmap_ea_yrs.tiff", avglost3_ea_yrs, width = 3, height = 4,
       dpi = 900)
ggsave("../figures/baseline_parms_figures/alldrugslost_heatmap_seq_yrs.tiff", avglost3_seq_yrs, width = 3, height = 4,
       dpi = 900)

