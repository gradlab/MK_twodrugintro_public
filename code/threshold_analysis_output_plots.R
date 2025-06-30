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


seq1perc_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-27seq_1perc_threshold/lockdf.csv")
seq5perc_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-28baselineparms_1000runs_seq/lockdf.csv")
seq10perc_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-27seq_10erc_threshold/lockdf.csv")

ea1perc_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-27ea_1perc_threshold/lockdf.csv")
ea5perc_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-28baselineparms_1000runs_ea/lockdf.csv")
ea10perc_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-27ea1_10perc_threshold/lockdf.csv")

seq1000_1perc_locksum_df <- make_lock_df(seq1perc_lockdf, times = seq(0,365*50, 1)) |> mutate(strategy = "sequential", threshold = 0.01)
write_csv(seq1000_1perc_locksum_df, "../output/baseline_parms_analysis/seq_1perc_locksumdf.csv")
seq1000_5perc_locksum_df <- make_lock_df(seq5perc_lockdf, times = seq(0,365*150, 1)) |> mutate(strategy = "sequential", threshold = 0.05)
write_csv(seq1000_5perc_locksum_df, "../output/baseline_parms_analysis/seq_5perc_locksumdf.csv")
seq1000_10perc_locksum_df <- make_lock_df(seq10perc_lockdf, times = seq(0,365*150, 1)) |> mutate(strategy = "sequential", threshold = 0.10)
write_csv(seq1000_10perc_locksum_df, "../output/baseline_parms_analysis/seq_10perc_locksumdf.csv")

ea1000_1perc_locksum_df <- make_lock_df(ea1perc_lockdf, times = seq(0,365*50, 1)) |> mutate(strategy = "equal allocation", threshold = 0.01)
write_csv(ea1000_1perc_locksum_df, "../output/baseline_parms_analysis/ea_1perc_locksumdf.csv")
ea1000_5perc_locksum_df <- make_lock_df(ea5perc_lockdf, times = seq(0,365*150, 1)) |> mutate(strategy = "equal allocation", threshold = 0.05)
write_csv(ea1000_5perc_locksum_df, "../output/baseline_parms_analysis/ea_5perc_locksumdf.csv")
ea1000_10perc_locksum_df <- make_lock_df(ea10perc_lockdf, times = seq(0,365*150, 1)) |> mutate(strategy = "equal allocation", threshold = 0.10)
write_csv(ea1000_10perc_locksum_df, "../output/baseline_parms_analysis/ea_10perc_locksumdf.csv")

avg_drugs_lost_threshold_50yrs <- rbind(seq1000_1perc_locksum_df, ea1000_1perc_locksum_df) |> rbind(seq1000_5perc_locksum_df) |>
  rbind(ea1000_5perc_locksum_df) |>
  rbind(seq1000_10perc_locksum_df) |> rbind(ea1000_10perc_locksum_df) |>
  mutate(Years = step / 365) |>
  filter(Years <= 50) |>
  mutate(threshold_strategy = paste(threshold, strategy)) |>
  ggplot(aes(x = Years, y = avg_drugs_lost, group = interaction(strategy, as.factor(threshold)), color = strategy, linetype = as.factor(threshold))) +
  geom_line() + theme_classic() +
  scale_color_manual(name = "Strategy", values = c("sequential" = '#08519c',"equal allocation" = '#d73027'),
                     labels = c("sequential" = "Sequential", "equal allocation" = "Equal Allocation")) +
  scale_linetype_manual(name = "Threshold", values = c("0.01" =	3, "0.05" = 1, "0.1" = 5), labels = c("0.01" = "1%", "0.05" = "5%", "0.1" = "10%")) +
  ylim(0,3) + ylab("Average Drugs Lost")

ggsave("../figures/baseline_parms_figures/threshold_compare_plot.tiff", avg_drugs_lost_threshold_50yrs, width = 6, height = 4,
       dpi = 900)




ea_heatmap_threshold <- rbind(ea1000_5perc_locksum_df, ea1000_1perc_locksum_df) |> rbind(ea1000_10perc_locksum_df) |>
  mutate(Years = step / 365) |>
  filter(Years %in% c(2, 5, 10, 30)) |>
  ggplot(aes(x = as.factor(threshold), y = as.factor(Years), fill = avg_drugs_lost)) + geom_tile(color = "white") +
  geom_text(aes(label = avg_drugs_lost), color = "black", size = 3) +
  scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 3), name = "Avg Drugs Lost") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Years") + xlab("Threshold") +
  scale_x_discrete(labels = c("1%", "5%", "10%"))

seq_heatmap_threshold <- rbind(seq1000_5perc_locksum_df, seq1000_1perc_locksum_df) |> rbind(ea1000_10perc_locksum_df) |>
  mutate(Years = step / 365) |>
  filter(Years %in% c(2, 5, 10, 30)) |>
  ggplot(aes(x = as.factor(threshold), y = as.factor(Years), fill = avg_drugs_lost)) + geom_tile(color = "white") +
  geom_text(aes(label = avg_drugs_lost), color = "black", size = 3) +
  scale_fill_gradient(low = "#deebf7", high = '#08519c', limits = c(0, 3), name = "Avg Drugs Lost") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Years") + xlab("Threshold") +
  scale_x_discrete(labels = c("1%", "5%", "10%"))

ggsave("../figures/baseline_parms_figures/threshold_heatmap_ea.tiff", ea_heatmap_threshold, width = 3.5, height = 5, dpi = 900 )
ggsave("../figures/baseline_parms_figures/threshold_heatmap_seq.tiff", seq_heatmap_threshold, width = 3.5, height = 5, dpi = 900)





