#load necessary packages
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


#run deterministic models under baseline parameters, make set of plots  
eadet_forplot <- run_stochastic_savebrief_6_3(model = equal_allocation_simplified_discrete_5_28, parms = baseline_parms_6_12, runs = "det",
                                              dt  =seq(0,365*150, 1))

seqdet_forplot <- run_stochastic_savebrief_6_3(model = sequential_simplified_discrete_5_28, parms = baseline_parms_6_12, runs = "det",
                                               dt  =seq(0,365*150, 1))

#save model output 
write_csv(eadet_forplot[[1]], "../output/deterministic_output/equal_allocation_prevs_deterministic.csv")
write_csv(eadet_forplot[[2]], "../output/deterministic_output/equal_allocation_compartments_deterministic.csv")
write_csv(seqdet_forplot[[1]], "../output/deterministic_output/sequential_prevs_deterministic.csv")
write_csv(seqdet_forplot[[2]], "../output/deterministic_output/sequential_compartments_deterministic.csv")


#make plots 
prevdf_ea <- eadet_forplot[[1]] |> mutate(strategy = "Equal Allocation")
prevdf_seq <- seqdet_forplot[[1]] |> mutate(strategy = "Sequential")

plotall_det_prevs_probs <- rbind(prevdf_ea, prevdf_seq) |> pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
  mutate(name_strategy = paste(name, strategy)) |> 
  mutate(years = step / 365) |>
  filter(years <= 100) |>
  ggplot(aes(x = years, y = value, color = name_strategy, linetype = name_strategy)) + geom_line() + theme_classic() +
  scale_color_manual(values = c("prevA Equal Allocation" = '#d73027', "E_a Equal Allocation" = '#d73027',
                                "prevB Equal Allocation" = '#fc8d59', "E_b Equal Allocation" = '#fc8d59',
                                "prevC Equal Allocation" = '#fdb863', "E_c Equal Allocation" = '#fdb863',
                                "prevA Sequential" = '#041c5a', "E_a Sequential" = '#041c5a',
                                "prevB Sequential" = '#2171b5', "E_b Sequential" = '#2171b5',
                                "prevC Sequential" = '#6baed6', "E_c Sequential" = '#6baed6')) +
  scale_linetype_manual(values =c("prevA Equal Allocation" = 1, "E_a Equal Allocation" = 2,
                                  "prevB Equal Allocation" = 1, "E_b Equal Allocation" = 2,
                                  "prevC Equal Allocation" = 1, "E_c Equal Allocation" = 2,
                                  "prevA Sequential" = 1, "E_a Sequential" = 2,
                                  "prevB Sequential" = 1, "E_b Sequential" = 2,
                                  "prevC Sequential" = 1, "E_c Sequential" = 2)) + xlab("Years") +
  ylab("")





seq_txprob <- prevdf_seq |> pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |> 
  filter(name %in% c('E_a', 'E_b', 'E_c')) |>
  mutate(Years = step / 365) |>
  filter(Years <= 50) |>
  ggplot(aes(x = Years, y = value, color = name)) + geom_line(alpha = 0.65, linetype = 2) + theme_classic() +
  scale_color_manual(name = "Treatment Probability", 
                     values = c("E_a" = '#041c5a', "E_b" = '#2171b5', "E_c" = '#6baed6'),
                     labels = c("E_a" = "Drug A", "E_b" = "Drug B", "E_c" = "Drug C")) +
  ylab("Probability of Treatment")



ea_txprob <- prevdf_ea |> pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |> 
  filter(name %in% c('E_a', 'E_b', 'E_c')) |>
  mutate(Years = step / 365) |>
  filter(Years <= 100) |>
  ggplot(aes(x = Years, y = value, color = name)) + geom_line(alpha = 0.65, linetype = 2) + theme_classic() +
  scale_color_manual(name = "Treatment Probability", 
                     values = c("E_a" = '#d73027', "E_b" = '#fc8d59', "E_c" = '#fdb863'),
                     labels = c("E_a" = "Drug A", "E_b" = "Drug B", "E_c" = "Drug C")) +
  ylab("Probability of Treatment")




seq_resprev <- prevdf_seq |> pivot_longer(cols = c(prevA, prevB, prevC)) |> 
  filter(name %in% c('prevA', 'prevB', 'prevC')) |>
  mutate(Years = step / 365) |>
  filter(Years <= 50) |>
  ggplot(aes(x = Years, y = value, color = name)) + geom_line(alpha = 0.65, linetype = 1) + theme_classic() +
  scale_color_manual(name = "Prevalence of Resistance To:", 
                     values = c("prevA" = '#041c5a', "prevB" = '#2171b5', "prevC" = '#6baed6'),
                     labels = c("prevA" = "Drug A", "prevB" = "Drug B", "prevC" = "Drug C")) +
  geom_hline(yintercept = 0.05, color = "darkgrey", linetype = 3) +
  ylab("Resistance Prevalence")



ea_resprev <- prevdf_ea |> pivot_longer(cols = c(prevA, prevB, prevC)) |> 
  filter(name %in% c('prevA', 'prevB', 'prevC')) |>
  mutate(Years = step / 365) |>
  filter(Years <= 100) |>
  ggplot(aes(x = Years, y = value, color = name)) + geom_line(alpha = 0.65, linetype = 1) + theme_classic() +
  scale_color_manual(name = "Prevalence of Resistance To:", 
                     values = c("prevA" = '#d73027', "prevB" = '#fc8d59', "prevC" = '#fdb863'),
                     labels = c("prevA" = "Drug A", "prevB" = "Drug B", "prevC" = "Drug C")) +
  geom_hline(yintercept = 0.05, color = "darkgrey", linetype = 3) + 
  ylab("Resistance Prevalence")

#save all figures
ggsave("../figures/deterministic_figures/with_legends/plotall_detplots.tiff", plotall_det_prevs_probs, width = 5, height = 3, dpi = 900)


ggsave("../figures/deterministic_figures/with_legends/seq_tx_probs.tiff", seq_txprob, width = 5, height = 2.3, dpi = 900)
ggsave("../figures/deterministic_figures/with_legends/ea_tx_probs.tiff", ea_txprob, width = 5, height = 2.3, dpi = 900)
ggsave("../figures/deterministic_figures/with_legends/seq_resprev.tiff", seq_resprev, width = 5, height = 2.3, dpi = 900)
ggsave("../figures/deterministic_figures/with_legends/ea_resprev.tiff", ea_resprev, width = 5, height = 2.3, dpi = 900)


ggsave("../figures/deterministic_figures/no_legends/seq_resprev_nolab.tiff", seq_resprev + theme(legend.position = "none"), width = 4, height = 2.3, dpi = 900)
ggsave("../figures/deterministic_figures/no_legends/seq_tx_probs_nolab.tiff", seq_txprob + theme(legend.position = "none"), width = 4, height = 2.3, dpi = 900)
ggsave("../figures/deterministic_figures/no_legends/ea_resprev_nolab.tiff", ea_resprev + theme(legend.position = "none"), width = 4, height = 2.3, dpi = 900)
ggsave("../figures/deterministic_figures/no_legends/ea_tx_probs_nolab.tiff", ea_txprob + theme(legend.position = "none"), width = 4, height = 2.3, dpi = 900)


