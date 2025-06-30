#load necessary packages
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(odin)
library(grid)


drugslost_df_ea <- read_csv("../output/2025-06-14ea_sens_fixed_613/ea_sensequal allocationlock_df.csv")
drugslost_df_seq <- read_csv("../output/2025-06-14seq_sens_fixed_613/seq_senssequentiallock_df.csv")
sensitive_lhs_filled <- read_csv("../output/6_11_sensitive_LHS_df_fixed.csv")

avglostdf_toplot <- rbind(drugslost_df_ea, drugslost_df_seq) |> select(-avg_all3_lost) |>
  pivot_wider(names_from = strategy, values_from = avg_drugs_lost) |>
  mutate(diff = sequential - `equal allocation`) |> mutate(sign = ifelse(diff < 0, 1, 0)) 

drugs_lost_scatter_5yrs <- avglostdf_toplot |> filter(step == 1825) |>
  ggplot(aes(y = `equal allocation`, x = `sequential`, color = as.factor(sign))) + geom_point(alpha = 0.5) + theme_classic() +
  ylim(0,3) + xlim(0,3) + scale_color_manual(name = "More drugs lost in:", values = c('1' = '#d73027', '0' = '#08519c'),
                                             labels = c("1" = "Equal Allocation", "0" = "Sequential")) + 
  xlab("Average # Drugs Lost Sequential") + ylab("Average # Drugs Lost Equal Allocation") +
  ggtitle("Time = 5 Years")
  
  drugs_lost_scatter_10yrs <- avglostdf_toplot |> filter(step == 3650) |>
  ggplot(aes(y = `equal allocation`, x = `sequential`, color = as.factor(sign))) + geom_point(alpha = 0.5) + theme_classic() +
  ylim(0,3) + xlim(0,3) + scale_color_manual(name = "More drugs lost in:", values = c('1' = '#d73027', '0' = '#08519c'),
                                             labels = c("1" = "Equal Allocation", "0" = "Sequential")) + 
  xlab("Average # Drugs Lost Sequential") + ylab("Average # Drugs Lost Equal Allocation") +
  ggtitle("Time = 10 Years")
  
  drugs_lost_scatter_20yrs <- avglostdf_toplot |> filter(step == 7300) |>
  ggplot(aes(y = `equal allocation`, x = `sequential`, color = as.factor(sign))) + geom_point(alpha = 0.5) + theme_classic() +
  ylim(0,3) + xlim(0,3) + scale_color_manual(name = "More drugs lost in:", values = c('1' = '#d73027', '0' = '#08519c'),
                                             labels = c("1" = "Equal Allocation", "0" = "Sequential")) + 
  xlab("Average # Drugs Lost Sequential") + ylab("Average # Drugs Lost Equal Allocation") +
  ggtitle("Time = 20 Years")
  
  ggsave("../figures/sens_analysis_figures/drugslost_scatter_5yrs.tiff", drugs_lost_scatter_5yrs, width = 4, height = 3, dpi = 900)
  ggsave("../figures/sens_analysis_figures/drugslost_scatter_20yrs.tiff", drugs_lost_scatter_20yrs, width = 4, height = 3, dpi = 900)
	ggsave("../figures/sens_analysis_figures/drugslost_scatter_10yrs.tiff", drugs_lost_scatter_10yrs, width = 4, height = 3, dpi = 900)
save_edge_cases_all22 <- avglostdf_toplot |> filter(step == 3650) |> arrange(diff) |> mutate(comb = parameter_comb) |> left_join(sensitive_lhs_filled) |>
  filter(diff < 0) 
  
  write_csv(save_edge_cases_all22, "../output/6_23_edgecase_sens_analy_all.csv")
  

   #load in prev dfs
  #load in DFs
eaA_df <- read_csv("../output/2025-06-14ea_sens_fixed_613/ea_sensequal allocationprevAdf.csv")
eaB_df <- read_csv("../output/2025-06-14ea_sens_fixed_613/ea_sensequal allocationprevBdf.csv")
eaC_df <- read_csv("../output/2025-06-14ea_sens_fixed_613/ea_sensequal allocationprevCdf.csv")

seqA_df <- read_csv("../output/2025-06-14seq_sens_fixed_613/seq_senssequentialprevAdf.csv")
seqB_df <- read_csv("../output/2025-06-14seq_sens_fixed_613/seq_senssequentialprevBdf.csv")
seqC_df <- read_csv("../output/2025-06-14seq_sens_fixed_613/seq_senssequentialprevCdf.csv")

seq_fulldf <- left_join(seqA_df, seqB_df) |> left_join(seqC_df)
ea_fulldf <- left_join(eaA_df, eaB_df) |> left_join(eaC_df)

all_combs_toplot <- rbind(seq_fulldf, ea_fulldf) |>
  pivot_wider(names_from = strategy, values_from = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
  mutate(prevAdiffsign = ifelse(prop_runs_prevA_above_sequential - `prop_runs_prevA_above_equal allocation` < 0, 1, 0)) |>
  mutate(prevBdiffsign = ifelse(prop_runs_prevB_above_sequential - `prop_runs_prevB_above_equal allocation` < 0, 1, 0)) |>
  mutate(prevCdiffsign = ifelse(prop_runs_prevC_above_sequential - `prop_runs_prevC_above_equal allocation` < 0, 1, 0))

  drugC_prevcompare_scatter_10yrs_all <- all_combs_toplot |> filter(step == 3650) |>
  ggplot(aes(x = prop_runs_prevC_above_sequential, y = `prop_runs_prevC_above_equal allocation`, color = as.factor(prevCdiffsign))) + geom_point(alpha = 0.5) + theme_classic() + scale_color_manual(name ="More runs hitting threshold in:",  values = c('1' = '#d73027', '0' = '#08519c'), labels = c("1" = "Equal Allocation", "0" = "Sequential")) +
  xlim(0,1) + ylim(0,1) + xlab("Prop hitting threshold drug C, Sequential") + ylab("Prop hitting threshold drug C, Equal Allocation") + ggtitle("Time = 10 Years")

	drugB_prevcompare_scatter_10yrs_all <- all_combs_toplot |> filter(step == 3650) |>
  ggplot(aes(x = prop_runs_prevB_above_sequential, y = `prop_runs_prevB_above_equal allocation`, color = as.factor(prevBdiffsign))) + geom_point(alpha = 0.5) + theme_classic() + scale_color_manual(name ="More runs hitting threshold in:",  values = c('1' = '#d73027', '0' = '#08519c'), labels = c("1" = "Equal Allocation", "0" = "Sequential")) +
  xlim(0,1) + ylim(0,1) + xlab("Prop hitting threshold drug B, Sequential") + ylab("Prop hitting threshold drug B, Equal Allocation") + ggtitle("Time = 10 Years")


  ggsave("../figures/sens_analysis_figures/drugC_prevcompare_scatter_10yrs_all.tiff", drugC_prevcompare_scatter_10yrs_all, width = 5.5, height = 3.5, dpi = 900)
	ggsave("../figures/sens_analysis_figures/drugB_prevcompare_scatter_10yrs_all.tiff", drugB_prevcompare_scatter_10yrs_all, width = 5.5, height = 3.5, dpi = 900)

