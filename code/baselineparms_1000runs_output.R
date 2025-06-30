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


seq_prevdf <- read_csv("../output/baseline_parms_analysis/2025-06-28baselineparms_1000runs_seq/prevdf.csv")
seq_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-28baselineparms_1000runs_seq/lockdf.csv")
ea_prevdf <- read_csv("../output/baseline_parms_analysis/2025-06-28baselineparms_1000runs_ea/prevdf.csv")
ea_lockdf <- read_csv("../output/baseline_parms_analysis/2025-06-28baselineparms_1000runs_ea/lockdf.csv")

seq1000_baselineparms_prop_df <- seq_prevdf |> group_by(run) |>
  mutate(prevA_above_threshold = cumsum(prevA >= 0.05) > 0 ,
         prevB_above_threshold = cumsum(prevB >= 0.05) > 0,
         prevC_above_threshold = cumsum(prevC >= 0.05) > 0) |>
  group_by(step) |>
  summarize(prop_runs_prevA_above = mean(prevA_above_threshold), prop_runs_prevB_above = mean(prevB_above_threshold),
            prop_runs_prevC_above = mean(prevC_above_threshold)) |> mutate(strategy = "sequential") |>
  mutate(Years = step / 365)

write_csv(seq1000_baselineparms_prop_df, "../output/baseline_parms_analysis/seq_baselineparms_propdf.csv")

seq1000_baselineparms_numlost_df <- seq_lockdf |> group_by(run) |>
  mutate(A_out = lock_threshold_A >= 1, B_out = lock_threshold_B >= 1, C_out = lock_threshold_C >= 1) |>
  mutate(num_lost = A_out + B_out + C_out) |>
  group_by(step) |>
  summarize(avg_drugs_lost = mean(num_lost)) |> mutate(strategy = "sequential") |>
  mutate(Years = step / 365)

write_csv(seq1000_baselineparms_numlost_df, "../output/baseline_parms_analysis/seq_baselineparms_avglost.csv")


seq1000_baselineparms_all3lost_df <- seq_lockdf |> group_by(run) |>
  mutate(A_out = lock_threshold_A >= 1, B_out = lock_threshold_B >= 1, C_out = lock_threshold_C >= 1) |>
  mutate(all3lost = (A_out + B_out + C_out == 3)) |>
  group_by(step) |>
  summarize(avg_lost3 = mean(all3lost)) |> mutate(strategy = "sequential") |>
  mutate(Years = step / 365)

write_csv(seq1000_baselineparms_all3lost_df, "../output/baseline_parms_analysis/seq_baselineparms_lostall3.csv")

#equal allocation
ea1000_baselineparms_prop_df <- ea_prevdf |> group_by(run) |>
  mutate(prevA_above_threshold = cumsum(prevA >= 0.05) > 0 ,
         prevB_above_threshold = cumsum(prevB >= 0.05) > 0,
         prevC_above_threshold = cumsum(prevC >= 0.05) > 0) |>
  group_by(step) |>
  summarize(prop_runs_prevA_above = mean(prevA_above_threshold), prop_runs_prevB_above = mean(prevB_above_threshold),
            prop_runs_prevC_above = mean(prevC_above_threshold)) |> mutate(strategy = "equal allocation") |>
	mutate(Years = step / 365)

write_csv(ea1000_baselineparms_prop_df, "../output/baseline_parms_analysis/ea_baselineparms_propdf.csv")

ea1000_baselineparms_numlost_df <- ea_lockdf |> group_by(run) |>
  mutate(A_out = lock_threshold_A >= 1, B_out = lock_threshold_B >= 1, C_out = lock_threshold_C >= 1) |>
  mutate(num_lost = A_out + B_out + C_out) |>
  group_by(step) |>
  summarize(avg_drugs_lost = mean(num_lost)) |> mutate(strategy = "equal allocation") |>
  mutate(Years = step / 365)
write_csv(ea1000_baselineparms_numlost_df, "../output/baseline_parms_analysis/ea_baselineparms_avglost.csv")

ea1000_baselineparms_all3lost_df <- ea_lockdf |> group_by(run) |>
  mutate(A_out = lock_threshold_A >= 1, B_out = lock_threshold_B >= 1, C_out = lock_threshold_C >= 1) |>
  mutate(all3lost = (A_out + B_out + C_out == 3)) |>
  group_by(step) |>
  summarize(avg_lost3 = mean(all3lost)) |> mutate(strategy = "equal allocation") |>
  mutate(Years = step / 365)

write_csv(ea1000_baselineparms_all3lost_df, "../output/baseline_parms_analysis/ea_baselineparms_lostall3.csv")

#make table of number of drugs lost
#calculate the number of simulations that have lost 1, 2, and 3 drugs for each strategy
seq_numslost_count_df <- seq_lockdf |> group_by(run) |>
  mutate(A_out = lock_threshold_A >= 1, B_out = lock_threshold_B >= 1, C_out = lock_threshold_C >= 1) |>
  mutate(num_lost = A_out + B_out + C_out)
seq_numslost <- seq_numslost_count_df |> group_by(step, num_lost) |>
  summarize(count = n(), .groups = "drop") |>
  complete(step, num_lost = 0:3, fill = list(count = 0)) |>
  pivot_wider(names_from = num_lost, values_from = count, values_fill = list(count = 0)) |>
  mutate(Years = step / 365) |> mutate(Strategy = "Sequential") |>
  select(Years, `0`, `1`, `2`, `3`, Strategy)
write_csv(seq_numslost, "../output/baseline_parms_analysis/seq_baselineparms_numslost.csv")

ea_numslost_count_df <- ea_lockdf |> group_by(run) |>
  mutate(A_out = lock_threshold_A >= 1, B_out = lock_threshold_B >= 1, C_out = lock_threshold_C >= 1) |>
  mutate(num_lost = A_out + B_out + C_out)
ea_numslost <- ea_numslost_count_df |>
  group_by(step, num_lost) |>
  summarize(count = n(), .groups = "drop") |>
  complete(step, num_lost = 0:3, fill = list(count = 0)) |>
  pivot_wider(
    names_from = num_lost,
    values_from = count,
    values_fill = list(count = 0)
  ) |>
  mutate(
    Years = step / 365,
    Strategy = "Equal Allocation"
  ) |> select(Years, `0`, `1`, `2`, `3`, Strategy)

write_csv(ea_numslost, "../output/baseline_parms_analysis/ea_baselineparms_numslost.csv")



comb_numslost_df <- rbind(seq_numslost, ea_numslost) |> arrange(Years)

write_csv(comb_numslost_df, "../output/baseline_parms_analysis/combined_baselineparms_numslost.csv")

comb_numslost_df_filtered <- comb_numslost_df |> filter( Years %in% c(2,5,10,30))

write_csv(comb_numslost_df_filtered, "../output/baseline_parms_analysis/combined_baselineparms_numslost_filtered.csv")


