#write a run script for easily running any model:
model_run_function <- function(model, parms, dt){
  mod <- do.call(model$new, parms)
  return(as.data.frame(mod$run(dt)))
}

#write a function that runs both deterministic models and plots the output
run_both_det_mods <- function(eamod_det, seqmod_det, baselineparms, dt){
  eadf <- model_run_function(model = eamod_det, parms = baseline_parms_5_26, dt = dt) |> mutate("strategy" = "equal allocation")
  seqdf <- model_run_function(model = seqmod_det, parms = baseline_parms_5_26, dt = dt) |> mutate("strategy" = "sequential")
  comb_df <- rbind(eadf, seqdf)
  
  #plot treatment probabilities and prevalences
  txprob_prevplot <- comb_df |> select(step, E_a, E_b, E_c, prevA, prevB, prevC, strategy) |>
    filter(step < 22000) |>
    pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic() + geom_hline(yintercept = 0.05, color = "darkgrey", linetype = 2)
  
  retxprobs <- comb_df |> select(step, L_a_b, L_a_c, L_b_a, L_b_c, L_c_a, L_c_b, strategy) |>
    filter(step < 22000) |>
    pivot_longer(cols = c(L_a_b, L_a_c, L_b_a, L_b_c, L_c_a, L_c_b)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic()
  
  lrt_probs <- comb_df |> select(step, v_abc, v_bc, v_ac, v_c, strategy) |>
    filter(step < 22000) |>
    pivot_longer(cols = c(v_abc, v_bc, v_ac, v_c)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic()
  
  totals_plot <- comb_df |> select(step, sum_N, sum_S, strategy) |>
    filter(step < 22000) |>
    pivot_longer(cols = c(sum_N, sum_S)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic()
  return(list(comb_df, txprob_prevplot, retxprobs, lrt_probs, totals_plot))
  
}

#make a function to show the relevant plots when I've run a specific set of parameters in a specific model
run_mod_and_plot <- function(model, parms, dt){
  runmod_df <- model_run_function(model, parms, dt)
  #plot prevalences
  prevplot <- runmod_df |> select(step, E_a, E_b, E_c, prevA, prevB, prevC) |> pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
    ggplot(aes(x = step, y = value, color = name)) + geom_line(alpha = 0.35) + theme_classic() + geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey")
  
  prevplot_zoom <- runmod_df |> select(step, E_a, E_b, E_c, prevA, prevB, prevC) |> filter(step < 20000) |>
    pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
    ggplot(aes(x = step, y = value, color = name)) + geom_line(alpha = 0.35) + theme_classic() + geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey")
  
  prevplot_sep <- runmod_df |> select(step, E_a, E_b, E_c, prevA, prevB, prevC) |>
    pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
    ggplot(aes(x = step, y = value, color = name)) + geom_line(alpha = 0.35) + theme_classic() + geom_hline(yintercept = 0.05, linetype = 2, color = "darkgrey") + facet_wrap(~name)
  #check mass balance
  Nplot <- runmod_df |> select(step, sum_S, sum_N) |>
    pivot_longer(cols = c(sum_S, sum_N)) |>
    ggplot(aes(x = step, y = value, color = name)) + geom_line(alpha = 0.5) + theme_classic()
  #plot all YZ compartments
  YZ_plots <- runmod_df |> select(step, sum_YZ_s, sum_YZ_a, sum_YZ_b, sum_YZ_c, sum_YZ_ab, sum_YZ_ac, sum_YZ_bc, sum_YZ_abc) |>
    pivot_longer(cols = c(sum_YZ_s, sum_YZ_a, sum_YZ_b, sum_YZ_c, sum_YZ_ab, sum_YZ_ac, sum_YZ_bc, sum_YZ_abc)) |>
    ggplot(aes(x = step, y = value, color = name)) + geom_line(alpha = 0.5) + theme_classic() + facet_wrap(~name)
  #overall prevalence of infection in system
  prev <- runmod_df |> select(step, prev) |>
    ggplot(aes(x = step, y = prev)) + geom_line() + theme_classic()
  #make sure tx probs sum to 1
  txprobssum <- runmod_df |> select(step, sum_tx_probs) |>
    ggplot(aes(x = step, y = sum_tx_probs)) + geom_line() + theme_classic()
  
  return(list(runmod_df, prevplot, prevplot_zoom, prevplot_sep, Nplot, YZ_plots,prev, txprobssum))
  
}

#write a function that runs both deterministic models and plots the output, updated
run_both_det_mods_5_29 <- function(eamod_det, seqmod_det, parms, dt){
  eadf <- model_run_function(model = eamod_det, parms = parms, dt = dt) |> mutate("strategy" = "equal allocation") |> select(-A_remain_lock, -B_remain_lock, -C_remain_lock)
  seqdf <- model_run_function(model = seqmod_det, parms = parms, dt = dt) |> mutate("strategy" = "sequential")
  comb_df <- rbind(eadf, seqdf)
  
  #plot treatment probabilities
  txprob_prevplot <- comb_df |> select(step, E_a, E_b, E_c, prevA, prevB, prevC, strategy) |>
    pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic() + geom_hline(yintercept = 0.05, color = "darkgrey", linetype = 2)
  
  txprob_prevplot_zoom <- comb_df |> select(step, E_a, E_b, E_c, prevA, prevB, prevC, strategy) |>
    filter(step < 22000) |>
    pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic() + geom_hline(yintercept = 0.05, color = "darkgrey", linetype = 2)
  
  #plot totals
  totals_plot <- comb_df |> select(step, sum_N, sum_S, strategy) |>
    pivot_longer(cols = c(sum_N, sum_S)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic()
  
  #plot compartments
  compartments_plot <- comb_df |> select(step, sum_YZ_a, sum_YZ_b, sum_YZ_c, sum_YZ_ac, sum_YZ_bc, sum_YZ_abc, strategy) |>
    pivot_longer(cols = c(sum_YZ_a, sum_YZ_b, sum_YZ_c, sum_YZ_ac, sum_YZ_bc, sum_YZ_abc)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic() + facet_wrap(~name)
  
  
  
  
  return(list(comb_df, txprob_prevplot, txprob_prevplot_zoom, totals_plot, compartments_plot))
  
}


#write a function quickly that will do 15 simulations of 1 stochastic model and save all output
run_stochastic_saveoutput_6_2 <- function(model, parms, runs, dt){
  mod <- do.call(model$new, parms)
  df <- data.frame()
  for(r in c(1:runs)){
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r)
    df <- rbind(df, tmpdf)
  }
  return(df)
}

#write a function that runs the deterministic versions and the stochastic versions and outputs various plots 

#function for running stochastic function but not keeping all output; just keep the compartments and the prevalences / prob of tx
run_stochastic_savebrief_6_3 <- function(model, parms, runs, dt){
  mod <- do.call(model$new, parms)
  df1 <- data.frame() #for keeping tx probs and prevs
  df2 <- data.frame() #for keeping compartments
  if(runs == "det"){
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = "det")
    df1 <- tmpdf |> select(step,  E_a, E_b, E_c, prevA, prevB, prevC, run)
    df2 <- tmpdf |> select(step, sum_N, sum_S, sum_Y_s, sum_Y_a, sum_Y_b, sum_Y_c,
                           sum_Y_ab, sum_Y_ac, sum_Y_bc, sum_Y_abc,sum_Z_s, sum_Z_a,
                           sum_Z_b, sum_Z_c, sum_Z_ab,sum_Z_ac, sum_Z_bc, sum_Z_abc, run)
  }
  else{
    for(r in c(1:runs)){
      tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r) |> select(step, E_a, E_b, E_c, prevA, prevB, prevC,
                                                                       sum_N, sum_S, sum_Y_s, sum_Y_a, sum_Y_b, sum_Y_c,
                                                                       sum_Y_ab, sum_Y_ac, sum_Y_bc, sum_Y_abc,
                                                                       sum_Z_s, sum_Z_a, sum_Z_b, sum_Z_c, sum_Z_ab,
                                                                       sum_Z_ac, sum_Z_bc, sum_Z_abc, run)
      tmpdf1 <- tmpdf |> select(step,  E_a, E_b, E_c, prevA, prevB, prevC, run)
      tmpdf2 <- tmpdf |> select(step, sum_N, sum_S, sum_Y_s, sum_Y_a, sum_Y_b, sum_Y_c,
                                sum_Y_ab, sum_Y_ac, sum_Y_bc, sum_Y_abc,sum_Z_s, sum_Z_a,
                                sum_Z_b, sum_Z_c, sum_Z_ab,sum_Z_ac, sum_Z_bc, sum_Z_abc, run)
      df1 <- rbind(df1, tmpdf1)
      df2 <- rbind(df2, tmpdf2)}
  }
  return(list(df1, df2))
}

#function to calculate TOLs
find_threshold_time_prevA <- function(data, threshold) {
  idx <- which(data$prevA >= threshold)
  if (length(idx) > 0) {
    return(data$step[min(idx)])
  } else {
    return(NA)  # Return NA if the threshold is never reached
  }
}

find_threshold_time_prevB <- function(data, threshold) {
  idx <- which(data$prevB >= threshold)
  if (length(idx) > 0) {
    return(data$step[min(idx)])
  } else {
    return(NA)  # Return NA if the threshold is never reached
  }
}

find_threshold_time_prevC <- function(data, threshold) {
  idx <- which(data$prevC >= threshold)
  if (length(idx) > 0) {
    return(data$step[min(idx)])
  } else {
    return(NA)  # Return NA if the threshold is never reached
  }
}

run_allmods_together_6_3 <- function(EA_det_mod, seq_det_mod, EA_stoch_mod, seq_stoch_mod, parms, runs, dt){
  #run EA deterministic
  EA_det_dfs <- run_stochastic_savebrief_6_3(model = EA_det_mod, parms = parms, runs = "det", dt)
  
  #run seq deterministic
  seq_det_dfs <- run_stochastic_savebrief_6_3(model = seq_det_mod, parms = parms, runs = "det", dt) 
  
  #run stochastic models
  EA_stoch_dfs <- run_stochastic_savebrief_6_3(model = EA_stoch_mod, parms = parms, runs = runs, dt = dt)
  seq_stoch_dfs <- run_stochastic_savebrief_6_3(model = seq_stoch_mod, parms = parms, runs = runs, dt = dt)
  
  #unlist and join 
  EA_det_prevs_df <- EA_det_dfs[[1]] |> mutate(strategy = "equal allocation", model = "deterministic")
  EA_det_comps_df <- EA_det_dfs[[2]] |> mutate(strategy = "equal allocation", model = "deterministic")
  seq_det_prevs_df <- seq_det_dfs[[1]] |> mutate(strategy = "seq", model = "deterministic")
  seq_det_comps_df <- seq_det_dfs[[2]] |> mutate(strategy = "seq", model = "deterministic") 
  EA_stoch_prevs_df <- EA_stoch_dfs[[1]] |> mutate(strategy = "equal allocation", model = "stochastic")
  EA_stoch_comps_df <- EA_stoch_dfs[[2]] |> mutate(strategy = "equal allocation", model = "stochastic")
  seq_stoch_prevs_df <- seq_stoch_dfs[[1]] |> mutate(strategy = "seq", model = "stochastic")
  seq_stoch_comps_df <- seq_stoch_dfs[[2]] |> mutate(strategy = "seq", model = "stochastic")
  
  #combine prevalence dfs
  comb_prev_df <- rbind(EA_det_prevs_df, seq_det_prevs_df, EA_stoch_prevs_df, seq_stoch_prevs_df)
  comb_comps_df <- rbind(EA_det_comps_df, seq_det_comps_df, EA_stoch_comps_df, seq_stoch_comps_df)
  
  #make TOL dfs
  
  return(list(comb_prev_df, comb_comps_df, parms))
}

#write a function that plots the compartments 
plot_strat_comp_traj <- function(df2, strat, compartment){
  #make plot of all trajectories
  alltrajplt <- df2 |> select(step, compartment, strategy, model, run) |>
    filter(strategy == strat) |> pivot_longer(cols = c(compartment)) |> 
    ggplot(aes(x = step, y = value, linetype = model, color = model, group = interaction(run, model), alpha = model)) + geom_line() + theme_classic() +
    scale_linetype_manual(values = c("deterministic"= 2, "stochastic" = 1)) +
    scale_alpha_manual(values = c("deterministic"= 1, "stochastic" = 0.25)) +
    scale_color_manual(values = c("deterministic" = "black", "stochastic" = "lightblue")) +
    ggtitle(compartment)
  
  det_sum_df <- df2 |> select(step, compartment, strategy, model, run) |> 
    filter(strategy == strat, model == "deterministic") |>
    pivot_longer(cols = c(compartment))
  
  avgplot <- df2 |> select(step, compartment, strategy, model, run) |>
    filter(strategy == strat, model == "stochastic") |>
    group_by(step) |>
    pivot_longer(cols = c(compartment)) |>
    summarize(avg = mean(value), se = sd(value)/sqrt(n())) |>
    mutate(model = "stochastic") |>
    ggplot(aes(x = step, y = avg, color = model), fill = "lightblue") + geom_line(fill = "lightblue") + geom_ribbon(aes(ymin = avg - 1.96*se, ymax = avg + 1.96*se), color = NA, fill = "lightblue", alpha = 0.5) + theme_classic() + 
    geom_line(aes(x = step, y = value, color = model), data = det_sum_df) +
    scale_color_manual(values = c("deterministic" = "black", "stochastic" = "lightblue")) +
    ggtitle(compartment)
  return(list(alltrajplt, avgplot))
}

plot_trajectories_compartments <- function(df2, foldername){
  complist <- colnames(df2)[2:19]
  #iterate through compartments and run EA and sequential 
  for(comp in complist){
    tmp1 <- plot_strat_comp_traj(df2 = df2, strat = "equal allocation",compartment = comp)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "alltrajs/" , "EA", comp, "all_traj_plot.tiff"), tmp1[[1]], width = 7, height = 4)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "avgtrajs/", "EA", comp, "avg_traj_plot.tiff"), tmp1[[2]], width = 7, height = 4)
    tmp2 <- plot_strat_comp_traj(df2 = df2, strat = "seq",compartment = comp)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "alltrajs/" , "seq", comp,  "all_traj_plot.tiff"), tmp2[[1]], width = 7, height = 4)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "avgtrajs/", "seq", comp, "avg_traj_plot.tiff"), tmp2[[2]], width = 7, height = 4)
  }
  
}

#plotting each trajectory is vector memory intensive, just plot average trajectories
plot_strat_avg_traj <- function(df2, strat, compartment){
  #make plot of all trajectories
  det_sum_df <- df2 |> select(step, compartment, strategy, model, run) |> 
    filter(strategy == strat, model == "deterministic") |>
    pivot_longer(cols = c(compartment))
  
  avgplot <- df2 |> select(step, compartment, strategy, model, run) |>
    filter(strategy == strat, model == "stochastic") |>
    group_by(step) |>
    pivot_longer(cols = c(compartment)) |>
    summarize(avg = mean(value), se = sd(value)/sqrt(n())) |>
    mutate(model = "stochastic") |>
    ggplot(aes(x = step, y = avg, color = model), fill = "lightblue") + geom_line(fill = "lightblue") + geom_ribbon(aes(ymin = avg - 1.96*se, ymax = avg + 1.96*se), color = NA, fill = "lightblue", alpha = 0.5) + theme_classic() + 
    geom_line(aes(x = step, y = value, color = model), data = det_sum_df) +
    scale_color_manual(values = c("deterministic" = "black", "stochastic" = "lightblue")) +
    ggtitle(compartment)
  return(avgplot)
}

plot_strat_avg_traj_naRM <- function(df2, strat, compartment){
  #make plot of all trajectories
  det_sum_df <- df2 |> select(step, compartment, strategy, model, run) |> 
    filter(strategy == strat, model == "deterministic") |>
    pivot_longer(cols = c(compartment))
  
  avgplot <- df2 |> select(step, compartment, strategy, model, run) |>
    filter(strategy == strat, model == "stochastic") |>
    group_by(step) |>
    pivot_longer(cols = c(compartment)) |>
    summarize(avg = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE)/sqrt(n())) |>
    mutate(model = "stochastic") |>
    ggplot(aes(x = step, y = avg, color = model), fill = "lightblue") + geom_line(fill = "lightblue") + geom_ribbon(aes(ymin = avg - 1.96*se, ymax = avg + 1.96*se), color = NA, fill = "lightblue", alpha = 0.5) + theme_classic() + 
    geom_line(aes(x = step, y = value, color = model), data = det_sum_df) +
    scale_color_manual(values = c("deterministic" = "black", "stochastic" = "lightblue")) +
    ggtitle(compartment)
  return(avgplot)
}



plot_avg_trajectories_compartments <- function(df2, foldername){
  complist <- colnames(df2)[2:19]
  #iterate through compartments and run EA and sequential 
  for(comp in complist){
    tmp1 <- plot_strat_avg_traj(df2 = df2, strat = "equal allocation",compartment = comp)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "avgtrajs/", "EA", comp, "avg_traj_plot.tiff"), tmp1, width = 7, height = 4)
    tmp2 <- plot_strat_avg_traj(df2 = df2, strat = "seq",compartment = comp)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "avgtrajs/", "seq", comp, "avg_traj_plot.tiff"), tmp2, width = 7, height = 4)
  }
  
}

plot_avg_trajectories_compartments_naRM <- function(df2, foldername){
  complist <- colnames(df2)[2:19]
  #iterate through compartments and run EA and sequential 
  for(comp in complist){
    tmp1 <- plot_strat_avg_traj_naRM(df2 = df2, strat = "equal allocation",compartment = comp)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "avgtrajs/", "EA", comp, "avg_traj_plot.tiff"), tmp1, width = 7, height = 4)
    tmp2 <- plot_strat_avg_traj(df2 = df2, strat = "seq",compartment = comp)
    ggsave(paste0("../figures/", Sys.Date(), foldername, "/", "avgtrajs/", "seq", comp, "avg_traj_plot.tiff"), tmp2, width = 7, height = 4)
  }
  
}



make_prop_plot <- function(df1, threshold){
  prop_df <- df1 |> group_by(run, strategy) |>
  filter(model == "stochastic") |>
  mutate(prevA_above_threshold = cumsum(prevA >= threshold) > 0 ,
         prevB_above_threshold = cumsum(prevB >= threshold) > 0,
         prevC_above_threshold = cumsum(prevC >= threshold) > 0) |>
  group_by(step,strategy) |>
  summarize(prop_runs_prevA_above = mean(prevA_above_threshold), prop_runs_prevB_above = mean(prevB_above_threshold), 
            prop_runs_prevC_above = mean(prevC_above_threshold)) |>
  pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above))
  
  prop_plot <- prop_df |> mutate(full_name = paste(strategy, name), model = "Stochastic") |>
    ggplot(aes(x = step, y = value, color = full_name, linetype = full_name)) + geom_line() +
    theme_classic() +
    scale_color_manual(values = c("equal allocation prop_runs_prevA_above" = '#d73027',
                                  "seq prop_runs_prevA_above" = '#041c5a',
                                  "equal allocation prop_runs_prevB_above" = '#fc8d59',
                                  "seq prop_runs_prevB_above" = '#2171b5',
                                  "equal allocation prop_runs_prevC_above" = '#fdb863',
                                  "seq prop_runs_prevC_above" = '#6baed6'),
                       labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A", 
                                  "Sequential, Drug B", "Sequential, Drug C"),
                       name = "Strategy and Drug Resistance Prevalence") +
    scale_linetype_manual(values = c("equal allocation prop_runs_prevA_above" = 1,
                                     "seq prop_runs_prevA_above" = 2,
                                     "equal allocation prop_runs_prevB_above" = 1,
                                     "seq prop_runs_prevB_above" = 2,
                                     "equal allocation prop_runs_prevC_above" = 1,
                                     "seq prop_runs_prevC_above" = 2),
                          labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A", 
                                     "Sequential, Drug B", "Sequential, Drug C"),
                          name = "Strategy and Drug Resistance Prevalence") +
    ylab("Proportion of simulations") +
    xlab("Days")
  return(prop_plot)
}


plot_deterministic_prevplot <- function(df1){
  df1 |> filter(model == "deterministic") |>
  pivot_longer(cols = c(E_a, E_b, E_c, prevA, prevB, prevC)) |>
  ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic()}

plot_det_compsplot <- function(df2){
  df2 |> filter(model == "deterministic") |>
    mutate(sum_YZ_s = sum_Y_s + sum_Z_s,
           sum_YZ_a = sum_Y_a + sum_Z_a,
           sum_YZ_b = sum_Y_b + sum_Z_b,
           sum_YZ_c = sum_Y_c + sum_Z_c,
           sum_YZ_ab = sum_Y_ab + sum_Z_ab,
           sum_YZ_ac = sum_Y_ac + sum_Z_ac,
           sum_YZ_bc = sum_Y_bc + sum_Z_bc,
           sum_YZ_abc = sum_Y_abc + sum_Z_abc) |>
    select(step, sum_N, sum_S, sum_YZ_s, sum_YZ_a, sum_YZ_b, sum_YZ_c, sum_YZ_ab, sum_YZ_ac, sum_YZ_bc, sum_YZ_abc, run, strategy, model) |>
    pivot_longer(cols = c(sum_N, sum_S, sum_YZ_s, sum_YZ_a, sum_YZ_b, sum_YZ_c, sum_YZ_ab, sum_YZ_ac, sum_YZ_bc, sum_YZ_abc)) |>
    ggplot(aes(x = step, y = value, color = name, linetype = strategy)) + geom_line(alpha = 0.5) + theme_classic()
}

#for running sensitivity analysis 
run_stochastic_brief_6_7 <- function(model, parms, runs, dt){
  mod <- do.call(model$new, parms)
  df <- data.frame()
  for(r in c(1:runs)){
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r) |> select(step, run, prevA, prevB, prevC)
    df <- rbind(df, tmpdf)
  }
  return(df)
}

generate_parms_unsens <- function(dfrow, threshold){
  b <- unlist(dfrow[3]) #transmission probability per partnership
  w_a <- 10^-4 #rate of resistance to A emerging on treatment with A #this is a probability, over course of treatment
  w_b <- 10^-5 #rate of resistance to B emerging on treatment with B
  w_c <- 10^-6 #rate of resistance to C emerging on treatment with C
  T_s <- 1/(unlist(dfrow[7])) #time to treatment for symptomatic infection
  T_m <- unlist(dfrow[8])/365  #screening rate (time to treatment for asymptomatic infection)
  T_sr <- T_s/unlist(dfrow[9]) #time to retreatment for symptomatic infection, if failure ##rethink this?
  kappa <- unlist(dfrow[10]) #probability of retreatment after treatment failure, symptomatic
  f_a <- 0.98 #relative fitness when resistant to A vs pan-susceptible 
  f_b <- 0.95 #relative fitness when resistant to B vs pan-susceptible
  f_c <- 0.95 #relative fitness when resistant to C vs pan-susceptible
  f_ab <- f_a*f_b #relative fitness when resistant to A&B vs pan-susceptible
  f_ac <- f_a*f_c #relative fitness when resistant to A&C vs pan-susceptible
  f_bc <- f_b*f_c #relative fitness when resistant to B&C vs pan-susceptible
  f_abc <- f_a*f_b*f_c #relative fitness when resistant to A,B&C vs pan-susceptible
  sigma <- unlist(dfrow[2])  #probability of symptomatic infection
  N_risk <- 3 #number of risk strata
  d <- 1/(unlist(dfrow[6])) #natural recovery rate from infection #this is probably shorter, Barbee study
  pop <- 10^6 #pop size
  pop.p <- c(0.3, 0.6, 0.1) #relative size of each risk group; low, M, high
  c_min <- unlist(dfrow[5]) #rate of partner change in lowest risk group
  activities <- c(1*c_min/365,
                  5*c_min/365,
                  20*c_min/365)  #sexual contacts per day
  epsilon <- 0.2228739 #mixing parameter (for non-assortativity)
  
  #risk stratified mixing 
  beta <- ((1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + epsilon*activities/(pop*pop.p)*diag(3))*b
  
  N_ini <- c(pop*pop.p[1],
             pop*pop.p[2],
             pop*pop.p[3])
  #distribute GC cases to have overall 3% prevalence
  x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)
  gc_lo <- round(N_ini[1]*x*0.029,0)
  gc_md <- round(N_ini[2]*x*0.154,0)
  gc_hi <- round(N_ini[3]*x*0.817,0)
  
  prev_symp <- 0.10
  Y_gc_lo <- gc_lo*prev_symp
  Z_gc_lo <- gc_lo*(1-prev_symp)
  Y_gc_md <- gc_md*prev_symp
  Z_gc_md <- gc_md*(1-prev_symp)
  Y_gc_hi <- gc_hi*prev_symp
  Z_gc_hi <- gc_hi*(1-prev_symp)
  
  Y_s_ini <- round(c(Y_gc_lo, Y_gc_md, Y_gc_hi))
  Z_s_ini <- round(c(Z_gc_lo, Z_gc_md, Z_gc_hi)) - c(0,1,1)
  
  Y_a_ini <- c(0,0,0)
  Z_a_ini <- c(0,1,1)
  
  Y_b_ini <- c(0,0,0)
  Z_b_ini <- c(0,0,0)
  
  Y_c_ini <- c(0,0,0)
  Z_c_ini <- c(0,0,0)
  
  Y_ab_ini <- c(0,0,0)
  Z_ab_ini <- c(0,0,0)
  
  Y_ac_ini <- c(0,0,0)
  Z_ac_ini <- c(0,0,0)
  
  Y_bc_ini <- c(0,0,0)
  Z_bc_ini <- c(0,0,0)
  
  Y_abc_ini <- c(0,0,0)
  Z_abc_ini <- c(0,0,0)
  
  S_ini <- c(N_ini[1] - gc_lo,
             N_ini[2] - gc_md,
             N_ini[3] - gc_hi)
  
  return(list(S_ini = S_ini,Y_s_ini = Y_s_ini, Z_s_ini = Z_s_ini, Y_a_ini = Y_a_ini, Z_a_ini = Z_a_ini, Y_b_ini = Y_b_ini,
              Z_b_ini = Z_b_ini, Y_c_ini = Y_c_ini, Z_c_ini = Z_c_ini, Y_ab_ini = Y_ab_ini, Z_ab_ini = Z_ab_ini, Y_ac_ini = Y_ac_ini,
              Z_ac_ini = Z_ac_ini, Y_bc_ini = Y_bc_ini,Z_bc_ini = Z_bc_ini,Y_abc_ini = Y_abc_ini,Z_abc_ini = Z_abc_ini, beta = beta,
              w_a = w_a, w_b = w_b, w_c = w_c, T_s = T_s, T_m = T_m, T_sr = T_sr, kappa = kappa, f_a = f_a, f_b = f_b, f_c = f_c, 
              f_ab = f_ab, f_ac = f_ac, f_bc = f_bc, f_abc = f_abc, sigma = sigma, N_risk = N_risk, d = d, threshold = threshold))
  
}


make_prop_df <- function(outputdf, threshold, times){
  prop_df <- outputdf |> group_by(run) |>
    mutate(prevA_above_threshold = cumsum(prevA >= threshold) > 0 ,
           prevB_above_threshold = cumsum(prevB >= threshold) > 0,
           prevC_above_threshold = cumsum(prevC >= threshold) > 0) |>
    group_by(step) |>
    summarize(prop_runs_prevA_above = mean(prevA_above_threshold), prop_runs_prevB_above = mean(prevB_above_threshold), 
              prop_runs_prevC_above = mean(prevC_above_threshold)) |>
    filter(step %in% times)
  return(prop_df)
}


run_unsens_analysis <- function(model, strategy, threshold, paramdf, runs, years, times, newfoldername, prefix){
  #make new folder to save things in:
  if (!dir.exists(paste0( "../output/", Sys.Date(), newfoldername))){
    dir.create(paste0("../output/", Sys.Date(), newfoldername), recursive = TRUE)
  }
  
  #initialize dataframes
  prevA_threshold_met_df <- data.frame()
  prevB_threshold_met_df <- data.frame()
  prevC_threshold_met_df <- data.frame()
  
  #iterate through the rows of the df
  for(row in c(1:nrow(paramdf))){
    dfrow = paramdf[row,]
    int_df <- make_prop_df(run_stochastic_brief_6_7(model = model, 
                                                    parms = generate_parms_unsens(dfrow, threshold = threshold), 
                                                    runs = runs, dt = seq(0,365*years, 1)), threshold= threshold, times = times)
    #add a column saying which parameter combination this is
    int_df_lab <- int_df |> mutate("parameter_comb" = row, "strategy" = strategy) |>
      select(step, parameter_comb, strategy, prop_runs_prevA_above,prop_runs_prevB_above, prop_runs_prevC_above)
    #save intermediate dataframe
    write_csv(int_df_lab, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "int_df_", row, ".csv"))
    
    #update running dataframes
    prevA_threshold_met_df <- rbind(prevA_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevA_above))
    prevB_threshold_met_df <- rbind(prevB_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevB_above))
    prevC_threshold_met_df <- rbind(prevC_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevC_above))
    
    write_csv(prevA_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevAdf.csv"))
    write_csv(prevB_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevBdf.csv"))
    write_csv(prevC_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevCdf.csv"))
  }
  
  write_csv(prevA_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevAdf.csv"))
  write_csv(prevB_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevBdf.csv"))
  write_csv(prevC_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevCdf.csv"))
  
  return(list(prevA_threshold_met_df, prevB_threshold_met_df, prevC_threshold_met_df))
  
  
}


#for sensitivity analysis
generate_parms_sens <- function(dfrow, threshold){
  set.seed(06092025)
  #Define parameters -- from 6/6 calibration
  b <- 0.5156510 #transmission probability per partnership
  w_a <- unlist(dfrow[5])#rate of resistance to A emerging on treatment with A #this is a probability, over course of treatment
  w_b <- unlist(dfrow[6]) #rate of resistance to B emerging on treatment with B
  w_c <- unlist(dfrow[7]) #rate of resistance to C emerging on treatment with C
  T_s <- 1/(13.8155932) #time to treatment for symptomatic infection
  T_m <- 0.3973783/365  #screening rate (time to treatment for asymptomatic infection)
  T_sr <- T_s/3 #time to retreatment for symptomatic infection, if failure ##rethink this?
  kappa <- 0.90 #probability of retreatment after treatment failure, symptomatic
  f_a <- unlist(dfrow[8]) #relative fitness when resistant to A vs pan-susceptible 
  f_b <- unlist(dfrow[9]) #relative fitness when resistant to B vs pan-susceptible
  f_c <- unlist(dfrow[10]) #relative fitness when resistant to C vs pan-susceptible
  f_ab <- f_a*f_b #relative fitness when resistant to A&B vs pan-susceptible
  f_ac <- f_a*f_c #relative fitness when resistant to A&C vs pan-susceptible
  f_bc <- f_b*f_c #relative fitness when resistant to B&C vs pan-susceptible
  f_abc <- f_a*f_b*f_c #relative fitness when resistant to A,B&C vs pan-susceptible
  sigma <- 0.5098301  #probability of symptomatic infection
  N_risk <- 3 #number of risk strata
  d <- 1/(95.1613537) #natural recovery rate from infection #this is probably shorter, Barbee study
  pop <- 10^6 #pop size
  pop.p <- c(0.3, 0.6, 0.1) #relative size of each risk group; low, M, high
  c_min <- 1.3905583 #rate of partner change in lowest risk group
  activities <- c(1*c_min/365, 
                  5*c_min/365, 
                  20*c_min/365)  #sexual contacts per day
  epsilon <- 0.2228739 #mixing parameter (for non-assortativity)
  
  #risk stratified mixing 
  beta <- ((1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + epsilon*activities/(pop*pop.p)*diag(3))*b
  
  N_ini <- c(pop*pop.p[1],
             pop*pop.p[2],
             pop*pop.p[3])
  
  #distribute GC cases to have overall 3% prevalence
  x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)
  gc_lo <- round(N_ini[1]*x*0.029,0)
  gc_md <- round(N_ini[2]*x*0.154,0)
  gc_hi <- round(N_ini[3]*x*0.817,0)
  
  prev_symp <- 0.10
  Y_gc_lo <- gc_lo*prev_symp
  Z_gc_lo <- gc_lo*(1-prev_symp)
  Y_gc_md <- gc_md*prev_symp
  Z_gc_md <- gc_md*(1-prev_symp)
  Y_gc_hi <- gc_hi*prev_symp
  Z_gc_hi <- gc_hi*(1-prev_symp)
  
  #determine the starting resistances first
  all_starting_a <- c(rmultinom(n = 1, size = unlist(dfrow[2]), c(0.0029, 0.0154, 0.0817, 0.0261,0.1386,0.7353))) #dividing up probability of symptom status and risk group
  Y_a_ini <- all_starting_a[1:3]
  Z_a_ini <- all_starting_a[4:6]
  
  all_starting_b <- c(rmultinom(n = 1, size = unlist(dfrow[3]), c(0.0029, 0.0154, 0.0817, 0.0261,0.1386,0.7353))) #dividing up probability of symptom status and risk group
  Y_b_ini <- all_starting_b[1:3]
  Z_b_ini <- all_starting_b[4:6]
  
  all_starting_c <- c(rmultinom(n = 1, size = unlist(dfrow[4]), c(0.0029, 0.0154, 0.0817, 0.0261,0.1386,0.7353))) #dividing up probability of symptom status and risk group
  Y_c_ini <- all_starting_c[1:3]
  Z_c_ini <- all_starting_c[4:6]
  
  
  Y_s_ini <- round(c(Y_gc_lo, Y_gc_md, Y_gc_hi)) - Y_a_ini - Y_b_ini - Y_c_ini
  Z_s_ini <- round(c(Z_gc_lo, Z_gc_md, Z_gc_hi)) - Z_a_ini - Z_b_ini - Z_c_ini
  
  Y_ab_ini <- c(0,0,0)
  Z_ab_ini <- c(0,0,0)
  
  Y_ac_ini <- c(0,0,0)
  Z_ac_ini <- c(0,0,0)
  
  Y_bc_ini <- c(0,0,0)
  Z_bc_ini <- c(0,0,0)
  
  Y_abc_ini <- c(0,0,0)
  Z_abc_ini <- c(0,0,0)
  
  S_ini <- c(N_ini[1] - gc_lo,
             N_ini[2] - gc_md,
             N_ini[3] - gc_hi)
  
  N = sum(Y_a_ini) + sum(Z_a_ini) + sum(Y_b_ini) + sum(Z_b_ini) + sum(Y_c_ini) + sum(Z_c_ini) + sum(Y_s_ini) + sum(Z_s_ini) + sum(S_ini)
  
  return(list(S_ini = S_ini,Y_s_ini = Y_s_ini, Z_s_ini = Z_s_ini, Y_a_ini = Y_a_ini, Z_a_ini = Z_a_ini, Y_b_ini = Y_b_ini,
              Z_b_ini = Z_b_ini, Y_c_ini = Y_c_ini, Z_c_ini = Z_c_ini, Y_ab_ini = Y_ab_ini, Z_ab_ini = Z_ab_ini, Y_ac_ini = Y_ac_ini,
              Z_ac_ini = Z_ac_ini, Y_bc_ini = Y_bc_ini,Z_bc_ini = Z_bc_ini,Y_abc_ini = Y_abc_ini,Z_abc_ini = Z_abc_ini, beta = beta,
              w_a = w_a, w_b = w_b, w_c = w_c, T_s = T_s, T_m = T_m, T_sr = T_sr, kappa = kappa, f_a = f_a, f_b = f_b, f_c = f_c, 
              f_ab = f_ab, f_ac = f_ac, f_bc = f_bc, f_abc = f_abc, sigma = sigma, N_risk = N_risk, d = d, threshold = threshold))
  
}

#update the run function to use random seeds
run_stochastic_brief_6_9 <- function(model, parms, runs, dt){
  seeds <- sample(1:10^6, size = runs, replace = FALSE)
  mod <- do.call(model$new, parms)
  df <- data.frame()
  for(r in c(1:runs)){
    set.seed(seeds[r])
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r) |> select(step, run, prevA, prevB, prevC)
    df <- rbind(df, tmpdf)
  }
  return(list(df, seeds))
}

run_sens_analysis <- function(model, strategy, threshold, paramdf, runs, years, times, newfoldername, prefix){
  #make new folder to save things in:
  if (!dir.exists(paste0( "../output/", Sys.Date(), newfoldername))){
    dir.create(paste0("../output/", Sys.Date(), newfoldername), recursive = TRUE)
  }
  
  #initialize dataframes
  prevA_threshold_met_df <- data.frame()
  prevB_threshold_met_df <- data.frame()
  prevC_threshold_met_df <- data.frame()
  
  
  #iterate through the rows of the df
  for(row in c(1:nrow(paramdf))){
    dfrow = paramdf[row,]
    #run model
    modrun <- run_stochastic_brief_6_9(model = model, parms = generate_parms_sens(dfrow, threshold = threshold),
                                       runs = runs, dt = seq(0,365*years, 1))
    int_df <- make_prop_df(modrun[[1]], threshold= threshold, times = times)
    #add a column saying which parameter combination this is
    int_df_lab <- int_df |> mutate("parameter_comb" = row, "strategy" = strategy) |>
      select(step, parameter_comb, strategy, prop_runs_prevA_above,prop_runs_prevB_above, prop_runs_prevC_above)
    #save intermediate dataframe
    write_csv(int_df_lab, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "int_df_", row, ".csv"))
    #save seed list
    write_csv(as.data.frame(modrun[[2]]), paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "seed_list_", row, ".csv"))
    
    #update running dataframes
    prevA_threshold_met_df <- rbind(prevA_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevA_above))
    prevB_threshold_met_df <- rbind(prevB_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevB_above))
    prevC_threshold_met_df <- rbind(prevC_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevC_above))
    
    #print out intermediate running dataframes
    write_csv(prevA_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "intermediate_running_prevA_df",".csv"))
    write_csv(prevB_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "intermediate_running_prevB_df",".csv"))
    write_csv(prevC_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "intermediate_running_prevC_df",".csv"))
  }
  
  write_csv(prevA_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevAdf.csv"))
  write_csv(prevB_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevBdf.csv"))
  write_csv(prevC_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevCdf.csv"))
  
  return(list(prevA_threshold_met_df, prevB_threshold_met_df, prevC_threshold_met_df))
  
  
}

#run stochastic models but with seeds presupplied, output propdf
run_stochastic_brief_6_11 <- function(model, parms, seeds, runs, threshold, dt){
  mod <- do.call(model$new, parms)
  df <- data.frame()
  for(r in c(1:runs)){
    set.seed(seeds[r])
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r) |> select(step, run, prevA, prevB, prevC)
    df <- rbind(df, tmpdf)
  }
  
  prop_df <- make_prop_df(df, threshold = threshold, times = dt)
  return(prop_df)
}


run_stochastic_brief_6_12 <- function(model, parms, runs, dt, newfoldername){
  if (!dir.exists(paste0( "../output/", Sys.Date(), newfoldername))){
    dir.create(paste0("../output/", Sys.Date(), newfoldername), recursive = TRUE)
  }
  seeds <- sample(1:10^6, size = runs, replace = FALSE)
  write_csv(as.data.frame(seeds), paste0( "../output/", Sys.Date(), newfoldername, "/seeds.csv"))
  mod <- do.call(model$new, parms)
  df1 <- data.frame()
  df2 <- data.frame()
  for(r in c(1:runs)){
    print(r)
    set.seed(seeds[r])
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r) |> select(step, run, prevA, prevB, prevC,
                                                                     lock_threshold_A, lock_threshold_B, lock_threshold_C)
    tmpdf1 <- tmpdf |> select(step, run, prevA, prevB, prevC)
    tmpdf2 <- tmpdf |> select(step, run, lock_threshold_A, lock_threshold_B, lock_threshold_C)
    
    df1 <- rbind(df1, tmpdf1)
    df2 <- rbind(df2, tmpdf2)
    write_csv(df1, paste0( "../output/", Sys.Date(), newfoldername, "/prevdf.csv"))
    write_csv(df2, paste0( "../output/", Sys.Date(), newfoldername, "/lockdf.csv"))
  }
  return(list(df1, df2, seeds))
}

run_stochastic_brief_6_13 <- function(model, parms, runs, dt){
  seeds <- sample(1:10^6, size = runs, replace = FALSE)
  mod <- do.call(model$new, parms)
  df1 <- data.frame()
  df2 <- data.frame()
  for(r in c(1:runs)){
    print(r)
    set.seed(seeds[r])
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r) |> select(step, run, prevA, prevB, prevC,
                                                                     lock_threshold_A, lock_threshold_B, lock_threshold_C)
    tmpdf1 <- tmpdf |> select(step, run, prevA, prevB, prevC)
    tmpdf2 <- tmpdf |> select(step, run, lock_threshold_A, lock_threshold_B, lock_threshold_C)
    
    df1 <- rbind(df1, tmpdf1)
    df2 <- rbind(df2, tmpdf2)
  }
  return(list(df1, df2, seeds))
}

generate_parms_sens_6_13 <- function(dfrow, threshold){
  set.seed(06092025)
  #Define parameters -- from 6/6 calibration
  b <- 0.5167924 #transmission probability per partnership
  w_a <- unlist(dfrow[5])#rate of resistance to A emerging on treatment with A #this is a probability, over course of treatment
  w_b <- unlist(dfrow[6]) #rate of resistance to B emerging on treatment with B
  w_c <- unlist(dfrow[7]) #rate of resistance to C emerging on treatment with C
  T_s <- 1/(14.9435648) #time to treatment for symptomatic infection
  T_m <- 0.4150241/365  #screening rate (time to treatment for asymptomatic infection)
  T_sr <- T_s/3 #time to retreatment for symptomatic infection, if failure ##rethink this?
  kappa <- 0.90 #probability of retreatment after treatment failure, symptomatic
  f_a <- unlist(dfrow[8]) #relative fitness when resistant to A vs pan-susceptible 
  f_b <- unlist(dfrow[9]) #relative fitness when resistant to B vs pan-susceptible
  f_c <- unlist(dfrow[10]) #relative fitness when resistant to C vs pan-susceptible
  f_ab <- f_a*f_b #relative fitness when resistant to A&B vs pan-susceptible
  f_ac <- f_a*f_c #relative fitness when resistant to A&C vs pan-susceptible
  f_bc <- f_b*f_c #relative fitness when resistant to B&C vs pan-susceptible
  f_abc <- f_a*f_b*f_c #relative fitness when resistant to A,B&C vs pan-susceptible
  sigma <- 0.4936291  #probability of symptomatic infection
  N_risk <- 3 #number of risk strata
  d <- 1/(91.8422808) #natural recovery rate from infection #this is probably shorter, Barbee study
  pop <- 10^6 #pop size
  pop.p <- c(0.3, 0.6, 0.1) #relative size of each risk group; low, M, high
  c_min <- 1.3785266 #rate of partner change in lowest risk group
  activities <- c(1*c_min/365, 
                  5*c_min/365, 
                  20*c_min/365)  #sexual contacts per day
  epsilon <- 0.2357691 #mixing parameter (for non-assortativity)
  
  #risk stratified mixing 
  beta <- ((1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + epsilon*activities/(pop*pop.p)*diag(3))*b
  
  N_ini <- c(pop*pop.p[1],
             pop*pop.p[2],
             pop*pop.p[3])
  
  #distribute GC cases to have overall 3% prevalence
  x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)
  gc_lo <- round(N_ini[1]*x*0.029,0)
  gc_md <- round(N_ini[2]*x*0.154,0)
  gc_hi <- round(N_ini[3]*x*0.817,0)
  
  prev_symp <- 0.10
  Y_gc_lo <- gc_lo*prev_symp
  Z_gc_lo <- gc_lo*(1-prev_symp)
  Y_gc_md <- gc_md*prev_symp
  Z_gc_md <- gc_md*(1-prev_symp)
  Y_gc_hi <- gc_hi*prev_symp
  Z_gc_hi <- gc_hi*(1-prev_symp)
  
  #determine the starting resistances first
  all_starting_a <- c(rmultinom(n = 1, size = unlist(dfrow[2]), c(0.0029, 0.0154, 0.0817, 0.0261,0.1386,0.7353))) #dividing up probability of symptom status and risk group
  Y_a_ini <- all_starting_a[1:3]
  Z_a_ini <- all_starting_a[4:6]
  
  all_starting_b <- c(rmultinom(n = 1, size = unlist(dfrow[3]), c(0.0029, 0.0154, 0.0817, 0.0261,0.1386,0.7353))) #dividing up probability of symptom status and risk group
  Y_b_ini <- all_starting_b[1:3]
  Z_b_ini <- all_starting_b[4:6]
  
  all_starting_c <- c(rmultinom(n = 1, size = unlist(dfrow[4]), c(0.0029, 0.0154, 0.0817, 0.0261,0.1386,0.7353))) #dividing up probability of symptom status and risk group
  Y_c_ini <- all_starting_c[1:3]
  Z_c_ini <- all_starting_c[4:6]
  
  
  Y_s_ini <- round(c(Y_gc_lo, Y_gc_md, Y_gc_hi)) - Y_a_ini - Y_b_ini - Y_c_ini
  Z_s_ini <- round(c(Z_gc_lo, Z_gc_md, Z_gc_hi)) - Z_a_ini - Z_b_ini - Z_c_ini
  
  Y_ab_ini <- c(0,0,0)
  Z_ab_ini <- c(0,0,0)
  
  Y_ac_ini <- c(0,0,0)
  Z_ac_ini <- c(0,0,0)
  
  Y_bc_ini <- c(0,0,0)
  Z_bc_ini <- c(0,0,0)
  
  Y_abc_ini <- c(0,0,0)
  Z_abc_ini <- c(0,0,0)
  
  S_ini <- c(N_ini[1] - gc_lo,
             N_ini[2] - gc_md,
             N_ini[3] - gc_hi)
  
  N = sum(Y_a_ini) + sum(Z_a_ini) + sum(Y_b_ini) + sum(Z_b_ini) + sum(Y_c_ini) + sum(Z_c_ini) + sum(Y_s_ini) + sum(Z_s_ini) + sum(S_ini)
  
  return(list(S_ini = S_ini,Y_s_ini = Y_s_ini, Z_s_ini = Z_s_ini, Y_a_ini = Y_a_ini, Z_a_ini = Z_a_ini, Y_b_ini = Y_b_ini,
              Z_b_ini = Z_b_ini, Y_c_ini = Y_c_ini, Z_c_ini = Z_c_ini, Y_ab_ini = Y_ab_ini, Z_ab_ini = Z_ab_ini, Y_ac_ini = Y_ac_ini,
              Z_ac_ini = Z_ac_ini, Y_bc_ini = Y_bc_ini,Z_bc_ini = Z_bc_ini,Y_abc_ini = Y_abc_ini,Z_abc_ini = Z_abc_ini, beta = beta,
              w_a = w_a, w_b = w_b, w_c = w_c, T_s = T_s, T_m = T_m, T_sr = T_sr, kappa = kappa, f_a = f_a, f_b = f_b, f_c = f_c, 
              f_ab = f_ab, f_ac = f_ac, f_bc = f_bc, f_abc = f_abc, sigma = sigma, N_risk = N_risk, d = d, threshold = threshold))
  
}


make_lock_df <- function(df2, times){
  df_summary <- df2 |> group_by(run) |>
    mutate(A_out = lock_threshold_A >= 1, B_out = lock_threshold_B >= 1, C_out = lock_threshold_C >= 1) |>
    mutate(num_lost = A_out + B_out + C_out) |> 
    mutate(all3lost = (A_out + B_out + C_out == 3)) |>
    group_by(step) |>
    summarize(avg_drugs_lost = mean(num_lost), avg_all3_lost = mean(all3lost))
  df_summary_filtered <- df_summary |> filter(step %in% times)
  return(df_summary_filtered)
}



run_sens_analysis_6_13 <- function(model, strategy, threshold, paramdf, runs, years, times, newfoldername, prefix){
  #make new folder to save things in:
  if (!dir.exists(paste0( "../output/", Sys.Date(), newfoldername))){
    dir.create(paste0("../output/", Sys.Date(), newfoldername), recursive = TRUE)
  }
  
  #initialize dataframes
  prevA_threshold_met_df <- data.frame()
  prevB_threshold_met_df <- data.frame()
  prevC_threshold_met_df <- data.frame()
  lock_df <- data.frame()
  
  
  #iterate through the rows of the df
  for(row in c(1:nrow(paramdf))){
    dfrow = paramdf[row,]
    #run model
    modrun <- run_stochastic_brief_6_13(model = model, parms = generate_parms_sens_6_13(dfrow, threshold = threshold),
                                       runs = runs, dt = seq(0,365*years, 1))
    int_df <- make_prop_df(modrun[[1]], threshold= threshold, times = times)
    #add a column saying which parameter combination this is
    int_df_lab <- int_df |> mutate("parameter_comb" = row, "strategy" = strategy) |>
      select(step, parameter_comb, strategy, prop_runs_prevA_above,prop_runs_prevB_above, prop_runs_prevC_above)
    
    int_lock_df <- make_lock_df(modrun[[2]], times = times) |> mutate("parameter_comb" = row, "strategy" = strategy)
    #save intermediate dataframes
    write_csv(int_df_lab, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "int_df_", row, ".csv"))
    write_csv(int_lock_df,  paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "int_lock_df_", row, ".csv"))
    #save seed list
    write_csv(as.data.frame(modrun[[3]]), paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "seed_list_", row, ".csv"))
    
    #update running dataframes
    prevA_threshold_met_df <- rbind(prevA_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevA_above))
    prevB_threshold_met_df <- rbind(prevB_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevB_above))
    prevC_threshold_met_df <- rbind(prevC_threshold_met_df, select(int_df_lab, step, parameter_comb, strategy, prop_runs_prevC_above))
    lock_df <- rbind(lock_df, int_lock_df)
    
    #print out intermediate running dataframes
    write_csv(prevA_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "intermediate_running_prevA_df",".csv"))
    write_csv(prevB_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "intermediate_running_prevB_df",".csv"))
    write_csv(prevC_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "intermediate_running_prevC_df",".csv"))
    write_csv(lock_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, "intermediate_lock_df",".csv"))
  }
  
  write_csv(prevA_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevAdf.csv"))
  write_csv(prevB_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevBdf.csv"))
  write_csv(prevC_threshold_met_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "prevCdf.csv"))
  write_csv(lock_df, paste0("../output/", Sys.Date(), newfoldername, "/", prefix, strategy, "lock_df.csv"))
  
  return(list(prevA_threshold_met_df, prevB_threshold_met_df, prevC_threshold_met_df, lock_df))
  
  
}

run_stochastic_brief_6_27 <- function(model, parms, runs, dt, newfoldername){
  if (!dir.exists(paste0( "../output/baseline_parms_analysis/", Sys.Date(), newfoldername))){
    dir.create(paste0("../output/baseline_parms_analysis/", Sys.Date(), newfoldername), recursive = TRUE)
  }
  seeds <- sample(1:10^6, size = runs, replace = FALSE)
  write_csv(as.data.frame(seeds), paste0( "../output/baseline_parms_analysis/", Sys.Date(), newfoldername, "/seeds.csv"))
  mod <- do.call(model$new, parms)
  df1 <- data.frame()
  df2 <- data.frame()
  for(r in c(1:runs)){
    print(r)
    set.seed(seeds[r])
    tmpdf <- as.data.frame(mod$run(dt)) |> mutate(run = r) |> select(step, run, prevA, prevB, prevC,
                                                                     lock_threshold_A, lock_threshold_B, lock_threshold_C)
    tmpdf1 <- tmpdf |> select(step, run, prevA, prevB, prevC)
    tmpdf2 <- tmpdf |> select(step, run, lock_threshold_A, lock_threshold_B, lock_threshold_C)
    
    df1 <- rbind(df1, tmpdf1)
    df2 <- rbind(df2, tmpdf2)
    write_csv(df1, paste0( "../output/baseline_parms_analysis/", Sys.Date(), newfoldername, "/prevdf.csv"))
    write_csv(df2, paste0( "../output/baseline_parms_analysis/", Sys.Date(), newfoldername, "/lockdf.csv"))
  }
  return(list(df1, df2, seeds))
}

make_prop_plot_years_solid <- function(prop_df, threshold){
  prop_plot <- prop_df |> mutate(full_name = paste(strategy, name), model = "Stochastic") |>
    mutate(Years = step / 365) |>
    ggplot(aes(x = Years, y = value, color = full_name, linetype = full_name)) + geom_line() +
    theme_classic() +
    scale_color_manual(values = c("equal allocation prop_runs_prevA_above" = '#d73027',
                                  "seq prop_runs_prevA_above" = '#041c5a',
                                  "equal allocation prop_runs_prevB_above" = '#fc8d59',
                                  "seq prop_runs_prevB_above" = '#2171b5',
                                  "equal allocation prop_runs_prevC_above" = '#fdb863',
                                  "seq prop_runs_prevC_above" = '#6baed6'),
                       labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A",
                                  "Sequential, Drug B", "Sequential, Drug C"),
                       name = "Strategy and Drug Resistance Prevalence") +
    scale_linetype_manual(values = c("equal allocation prop_runs_prevA_above" = 1,
                                     "seq prop_runs_prevA_above" = 1,
                                     "equal allocation prop_runs_prevB_above" = 1,
                                     "seq prop_runs_prevB_above" = 1,
                                     "equal allocation prop_runs_prevC_above" = 1,
                                     "seq prop_runs_prevC_above" = 1),
                          labels = c("Equal allocation, Drug A", "Equal allocation, Drug B", "Equal Allocation, Drug C",  "Sequential, Drug A",
                                     "Sequential, Drug B", "Sequential, Drug C"),
                          name = "Strategy and Drug Resistance Prevalence") +
    ylab("Proportion of simulations") +  ylim(0,1)
  return(prop_plot)
}





make_heat_maps_from_fxn_output_yrs <- function(seq_A, seq_B, seq_C, ea_A, ea_B, ea_C, foldername){
  if (!dir.exists(paste0( "../figures/", Sys.Date(), foldername))){
    dir.create(paste0("../figures/", Sys.Date(), foldername), recursive = TRUE)
  }

  #data joining and manipulation

  seq_fulldf <- left_join(seq_A, seq_B) |> left_join(seq_C)
  seq_full_1000_df <- seq_fulldf |> filter(step  == 1000)
  seq_full_5000_df <- seq_fulldf |> filter(step  == 5000)
  seq_full_10000_df <- seq_fulldf |> filter(step  == 10000)

  seq_full_2yr_df <- seq_fulldf |> filter(step  == 730)
  seq_full_5yr_df <- seq_fulldf |> filter(step  == 1825)
  seq_full_10yr_df <- seq_fulldf |> filter(step  == 3650)

  ea_fulldf <- left_join(ea_A, ea_B) |> left_join(ea_C)
  ea_full_1000_df <- ea_fulldf |> filter(step  == 1000)
  ea_full_5000_df <- ea_fulldf |> filter(step  == 5000)
  ea_full_10000_df <- ea_fulldf |> filter(step  == 10000)

  ea_full_2yr_df <- ea_fulldf |> filter(step  == 730)
  ea_full_5yr_df <- ea_fulldf |> filter(step  == 1825)
  ea_full_10yr_df <- ea_fulldf |> filter(step  == 3650)

  difference_df_all <- rbind(seq_fulldf, ea_fulldf) |>
    pivot_wider(names_from = strategy, values_from = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |> mutate(prop_A_seq_minus_ea = prop_runs_prevA_above_sequential - `prop_runs_prevA_above_equal allocation`) |>
    mutate(prop_B_seq_minus_ea = prop_runs_prevB_above_sequential - `prop_runs_prevB_above_equal allocation`) |>
    mutate(prop_C_seq_minus_ea = prop_runs_prevC_above_sequential - `prop_runs_prevC_above_equal allocation`) |>
    select(step, parameter_comb, prop_A_seq_minus_ea, prop_B_seq_minus_ea, prop_C_seq_minus_ea)
  difference_df_1000 <- difference_df_all |> filter(step == 1000)
  difference_df_5000 <- difference_df_all |> filter(step == 5000)
  difference_df_10000 <- difference_df_all |> filter(step == 10000)

  difference_df_2yr <- difference_df_all |> filter(step == 730)
  difference_df_5yr <- difference_df_all |> filter(step == 1825)
  difference_df_10yr <- difference_df_all |> filter(step == 3650)

  #plotting seq only heatmaps

  seqplot_full_1000 <- seq_full_1000_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/seq_1000.tiff"), seqplot_full_1000, width = 4, height = 8, dpi = 900)

  seqplot_full_5000 <- seq_full_5000_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/seq_5000.tiff"), seqplot_full_5000, width = 4, height = 8, dpi = 900)

  seqplot_full_10000 <- seq_full_10000_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/seq_10000.tiff"), seqplot_full_10000, width = 4, height = 8, dpi = 900)

  seqplot_full_2yr <- seq_full_2yr_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/seq_2yr.tiff"), seqplot_full_2yr, width = 4, height = 8, dpi = 900)

  seqplot_full_5yr <- seq_full_5yr_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/seq_5yr.tiff"), seqplot_full_5yr, width = 4, height = 8, dpi = 900)

  seqplot_full_10yr <- seq_full_10yr_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/seq_10yr.tiff"), seqplot_full_10yr, width = 4, height = 8, dpi = 900)



  #plotting ea only heatmaps

  eaplot_full_1000 <- ea_full_1000_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/ea_1000.tiff"), eaplot_full_1000, width = 4, height = 8, dpi = 900)

  eaplot_full_5000 <- ea_full_5000_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/ea_5000.tiff"), eaplot_full_5000, width = 4, height = 8, dpi = 900)

  eaplot_full_10000 <- ea_full_10000_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/ea_10000.tiff"), eaplot_full_10000, width = 4, height = 8, dpi = 900)


  eaplot_full_2yr <- ea_full_2yr_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/ea_2yr.tiff"), eaplot_full_2yr, width = 4, height = 8, dpi = 900)

  eaplot_full_5yr <- ea_full_5yr_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/ea_5yr.tiff"), eaplot_full_5yr, width = 4, height = 8, dpi = 900)

  eaplot_full_10yr <- ea_full_10yr_df |>
    pivot_longer(cols = c(prop_runs_prevA_above, prop_runs_prevB_above, prop_runs_prevC_above)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient(low = "#fde0dd", high = '#d73027', limits = c(0, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ggtitle("Proportion of runs hitting 5% Resistance Prevalence Threshold") + ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/ea_10yr.tiff"), eaplot_full_10yr, width = 4, height = 8, dpi = 900)


  #plotting combined difference heat maps

  difference_1000 <- difference_df_1000 |>   pivot_longer(cols = c(prop_A_seq_minus_ea, prop_B_seq_minus_ea, prop_C_seq_minus_ea)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = '#08519c', limits = c(-1, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/diffmap_1000.tiff"), difference_1000, width = 4, height = 8)


  difference_5000 <- difference_df_5000 |>   pivot_longer(cols = c(prop_A_seq_minus_ea, prop_B_seq_minus_ea, prop_C_seq_minus_ea)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = '#08519c', limits = c(-1, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/diffmap_5000.tiff"), difference_5000, width = 4, height = 8, dpi = 900)

  difference_10000 <- difference_df_10000 |>   pivot_longer(cols = c(prop_A_seq_minus_ea, prop_B_seq_minus_ea, prop_C_seq_minus_ea)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = '#08519c', limits = c(-1, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/diffmap_10000.tiff"), difference_10000, width = 4, height = 8, dpi = 900)

  difference_2yr <- difference_df_2yr |>   pivot_longer(cols = c(prop_A_seq_minus_ea, prop_B_seq_minus_ea, prop_C_seq_minus_ea)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = '#08519c', limits = c(-1, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/diffmap_2yr.tiff"), difference_2yr, width = 4, height = 8, dpi = 900)


  difference_5yr <- difference_df_5yr |>   pivot_longer(cols = c(prop_A_seq_minus_ea, prop_B_seq_minus_ea, prop_C_seq_minus_ea)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = '#08519c', limits = c(-1, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/diffmap_5yr.tiff"), difference_5yr, width = 4, height = 8, dpi = 900)

  difference_10yr <- difference_df_10yr |>   pivot_longer(cols = c(prop_A_seq_minus_ea, prop_B_seq_minus_ea, prop_C_seq_minus_ea)) |>
    ggplot(aes(x = name, y = parameter_comb, fill = value)) + geom_tile(color = "white") +
    # geom_text(aes(label = value), color = "black", size = 3) +
    scale_fill_gradient2(low = "#d73027", mid = "white", high = '#08519c', limits = c(-1, 1), name = "Proportion of runs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Drug A", "Drug B", "Drug C")) +
    ylab("LHS Parameter Combination") + xlab("")

  ggsave(paste0("../figures/", Sys.Date(), foldername, "/diffmap_10yr.tiff"), difference_10yr, width = 4, height = 8, dpi = 900)

  return(list(seq_fulldf, ea_fulldf, difference_df_all))

}
