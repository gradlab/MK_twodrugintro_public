library(readr)
library(lhs)


#start with parameters which the model SHOULDN'T be sensitive to
#these are: sigma, b, epsilon, c_min, d, T_s, T_m , do the factor the T_s is divdied by to get T_sr
#sigma, b, epsilon are probabilties so keep as they are
#c_min is partner change rate; inverse of avg duration of partnership; make range from new partnership every day to avg 5 years; so 0.2 to 365
#d is avg duration of natural recovery, the d parameter is actually 1 / (#days), but just draw number of days between 1 and 365
#T_s is treatment rate if initial treatment is successful; calibrated value is about 14 days (2 weeks), do 1 day to 90 days same issue as d
#T_m is essentially a probabiltiy so keep as is
#T_sr right now is number of times longer than T_s it is if retreated, try all the way from 1 to 10
#kappa is the probability of being retreated once failing treatment; it doens't end up meaning precisely this but it is used in generating rate of recovery from retreatment

unsensitive_lhs <- as.data.frame(randomLHS(1000, 10))
colnames(unsensitive_lhs) <- c('comb', 'sigma', 'b', 'epsilon', 'c_min', 'd', 'T_s', 'T_m', 'T_sr', 'kappa')
unsensitive_lhs_filled <- unsensitive_lhs |> mutate(comb = c(1:1000),
                                                    sigma = sigma,
                                                    b = b,
                                                    epsilon = epsilon,
                                                    c_min = 0.2 + c_min*(365-0.2),
                                                    d = 1 + d*364,
                                                    T_s = 1 + T_s*89,
                                                    T_m = T_m,
                                                    T_sr = 1 + T_sr*9,
                                                    kappa = kappa)

#now start with parameters which I think the model might be sensitive to
#use log uniform scale for small rates
#parameters: a_res_start, b_res_start, c_res_start, w_a, w_b, w_c, f_a, f_b, f_c
sensitive_lhs <- as.data.frame(randomLHS(1000, 10))
colnames(sensitive_lhs) <- c('comb', 'a_res_start', 'b_res_start', 'c_res_start', 'w_a', 'w_b', 'w_c', 'f_a', 'f_b', "f_c")
sensitive_lhs_filled<- sensitive_lhs |>
  mutate(comb = c(1:1000),
         a_res_start = round(0 + a_res_start*10),
         b_res_start = round(0 + b_res_start*10),
         c_res_start = round(0 + c_res_start*10),
         w_a = 10^(-7 + w_a*4),
         w_b = 10^(-7 + w_b*4), #fixed 6-11
         w_c = 10^(-7 + w_c*4), #fixed 6-11
         f_a = 0.5 + f_a*0.5, #fixed 6-11
         f_b = 0.5 + f_b*0.5,
         f_c = 0.5 + f_c*0.5)

write_csv(unsensitive_lhs_filled, "../output/6_7_unsensitive_LHS_df.csv")
write_csv(sensitive_lhs_filled, "../output/6_11_sensitive_LHS_df_fixed.csv")
