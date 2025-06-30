#adapt stochastic models to discrete time with binary switches and make stochastic

#old model that I'm working off of
sequential_1_18_retx_mod <- odin({
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) + f_a*(Y_a[i]+Z_a[i]) + f_b*(Y_b[i] + Z_b[i]) + f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + f_ac*(Y_ac[i] + Z_ac[i]) + f_bc*(Y_bc[i] + Z_bc[i]) + f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out equation for S
  deriv(S[]) <-  -sum(BI_S[,i]) * S[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] +
    E_b*(1 - w_b)*T_s*Y_a[i] + E_c*(1 - w_c)*T_s*Y_a[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] +
    E_c*(1 - w_c)*T_s*Y_b[i] + E_a*(1 - w_a)*T_s*Y_b[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] +
    E_b*(1 - w_b)*T_s*Y_c[i] + E_a*(1 - w_a)*T_s*Y_c[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] +
    E_c*(1 - w_c)*T_s*Y_ab[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] +
    E_b*(1 - w_b)*T_s*Y_ac[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] +
    E_a*(1 - w_a)*T_s*Y_bc[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] +
    E_b*(1 - w_b)*T_m*Z_a[i] + E_c*(1 - w_c)*T_m*Z_a[i] +
    E_a*(1 - w_a)*T_m*Z_b[i] + E_c*(1 - w_c)*T_m*Z_b[i] +
    E_b*(1 - w_b)*T_m*Z_c[i] + E_a*(1 - w_a)*T_m*Z_c[i] +
    E_c*(1 - w_c)*T_m*Z_ab[i] + E_b*(1 - w_b)*T_m*Z_ac[i] + E_a*(1 - w_a)*T_m*Z_bc[i] +
    # #add last resort treatments that should go on and off depending on which drugs are being used
    # E_a*kappa*T_sr*Y_a[i] + E_b*kappa*T_sr*Y_b[i] + E_b*kappa*T_sr*Y_ab[i] + E_c*kappa*T_sr*Y_c[i] + E_c*kappa*T_sr*Y_ac[i] + E_c*kappa*T_sr*Y_bc[i] +
    v_c*kappa*T_sr*Y_c[i] + v_ac*kappa*T_sr*Y_ac[i] + v_bc*kappa*T_sr*Y_bc[i] + #added LRT 1/18
    v_abc*kappa*T_sr*Y_abc[i] + 
    (d + r)*(N[i] - S[i]) 
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  
  #write out equations for dY_s and dZ_s
  deriv(Y_s[]) <- sigma*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_s*Y_s[i] - E_b*w_b*T_s*Y_s[i] - E_c*w_c*T_s*Y_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] - (d + r)*Y_s[i]
  deriv(Z_s[]) <- (1 - sigma)*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_m*Z_s[i] - E_b*w_b*T_m*Z_s[i] - E_c*w_c*T_m*Z_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] - (d + r)*Z_s[i]
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out equations for dY_a and dZ_a
  deriv(Y_a[]) <- sigma*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_a[i] - E_c*(1-w_c)*T_s*Y_a[i] - 
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] - 
    E_b*w_b*T_s*Y_a[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] - E_c*w_c*T_s*Y_a[i] - E_a*L_c_a*w_c*kappa*T_s*Y_a[i] - 
    E_a*kappa*T_sr*Y_a[i] - (d + r)*Y_a[i] #added last resort treatment only active when using  A (only A scenario)
  
  deriv(Z_a[]) <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_a[i] - E_c*(1-w_c)*T_m*Z_a[i] - 
    E_b*w_b*T_m*Z_a[i] - E_c*w_c*T_m*Z_a[i] - (d + r)*Z_a[i]
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out equations for dY_b and dZ_b
  deriv(Y_b[]) <- sigma*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_s*Y_s[i] - E_c*(1 - w_c)*T_s*Y_b[i] - E_a*(1-w_a)*T_s*Y_b[i] - 
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] - 
    E_a*w_a*T_s*Y_b[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] - E_c*w_c*T_s*Y_b[i] - E_b*L_c_b*w_c*kappa*T_s*Y_b[i] - 
    E_b*kappa*T_sr*Y_b[i] - (d + r)*Y_b[i] #added last resort treatment when E_b is active
  
  deriv(Z_b[]) <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_m*Z_s[i] - E_c*(1 - w_c)*T_m*Z_b[i] - E_a*(1-w_a)*T_m*Z_b[i] - 
    E_c*w_c*T_m*Z_b[i] - E_a*w_a*T_m*Z_b[i] - (d + r)*Z_b[i]
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out equations for dY_c and dZ_c
  deriv(Y_c[]) <- sigma*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_c[i] - E_a*(1-w_a)*T_s*Y_c[i] - 
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] - 
    E_a*w_a*T_s*Y_c[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] - E_b*w_b*T_s*Y_c[i] - E_c*L_b_c*w_b*kappa*T_s*Y_c[i] - 
    v_c*kappa*T_sr*Y_c[i] - (d + r)*Y_c[i] #added / modified last resort treatment 1/18
  deriv(Z_c[]) <- (1 - sigma)*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_c[i] - E_a*(1-w_a)*T_m*Z_c[i] - 
    E_a*w_a*T_m*Z_c[i] - E_b*w_b*T_m*Z_c[i] - (d + r)*Z_c[i]
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out equations for dY_ab and dZ_ab
  deriv(Y_ab[]) <- sigma*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_s*Y_a[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_b[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] -
    E_c*(1 - w_c)*T_s*Y_ab[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] - 
    E_c*w_c*T_s*Y_ab[i] - E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] - E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] -
    (d+r)*Y_ab[i] #removed last resort treatment term 1/18 
  deriv(Z_ab[]) <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_m*Z_a[i] + E_a*w_a*T_m*Z_b[i] - E_c*(1 - w_c)*T_m*Z_ab[i] - E_c*w_c*T_m*Z_ab[i] - (d + r)*Z_ab[i]
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out equations for dY_ac and dZ_ac
  deriv(Y_ac[]) <- sigma*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_s*Y_a[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_c[i] + E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] -
    E_b*(1 - w_b)*T_s*Y_ac[i] - E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] - 
    E_b*w_b*T_s*Y_ac[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] - E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] -
    v_ac*kappa*T_sr*Y_ac[i] - (d+r)*Y_ac[i] #added / modified last resort treatment 1/18
  deriv(Z_ac[]) <- (1-sigma)*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_m*Z_a[i] + E_a*w_a*T_m*Z_c[i] - E_b*(1 - w_b)*T_m*Z_ac[i] - E_b*w_b*T_m*Z_ac[i] - (d + r)*Z_ac[i]
  
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out equations for dY_ac and dZ_ac
  deriv(Y_bc[]) <- sigma*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_s*Y_b[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_b[i] + E_b*w_b*T_s*Y_c[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_c[i] -
    E_a*(1 - w_a)*T_s*Y_bc[i] - E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] - E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] - 
    E_a*w_a*T_s*Y_bc[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] -
    v_bc*kappa*T_sr*Y_bc[i] - (d+r)*Y_bc[i] #added/modified LRT term 1/18
  deriv(Z_bc[]) <- (1-sigma)*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_m*Z_b[i] + E_b*w_b*T_m*Z_c[i] - 
    E_a*(1 - w_a)*T_m*Z_bc[i]  - E_a*w_a*T_m*Z_bc[i] - (d + r)*Z_bc[i] 
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  #write out equations for dY_abc and dZ_abc
  deriv(Y_abc[]) <- sigma*f_abc*sum(BI_YZabc[,i])*S[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] + E_c*w_c*T_s*Y_ab[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] +
    E_b*w_b*T_s*Y_ac[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] + E_a*w_a*T_s*Y_bc[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] +
    E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] - v_abc*kappa*T_sr*Y_abc[i] - (d + r)*Y_abc[i] #add LRT prob v_abc 1/18
  deriv(Z_abc[]) <- (1 - sigma)*f_abc*sum(BI_YZabc[,i])*S[i] + E_c*w_c*T_m*Z_ab[i] + E_b*w_b*T_m*Z_ac[i] + E_a*w_a*T_m*Z_bc[i] - (d + r)*Z_abc[i]
  
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  #adding this 5/27
  sum_N <- sum(N)
  
  
  #prevalences: fixed as of 9/25
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  denom <- sum(N) - sum(S)
  
  
  #set lock values:
  #think these are vestigial now
  deriv(lockA) <- ((prevA >= threshold)*(1-lockA)) 
  deriv(lockB) <- ((prevB >= threshold)* (1-lockB))
  
  # Smooth transitions of probabilities of treatment with drug
  Q_A <- transition_timeA + duration / 2
  tempA <- (t / Q_A-1) / tune_parameter
  E_a <- 1 -  1 / (1 + exp(-tempA))
  Q_C <- transition_timeC + duration/2
  tempC <- (t/Q_C-1)/ tune_parameter
  E_c <- 1 / (1 + exp(-tempC))
  E_b <- 1  - E_a - E_c
  
  #Set retreatment terms; for 1-18 sequential model are fixed
  L_b_a <- 0.5
  L_b_c <- 0
  L_c_a <- 0.5
  L_a_b <- 0
  L_a_c <- 0
  L_c_b <- 1
  
  #last resort treatment terms
  v_abc <- 1
  v_bc <- 1 - E_a 
  v_ac <- E_c
  v_c <- E_c
  
  
  #normalize these so that people always get treated with one drug (sum is always 1)
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lockA) <- 0
  initial(lockB) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  transition_timeA <- user()
  transition_timeC <- user()
  duration <- user()
  tune_parameter <- user()
  b <- user()
  # E_a <- user()
  # E_b <- user()
  # E_c <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  r <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(denom) <- denom
  output(A_transition_time) <- transition_timeA
  output(L_b_a) <- L_b_a
  output(L_b_c) <- L_b_c
  output(L_c_a) <- L_c_a
  output(L_a_b) <- L_a_b
  output(L_a_c) <- L_a_c
  output(L_c_b) <- L_c_b 
  output(v_abc) <- v_abc
  output(v_bc) <- v_bc
  output(v_ac) <- v_ac
  output(v_c) <- v_c
  output(sum_N) <- sum_N
  
})

#start with most up-to-date sequential model and make it discrete
sequential_discrete_5_26 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) + f_a*(Y_a[i]+Z_a[i]) + f_b*(Y_b[i] + Z_b[i]) + f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + f_ac*(Y_ac[i] + Z_ac[i]) + f_bc*(Y_bc[i] + Z_bc[i]) + f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  
  #write out equation for S
  update(S[]) <- S[i] -sum(BI_S[,i]) * S[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] +
    E_b*(1 - w_b)*T_s*Y_a[i] + E_c*(1 - w_c)*T_s*Y_a[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] +
    E_c*(1 - w_c)*T_s*Y_b[i] + E_a*(1 - w_a)*T_s*Y_b[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] +
    E_b*(1 - w_b)*T_s*Y_c[i] + E_a*(1 - w_a)*T_s*Y_c[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] +
    E_c*(1 - w_c)*T_s*Y_ab[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] +
    E_b*(1 - w_b)*T_s*Y_ac[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] +
    E_a*(1 - w_a)*T_s*Y_bc[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] +
    E_b*(1 - w_b)*T_m*Z_a[i] + E_c*(1 - w_c)*T_m*Z_a[i] +
    E_a*(1 - w_a)*T_m*Z_b[i] + E_c*(1 - w_c)*T_m*Z_b[i] +
    E_b*(1 - w_b)*T_m*Z_c[i] + E_a*(1 - w_a)*T_m*Z_c[i] +
    E_c*(1 - w_c)*T_m*Z_ab[i] + E_b*(1 - w_b)*T_m*Z_ac[i] + E_a*(1 - w_a)*T_m*Z_bc[i] +
    # #add last resort treatments that should go on and off depending on which drugs are being used
    v_c*kappa*T_sr*Y_c[i] + v_ac*kappa*T_sr*Y_ac[i] + v_bc*kappa*T_sr*Y_bc[i] + #added LRT 1/18
    v_abc*kappa*T_sr*Y_abc[i] + 
    (d)*(N[i] - S[i]) 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + sigma*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_s*Y_s[i] - E_b*w_b*T_s*Y_s[i] - E_c*w_c*T_s*Y_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] - (d )*Y_s[i]
  update(Z_s[]) <- Z_s[i] + (1 - sigma)*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_m*Z_s[i] - E_b*w_b*T_m*Z_s[i] - E_c*w_c*T_m*Z_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] - (d)*Z_s[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + sigma*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_a[i] - E_c*(1-w_c)*T_s*Y_a[i] - 
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] - 
    E_b*w_b*T_s*Y_a[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] - E_c*w_c*T_s*Y_a[i] - E_a*L_c_a*w_c*kappa*T_s*Y_a[i] - 
    E_a*kappa*T_sr*Y_a[i] - (d)*Y_a[i] #added last resort treatment only active when using  A (only A scenario)
  
  update(Z_a[]) <- Z_a[i] + (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_a[i] - E_c*(1-w_c)*T_m*Z_a[i] - 
    E_b*w_b*T_m*Z_a[i] - E_c*w_c*T_m*Z_a[i] - (d)*Z_a[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + sigma*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_s*Y_s[i] - E_c*(1 - w_c)*T_s*Y_b[i] - E_a*(1-w_a)*T_s*Y_b[i] - 
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] - 
    E_a*w_a*T_s*Y_b[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] - E_c*w_c*T_s*Y_b[i] - E_b*L_c_b*w_c*kappa*T_s*Y_b[i] - 
    E_b*kappa*T_sr*Y_b[i] - (d)*Y_b[i] #added last resort treatment when E_b is active
  
  update(Z_b[]) <- Z_b[i] + (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_m*Z_s[i] - E_c*(1 - w_c)*T_m*Z_b[i] - E_a*(1-w_a)*T_m*Z_b[i] - 
    E_c*w_c*T_m*Z_b[i] - E_a*w_a*T_m*Z_b[i] - (d)*Z_b[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + sigma*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_c[i] - E_a*(1-w_a)*T_s*Y_c[i] - 
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] - 
    E_a*w_a*T_s*Y_c[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] - E_b*w_b*T_s*Y_c[i] - E_c*L_b_c*w_b*kappa*T_s*Y_c[i] - 
    v_c*kappa*T_sr*Y_c[i] - (d)*Y_c[i] #added / modified last resort treatment 1/18
  
  update(Z_c[]) <- Z_c[i] + (1 - sigma)*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_c[i] - E_a*(1-w_a)*T_m*Z_c[i] - 
    E_a*w_a*T_m*Z_c[i] - E_b*w_b*T_m*Z_c[i] - (d)*Z_c[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + sigma*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_s*Y_a[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_b[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] -
    E_c*(1 - w_c)*T_s*Y_ab[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] - 
    E_c*w_c*T_s*Y_ab[i] - E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] - E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] -
    (d)*Y_ab[i] #removed last resort treatment term 1/18 
  
  update(Z_ab[]) <- Z_ab[i] + (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_m*Z_a[i] + E_a*w_a*T_m*Z_b[i] - E_c*(1 - w_c)*T_m*Z_ab[i] - E_c*w_c*T_m*Z_ab[i] - (d)*Z_ab[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + sigma*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_s*Y_a[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_c[i] + E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] -
    E_b*(1 - w_b)*T_s*Y_ac[i] - E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] - 
    E_b*w_b*T_s*Y_ac[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] - E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] -
    v_ac*kappa*T_sr*Y_ac[i] - (d)*Y_ac[i] #added / modified last resort treatment 1/18
  
  update(Z_ac[]) <- Z_ac[i] + (1-sigma)*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_m*Z_a[i] + E_a*w_a*T_m*Z_c[i] - E_b*(1 - w_b)*T_m*Z_ac[i] - E_b*w_b*T_m*Z_ac[i] - (d)*Z_ac[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + sigma*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_s*Y_b[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_b[i] + E_b*w_b*T_s*Y_c[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_c[i] -
    E_a*(1 - w_a)*T_s*Y_bc[i] - E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] - E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] - 
    E_a*w_a*T_s*Y_bc[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] -
    v_bc*kappa*T_sr*Y_bc[i] - (d)*Y_bc[i] #added/modified LRT term 1/18
  
  update(Z_bc[]) <- Z_bc[i] + (1-sigma)*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_m*Z_b[i] + E_b*w_b*T_m*Z_c[i] - 
    E_a*(1 - w_a)*T_m*Z_bc[i]  - E_a*w_a*T_m*Z_bc[i] - (d)*Z_bc[i] 

  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + sigma*f_abc*sum(BI_YZabc[,i])*S[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] + E_c*w_c*T_s*Y_ab[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] +
    E_b*w_b*T_s*Y_ac[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] + E_a*w_a*T_s*Y_bc[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] +
    E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] - v_abc*kappa*T_sr*Y_abc[i] - (d)*Y_abc[i] #add LRT prob v_abc 1/18
  
  update(Z_abc[]) <- Z_abc[i] + (1 - sigma)*f_abc*sum(BI_YZabc[,i])*S[i] + E_c*w_c*T_m*Z_ab[i] + E_b*w_b*T_m*Z_ac[i] + E_a*w_a*T_m*Z_bc[i] - (d)*Z_abc[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  denom <- sum(N) - sum(S)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  # Binary switches between treatment
  E_a <- (lock_threshold_A < 1)
  E_c <- 1 - (lock_threshold_B < 1)
  E_b <- 1  - E_a - E_c
  
  #Set retreatment terms; for 1-18 sequential model are fixed
  L_b_a <- 0.5
  L_b_c <- 0
  L_c_a <- 0.5
  L_a_b <- 0
  L_a_c <- 0
  L_c_b <- 1
  
  #last resort treatment terms
  v_abc <- 1
  v_bc <- 1 - E_a 
  v_ac <- E_c
  v_c <- E_c
  
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(denom) <- denom
  output(L_b_a) <- L_b_a
  output(L_b_c) <- L_b_c
  output(L_c_a) <- L_c_a
  output(L_a_b) <- L_a_b
  output(L_a_c) <- L_a_c
  output(L_c_b) <- L_c_b 
  output(v_abc) <- v_abc
  output(v_bc) <- v_bc
  output(v_ac) <- v_ac
  output(v_c) <- v_c
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(sum_N) <- sum_N
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  
})

#old EA model to work off of
ea_1_18_retx_mod <- odin({
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) + f_a*(Y_a[i]+Z_a[i]) + f_b*(Y_b[i] + Z_b[i]) + f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + f_ac*(Y_ac[i] + Z_ac[i]) + f_bc*(Y_bc[i] + Z_bc[i]) + f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out equation for S
  deriv(S[]) <-  -sum(BI_S[,i]) * S[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] +
    E_b*(1 - w_b)*T_s*Y_a[i] + E_c*(1 - w_c)*T_s*Y_a[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] +
    E_c*(1 - w_c)*T_s*Y_b[i] + E_a*(1 - w_a)*T_s*Y_b[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] +
    E_b*(1 - w_b)*T_s*Y_c[i] + E_a*(1 - w_a)*T_s*Y_c[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] +
    E_c*(1 - w_c)*T_s*Y_ab[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] +
    E_b*(1 - w_b)*T_s*Y_ac[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] +
    E_a*(1 - w_a)*T_s*Y_bc[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] +
    E_b*(1 - w_b)*T_m*Z_a[i] + E_c*(1 - w_c)*T_m*Z_a[i] +
    E_a*(1 - w_a)*T_m*Z_b[i] + E_c*(1 - w_c)*T_m*Z_b[i] +
    E_b*(1 - w_b)*T_m*Z_c[i] + E_a*(1 - w_a)*T_m*Z_c[i] +
    E_c*(1 - w_c)*T_m*Z_ab[i] + E_b*(1 - w_b)*T_m*Z_ac[i] + E_a*(1 - w_a)*T_m*Z_bc[i] +
    # #add last resort treatments that should go on and off depending on which drugs are being used
    # E_a*kappa*T_sr*Y_a[i] + E_b*kappa*T_sr*Y_b[i] + E_b*kappa*T_sr*Y_ab[i] + E_c*kappa*T_sr*Y_c[i] + E_c*kappa*T_sr*Y_ac[i] + E_c*kappa*T_sr*Y_bc[i] +
    v_c*kappa*T_sr*Y_c[i] + v_ac*kappa*T_sr*Y_ac[i] + v_bc*kappa*T_sr*Y_bc[i] + #added LRT 1/18
    v_abc*kappa*T_sr*Y_abc[i] + 
    (d + r)*(N[i] - S[i]) 
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  
  #write out equations for dY_s and dZ_s
  deriv(Y_s[]) <- sigma*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_s*Y_s[i] - E_b*w_b*T_s*Y_s[i] - E_c*w_c*T_s*Y_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] - (d + r)*Y_s[i]
  deriv(Z_s[]) <- (1 - sigma)*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_m*Z_s[i] - E_b*w_b*T_m*Z_s[i] - E_c*w_c*T_m*Z_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] - (d + r)*Z_s[i]
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out equations for dY_a and dZ_a
  deriv(Y_a[]) <- sigma*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_a[i] - E_c*(1-w_c)*T_s*Y_a[i] - 
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] - 
    E_b*w_b*T_s*Y_a[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] - E_c*w_c*T_s*Y_a[i] - E_a*L_c_a*w_c*kappa*T_s*Y_a[i] - 
    E_a*kappa*T_sr*Y_a[i] - (d + r)*Y_a[i] #added last resort treatment only active when using  A (only A scenario)
  
  deriv(Z_a[]) <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_a[i] - E_c*(1-w_c)*T_m*Z_a[i] - 
    E_b*w_b*T_m*Z_a[i] - E_c*w_c*T_m*Z_a[i] - (d + r)*Z_a[i]
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out equations for dY_b and dZ_b
  deriv(Y_b[]) <- sigma*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_s*Y_s[i] - E_c*(1 - w_c)*T_s*Y_b[i] - E_a*(1-w_a)*T_s*Y_b[i] - 
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] - 
    E_a*w_a*T_s*Y_b[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] - E_c*w_c*T_s*Y_b[i] - E_b*L_c_b*w_c*kappa*T_s*Y_b[i] - 
    E_b*kappa*T_sr*Y_b[i] - (d + r)*Y_b[i] #added last resort treatment when E_b is active
  
  deriv(Z_b[]) <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_m*Z_s[i] - E_c*(1 - w_c)*T_m*Z_b[i] - E_a*(1-w_a)*T_m*Z_b[i] - 
    E_c*w_c*T_m*Z_b[i] - E_a*w_a*T_m*Z_b[i] - (d + r)*Z_b[i]
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out equations for dY_c and dZ_c
  deriv(Y_c[]) <- sigma*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_c[i] - E_a*(1-w_a)*T_s*Y_c[i] - 
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] - 
    E_a*w_a*T_s*Y_c[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] - E_b*w_b*T_s*Y_c[i] - E_c*L_b_c*w_b*kappa*T_s*Y_c[i] - 
    v_c*kappa*T_sr*Y_c[i] - (d + r)*Y_c[i] #modified last resort treatment when C is active 1/18 with v_c term
  deriv(Z_c[]) <- (1 - sigma)*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_c[i] - E_a*(1-w_a)*T_m*Z_c[i] - 
    E_a*w_a*T_m*Z_c[i] - E_b*w_b*T_m*Z_c[i] - (d + r)*Z_c[i]
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out equations for dY_ab and dZ_ab
  deriv(Y_ab[]) <- sigma*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_s*Y_a[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_b[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] -
    E_c*(1 - w_c)*T_s*Y_ab[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] - 
    E_c*w_c*T_s*Y_ab[i] - E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] - E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] -
    (d+r)*Y_ab[i] #removed responsive LRT 1/18
  deriv(Z_ab[]) <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_m*Z_a[i] + E_a*w_a*T_m*Z_b[i] - E_c*(1 - w_c)*T_m*Z_ab[i] - E_c*w_c*T_m*Z_ab[i] - (d + r)*Z_ab[i]
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out equations for dY_ac and dZ_ac
  deriv(Y_ac[]) <- sigma*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_s*Y_a[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_c[i] + E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] -
    E_b*(1 - w_b)*T_s*Y_ac[i] - E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] - 
    E_b*w_b*T_s*Y_ac[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] - E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] -
    v_ac*kappa*T_sr*Y_ac[i] - (d+r)*Y_ac[i] #modified last resort treatment with v_ac 1/18
  deriv(Z_ac[]) <- (1-sigma)*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_m*Z_a[i] + E_a*w_a*T_m*Z_c[i] - E_b*(1 - w_b)*T_m*Z_ac[i] - E_b*w_b*T_m*Z_ac[i] - (d + r)*Z_ac[i]
  
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out equations for dY_ac and dZ_ac
  deriv(Y_bc[]) <- sigma*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_s*Y_b[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_b[i] + E_b*w_b*T_s*Y_c[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_c[i] -
    E_a*(1 - w_a)*T_s*Y_bc[i] - E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] - E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] - 
    E_a*w_a*T_s*Y_bc[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] -
    v_bc*kappa*T_sr*Y_bc[i] - (d+r)*Y_bc[i] #added LRT v_bc 1/18 #FIXED 5/26
  
  deriv(Z_bc[]) <- (1-sigma)*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_m*Z_b[i] + E_b*w_b*T_m*Z_c[i] - 
    E_a*(1 - w_a)*T_m*Z_bc[i] - E_a*w_a*T_m*Z_bc[i]   -  (d + r)*Z_bc[i]  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  #write out equations for dY_abc and dZ_abc
  deriv(Y_abc[]) <- sigma*f_abc*sum(BI_YZabc[,i])*S[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] + E_c*w_c*T_s*Y_ab[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] +
    E_b*w_b*T_s*Y_ac[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] + E_a*w_a*T_s*Y_bc[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] +
    E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] - v_abc*kappa*T_sr*Y_abc[i] - (d + r)*Y_abc[i] #FOUND MISTAKE 9/19 #1/18 adding in v_abc LRT term
  deriv(Z_abc[]) <- (1 - sigma)*f_abc*sum(BI_YZabc[,i])*S[i] + E_c*w_c*T_m*Z_ab[i] + E_b*w_b*T_m*Z_ac[i] + E_a*w_a*T_m*Z_bc[i] - (d + r)*Z_abc[i]
  
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  
  #prevalences: fixed as of 9/25
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  denom <- sum(N) - sum(S)
  
  
  #set lock values: #probably don't need these? 
  deriv(lockA) <- ((prevA >= threshold)*(1-lockA)) 
  deriv(lockB) <- ((prevB >= threshold)* (1-lockB))
  
  # Smooth transitions of probabilities of treatment with drug
  
  #need A transition that goes from 0.33 down to 0 at first changepoint
  Q_A <- transition_timeA + duration / 2
  tempA <- (t / Q_A-1) / tune_parameter
  E_a <- startA -  startA / (1 + exp(-tempA))
  
  
  #make control1, which goes from 1/3 to 1/2 and stays there
  Q_con1 <- transition_timeA + duration / 2
  tempB <- (t / Q_con1 -1) / tune_parameter
  control1 <- startB + (midB - startB) / (1 + exp(-tempB))
  
  #make control2, which goes from 0 to 0.5 when C turns on 
  Q_con2 <- transition_timeC + duration / 2
  tempC <- (t / Q_con2 - 1) / tune_parameter
  control2 <- midB / (1 + exp(-tempC))
  
  #set E_b and E_c based on the controls
  E_b <- control1 - control2
  E_c <- control1 + control2
  
  #retreatment probabilities
  #smooth transitions of probabilities of retreatment 
  #retx1 takes care of transitions that occur when A phases out
  Q_retx1 <- transition_timeA + duration/2
  temp_retx1 <- (t/Q_retx1 -1)/ tune_parameter
  retx1 <-  0.5 + 0.5 / (1 + exp(-temp_retx1))
  L_b_a <- 0.5
  L_c_a <- 0.5
  L_c_b <- retx1
  L_a_b <- 1 - retx1
  L_a_c <- 1 - retx1
  #retx2 needed to take care of the 1 transition that occurs when B phases out
  Q_retx2 <- transition_timeC + duration/2
  temp_retx2 <- (t/Q_retx2 - 1) / tune_parameter
  retx2 <- 1/(1 + exp(-temp_retx2))
  L_b_c <- retx1 - retx2
  
  #last resort treatment terms
  v_abc <- 1
  v_bc <- 1/ ( 1 + exp(-temp_retx1))
  v_ac <- retx2
  v_c <- retx2
  
  
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lockA) <- 0
  initial(lockB) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  transition_timeA <- user() #this will be calculated by the model
  startA <- user() #this should start at 1/3 under equal allocation 
  #  transition_timeB <- user() #this will be calculated by the model
  startB <- user() #this should start at 1/3 under equal allocation
  midB <- user() #this should be 0.5 in the equal allocation model
  transition_timeC <- user()
  duration <- user()
  tune_parameter <- user()
  b <- user()
  # E_a <- user()
  # E_b <- user()
  # E_c <- user()
  # L_b_a <- user()
  # L_b_c <- user()
  # L_c_a <- user()
  # L_a_b <- user()
  # L_a_c <- user()
  # L_c_b <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  r <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(denom) <- denom
  output(A_transition_time) <- transition_timeA
  output(L_b_a) <- L_b_a
  output(L_b_c) <- L_b_c
  output(L_c_a) <- L_c_a
  output(L_a_b) <- L_a_b
  output(L_a_c) <- L_a_c
  output(L_c_b) <- L_c_b
  output(v_abc) <- v_abc
  output(v_bc) <- v_bc
  output(v_ac) <- v_ac
  output(v_c) <- v_c
})

ea_discrete_5_26 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) + f_a*(Y_a[i]+Z_a[i]) + f_b*(Y_b[i] + Z_b[i]) + f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + f_ac*(Y_ac[i] + Z_ac[i]) + f_bc*(Y_bc[i] + Z_bc[i]) + f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  
  #write out equation for S
  update(S[]) <- S[i] -sum(BI_S[,i]) * S[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] +
    E_b*(1 - w_b)*T_s*Y_a[i] + E_c*(1 - w_c)*T_s*Y_a[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] +
    E_c*(1 - w_c)*T_s*Y_b[i] + E_a*(1 - w_a)*T_s*Y_b[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] +
    E_b*(1 - w_b)*T_s*Y_c[i] + E_a*(1 - w_a)*T_s*Y_c[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] +
    E_c*(1 - w_c)*T_s*Y_ab[i] + E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] + E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] +
    E_b*(1 - w_b)*T_s*Y_ac[i] + E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] + E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] +
    E_a*(1 - w_a)*T_s*Y_bc[i] + E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] + E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] +
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] +
    E_b*(1 - w_b)*T_m*Z_a[i] + E_c*(1 - w_c)*T_m*Z_a[i] +
    E_a*(1 - w_a)*T_m*Z_b[i] + E_c*(1 - w_c)*T_m*Z_b[i] +
    E_b*(1 - w_b)*T_m*Z_c[i] + E_a*(1 - w_a)*T_m*Z_c[i] +
    E_c*(1 - w_c)*T_m*Z_ab[i] + E_b*(1 - w_b)*T_m*Z_ac[i] + E_a*(1 - w_a)*T_m*Z_bc[i] +
    # #add last resort treatments that should go on and off depending on which drugs are being used
    v_c*kappa*T_sr*Y_c[i] + v_ac*kappa*T_sr*Y_ac[i] + v_bc*kappa*T_sr*Y_bc[i] + #added LRT 1/18
    v_abc*kappa*T_sr*Y_abc[i] + 
    (d)*(N[i] - S[i]) 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + sigma*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_s*Y_s[i] - E_b*w_b*T_s*Y_s[i] - E_c*w_c*T_s*Y_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] - (d)*Y_s[i]
  update(Z_s[]) <- Z_s[i] + (1 - sigma)*sum(BI_YZs[,i])*S[i] - E_a*w_a*T_m*Z_s[i] - E_b*w_b*T_m*Z_s[i] - E_c*w_c*T_m*Z_s[i] -
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] - (d)*Z_s[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + sigma*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_a[i] - E_c*(1-w_c)*T_s*Y_a[i] - 
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] - 
    E_b*w_b*T_s*Y_a[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] - E_c*w_c*T_s*Y_a[i] - E_a*L_c_a*w_c*kappa*T_s*Y_a[i] - 
    E_a*kappa*T_sr*Y_a[i] - (d)*Y_a[i] #added last resort treatment only active when using  A (only A scenario)
  
  update(Z_a[]) <- Z_a[i] + (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i] + E_a*w_a*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_a[i] - E_c*(1-w_c)*T_m*Z_a[i] - 
    E_b*w_b*T_m*Z_a[i] - E_c*w_c*T_m*Z_a[i] - (d )*Z_a[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + sigma*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_s*Y_s[i] - E_c*(1 - w_c)*T_s*Y_b[i] - E_a*(1-w_a)*T_s*Y_b[i] - 
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] - 
    E_a*w_a*T_s*Y_b[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] - E_c*w_c*T_s*Y_b[i] - E_b*L_c_b*w_c*kappa*T_s*Y_b[i] - 
    E_b*kappa*T_sr*Y_b[i] - (d)*Y_b[i] #added last resort treatment when E_b is active
  
  update(Z_b[]) <- Z_b[i] + (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i] + E_b*w_b*T_m*Z_s[i] - E_c*(1 - w_c)*T_m*Z_b[i] - E_a*(1-w_a)*T_m*Z_b[i] - 
    E_c*w_c*T_m*Z_b[i] - E_a*w_a*T_m*Z_b[i] - (d)*Z_b[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + sigma*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_s*Y_s[i] - E_b*(1 - w_b)*T_s*Y_c[i] - E_a*(1-w_a)*T_s*Y_c[i] - 
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] - 
    E_a*w_a*T_s*Y_c[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] - E_b*w_b*T_s*Y_c[i] - E_c*L_b_c*w_b*kappa*T_s*Y_c[i] - 
    v_c*kappa*T_sr*Y_c[i] - (d)*Y_c[i] #added / modified last resort treatment 1/18
  
  update(Z_c[]) <- Z_c[i] + (1 - sigma)*f_c*sum(BI_YZc[,i])*S[i] + E_c*w_c*T_m*Z_s[i] - E_b*(1 - w_b)*T_m*Z_c[i] - E_a*(1-w_a)*T_m*Z_c[i] - 
    E_a*w_a*T_m*Z_c[i] - E_b*w_b*T_m*Z_c[i] - (d)*Z_c[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + sigma*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_s*Y_a[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_b[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] -
    E_c*(1 - w_c)*T_s*Y_ab[i] - E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] - E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] - 
    E_c*w_c*T_s*Y_ab[i] - E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] - E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] -
    (d)*Y_ab[i] #removed last resort treatment term 1/18 
  
  update(Z_ab[]) <- Z_ab[i] + (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i] + E_b*w_b*T_m*Z_a[i] + E_a*w_a*T_m*Z_b[i] - E_c*(1 - w_c)*T_m*Z_ab[i] - E_c*w_c*T_m*Z_ab[i] - (d)*Z_ab[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + sigma*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_s*Y_a[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_a[i] + E_a*w_a*T_s*Y_c[i] + E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] -
    E_b*(1 - w_b)*T_s*Y_ac[i] - E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] - E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] - 
    E_b*w_b*T_s*Y_ac[i] - E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] - E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] -
    v_ac*kappa*T_sr*Y_ac[i] - (d)*Y_ac[i] #added / modified last resort treatment 1/18
  
  update(Z_ac[]) <- Z_ac[i] + (1-sigma)*f_ac*sum(BI_YZac[,i])*S[i] + E_c*w_c*T_m*Z_a[i] + E_a*w_a*T_m*Z_c[i] - E_b*(1 - w_b)*T_m*Z_ac[i] - E_b*w_b*T_m*Z_ac[i] - (d )*Z_ac[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + sigma*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_s*Y_b[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_b[i] + E_b*w_b*T_s*Y_c[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_c[i] -
    E_a*(1 - w_a)*T_s*Y_bc[i] - E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] - E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] - 
    E_a*w_a*T_s*Y_bc[i] - E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] - E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] -
    v_bc*kappa*T_sr*Y_bc[i] - (d)*Y_bc[i] #added/modified LRT term 1/18 
  
  update(Z_bc[]) <- Z_bc[i] + (1-sigma)*f_bc*sum(BI_YZbc[,i])*S[i] + E_c*w_c*T_m*Z_b[i] + E_b*w_b*T_m*Z_c[i] - 
    E_a*(1 - w_a)*T_m*Z_bc[i]  - E_a*w_a*T_m*Z_bc[i] - (d)*Z_bc[i] 
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + sigma*f_abc*sum(BI_YZabc[,i])*S[i] + E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] + E_c*w_c*T_s*Y_ab[i] + E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] +
    E_b*w_b*T_s*Y_ac[i] + E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] + E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] + E_a*w_a*T_s*Y_bc[i] + E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] +
    E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] - v_abc*kappa*T_sr*Y_abc[i] - (d)*Y_abc[i] #add LRT prob v_abc 1/18
  
  update(Z_abc[]) <- Z_abc[i] + (1 - sigma)*f_abc*sum(BI_YZabc[,i])*S[i] + E_c*w_c*T_m*Z_ab[i] + E_b*w_b*T_m*Z_ac[i] + E_a*w_a*T_m*Z_bc[i] - (d)*Z_abc[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  denom <- sum(N) - sum(S)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #create binary treatment switches; track when each drug is lost and set locks
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #need A transition that goes from 0.33 down to 0 at first changepoint
  E_a <- 1/3*(lock_threshold_A < 1)
  E_b <- 1/3*(lock_threshold_A < 1) + 0.5*(lock_threshold_A >= 1) - 0.5*(lock_threshold_B >= 1) 
  E_c <- 1/3*(lock_threshold_A < 1) + 0.5*(lock_threshold_A >= 1) + 0.5*(lock_threshold_B >= 1) 
  
  
  #retreatment probabilities
  L_b_a <- 0.5 #will no longer be applicable when A phases out
  L_c_a <- 0.5 #will no longer be applicable when A phases out
  L_c_b <- 0.5 + 0.5*(lock_threshold_A >= 1) #should start at 0.5 and go up to 1 then will not be applicable when B phases out
  L_a_b <- 0.5 - 0.5*(lock_threshold_A >= 1)  #should start at 0.5 then go to 0 when A is phased out
  L_a_c <- 0.5 - 0.5*(lock_threshold_A >= 1)  #should start at 0.5 then go to 0 when A is phased out
  L_b_c <- 0.5 + 0.5*(lock_threshold_A >= 1) - 1*(lock_threshold_B >= 1) #this one starts at 0.5, goes to 1, then needs to go down to 0 when B is phased out
  
  #last resort treatment terms
  v_abc <- 1 #always on, need way out of Y_abc 
  v_bc <- 1*(lock_threshold_A >= 1) #turns on when A goes away so there's a way out of Y_bc
  v_ac <- (lock_threshold_B >= 1) #turns on when B goes away
  v_c <- (lock_threshold_B >= 1) #turns on when B goes away
  
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(denom) <- denom
  output(L_b_a) <- L_b_a
  output(L_b_c) <- L_b_c
  output(L_c_a) <- L_c_a
  output(L_a_b) <- L_a_b
  output(L_a_c) <- L_a_c
  output(L_c_b) <- L_c_b
  output(v_abc) <- v_abc
  output(v_bc) <- v_bc
  output(v_ac) <- v_ac
  output(v_c) <- v_c
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(sum_N) <- sum_N
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  
})


#make discrete sequential model stochastic
sequential_stochastic_5_27 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) + f_a*(Y_a[i]+Z_a[i]) + f_b*(Y_b[i] + Z_b[i]) + f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + f_ac*(Y_ac[i] + Z_ac[i]) + f_bc*(Y_bc[i] + Z_bc[i]) + f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] -sum(BI_S[,i]) * S[i] + #leaving S and going into infectious compartment
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] + #leaving Y_s after successful treatment; this term is in Y_s
    E_b*(1 - w_b)*T_s*Y_a[i] + #Y_a successful treatment with B, matched in Y_a
    E_c*(1 - w_c)*T_s*Y_a[i] + #Y_a successful treatment with C, matched in Y_a
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] + #Y_a successful retreatment with B after A, matched in  Y_a
    E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] + #Y_a successful retreatment with C after A, matched in Y_a
    E_c*(1 - w_c)*T_s*Y_b[i] + #Y_b successful initial treatment with C, matched in Y_b
    E_a*(1 - w_a)*T_s*Y_b[i] + #Y_b successful initial treatment with A, matched in Y_b
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] + #Y_b successful retreatment with A after B, matched in Y_b
    E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] + #Y_b successful retreatment with C after B, matched in Y_b
    E_b*(1 - w_b)*T_s*Y_c[i] + #Y_c successful treatment with B, matched in Y_c
    E_a*(1 - w_a)*T_s*Y_c[i] + #Y_c successful treatment iwth A, matched in Y_c
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] + #Y_c successful retreatment with A after C, matched in Y_c
    E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] + #Y_c successful retreatment with B after C, matched in Y_c
    E_c*(1 - w_c)*T_s*Y_ab[i] + #Y_ab successful treatment with C, matched in Y_ab
    E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] + #Y_ab successful retreatment with C after A, matched in Y_ab
    E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] + #Y_ab successful retreatment with C after B, matched in Y_ab
    E_b*(1 - w_b)*T_s*Y_ac[i] + #Y_ac successful treatment with B, matched in Y_ac
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] + #Y_ac successful retreatment with B after A, matched in Y_ac
    E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] + #Y_ac successful retreatment with B after C, matched in Y_ac
    E_a*(1 - w_a)*T_s*Y_bc[i] + #Y_bc successful treatment with A, matched in Y_bc
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] + #Y_bc successful retreatment with A after B, matched in Y_bc
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] + #Y_bc successful retreatment with A after C, matched in Y_bc
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] + #Successful treatment with A, B, or C from Z_s, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_a[i] + #Z_a successful treatment with B, matched in Z_a
    E_c*(1 - w_c)*T_m*Z_a[i] + #Z_a successful treatment with C, matched in Z_a
    E_a*(1 - w_a)*T_m*Z_b[i] + #Z_b successful treatment with A, matched in Z_b
    E_c*(1 - w_c)*T_m*Z_b[i] + #Z_b successful treatment with C, matched in Z_b
    E_b*(1 - w_b)*T_m*Z_c[i] + #Z_c successful treatment with B, matched in Z_c
    E_a*(1 - w_a)*T_m*Z_c[i] + #Z_c successful treatment with A, matched in Z_c
    E_c*(1 - w_c)*T_m*Z_ab[i] + #Z_ab successful treatment with C, matched in Z_ab
    E_b*(1 - w_b)*T_m*Z_ac[i] + #Z_ac successful treatment with B, matched in Z_ac
    E_a*(1 - w_a)*T_m*Z_bc[i] + #Z_bc successful treatment with A, matched in Z_bc
    # #add last resort treatments that should go on and off depending on which drugs are being used, 
    #added 1/18, adjusted 5/27 in troubleshooting
    v_c*kappa*T_sr*Y_c[i] + #Y_c last resort retreatment, matched in Y_c #CHECK THAT THIS IS DOING WHAT I WANT IT TO DO
    v_ac*kappa*T_sr*Y_ac[i] + #Y_ac last resort retreatment, matched in Y_ac
    v_bc*kappa*T_sr*Y_bc[i] + #Y_bc last resort retreatment, matched in Y_bc
    v_abc*kappa*T_sr*Y_abc[i] + #Y_abc last resort retreatment, matched in Y_abc
    d * Y_s[i] + #Y_s natural recovery, matched in Y_s
    d * Z_s[i] + #Z_s natural recovery, matched in Z_s
    d * Y_a[i] + #Y_a natural recovery, matched in Y_a
    d * Z_a[i] + #Z_a natural recovery, matched in Z_a
    d * Y_b[i] + #Y_b natural recovery, matched in Y_b
    d* Z_b[i] + #Z_b natural recovery, matched in Z_b
    d * Y_c[i] + #Y_c natural recovery, matched in Y_c
    d * Z_c[i] + #Z_c natural recovery, matched in Z_c
    d * Y_ab[i] + #Y_ab natural recovery, matched in Y_ab
    d * Z_ab[i] + #Z_ab natural recovery, matched in Z_ab
    d * Y_ac[i] + #Y_ac natural recovery, matched in Y_ac
    d * Z_ac[i] + #Z_ac natural recovery, matched in Z_ac
    d * Y_bc[i] + #Y_bc natural recovery, matched in Y_bc
    d * Z_bc[i] + #Z_bc natural recovery, matched in Z_bc
    d * Y_abc[i] + #Y_abc natural recovery, matched in Y_abc
    d * Z_abc[i] #Z_abc natural recovery, matched in Z_abc
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + sigma*sum(BI_YZs[,i])*S[i] - #infection 
    E_a*w_a*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with A resulting in resistance, matched in Y_a
    E_b*w_b*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with B resulting in resistance, matched in Y_b
    E_c*w_c*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with C resulting in resistance, matched in Y_c
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] - #successfully treated on the first go, return to S, is matched
    d * Y_s[i] #Y_s natural recovery, matched in S
  
  update(Z_s[]) <- Z_s[i] + (1 - sigma)*sum(BI_YZs[,i])*S[i] - #infection
    E_a*w_a*T_m*Z_s[i] - #Z_s unsuccessful treatment  with A resulting in resistance, matched in Z_a
    E_b*w_b*T_m*Z_s[i] - #Z_s unsuccessful treatment with B resulting in resistance, matched in Z_b
    E_c*w_c*T_m*Z_s[i] - #Z_s unsuccessful treatment with C resulting in resistance, matched in Z_c
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] - #Successful treatment with A, B, or C from Z_s, matched in Z_s
    d * Z_s[i] #Z_s natural recovery, matched in S
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + sigma*f_a*sum(BI_YZa[,i])*S[i] +
    E_a*w_a*T_s*Y_s[i] - #Y_s unsuccessful initial treamtment with A resulting in resistance, matched in Y_s
    E_b*(1 - w_b)*T_s*Y_a[i] - #successful treatment with B, goes back to S and is matched
    E_c*(1 - w_c)*T_s*Y_a[i] - #Y_a successful treatment with C, matched in S
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_a[i] - #Y_a successful retreatment with B after A, matched in S
    E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_a[i] - #Y_a successful retreatment with C after A, matched in S
    E_b*w_b*T_s*Y_a[i] - #Y_a unsuccessful treatment with B resulting in resistance, matched in Y_ab
    E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] - #Y_a unsuccessful retreatment with B after A resulting in resistance, matched in Y_ab
    E_c*w_c*T_s*Y_a[i] - #Y_a unsuccessful treatment with C resulting in resistance, matched in Y_ac
    E_a*L_c_a*w_c*kappa*T_sr*Y_a[i] - #Y_a unsuccessful retreatment with C after A resulting in resistance, matched in Y_ac  #fixed typo in this term on 5/27
    #E_a*kappa*T_sr*Y_a[i] - #added last resort treatment only active when using  A (only A scenario); THIS DOESN'T MAKE SENSE SHOULDN'T BE HERE
    d * Y_a[i] #Y_a natural recovery, matched in S
  
  update(Z_a[]) <- Z_a[i] + (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i] + #infection
    E_a*w_a*T_m*Z_s[i] - #Z_s unsuccessful treatment with A resulting in resistance, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_a[i] -#Z_a successful treatment with B, matched in S
    E_c*(1 - w_c)*T_m*Z_a[i] - #Z_a successful treatment with C, matched in S
    E_b*w_b*T_m*Z_a[i] - #Z_a unsuccessful treatment with B resulting in resistance, matched in Z_ab
    E_c*w_c*T_m*Z_a[i] - #Z_a unsuccessful treatment with C resulting in resistance, matched in Z_ac
    d * Z_a[i] #Z_a natural recovery, matched in S
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + sigma*f_b*sum(BI_YZb[,i])*S[i] +
    E_b*w_b*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with B resulting in resistance, matched in Y_s
    E_c*(1 - w_c)*T_s*Y_b[i] - #Y_b successful initial treatment with C, matched in S
    E_a*(1 - w_a)*T_s*Y_b[i] - #Y_b successful initial treatment with A, matched in S
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_b[i] - #Y_b successful retreatment with A after B, matched in S
    E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_b[i] - #Y_b successful retreatment with C after B, matched in S
    E_a*w_a*T_s*Y_b[i] - #Y_b unsuccessful treatment with A resulting in resistance, matched in Y_ab
    E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] - #Y_b unsuccessful retreatment with A resulting in resistance, matched in Y_ab
    E_c*w_c*T_s*Y_b[i] - #Y_b unsuccessful treatment with C resulting in resistance, matched in Y_bc
    E_b*L_c_b*w_c*kappa*T_sr*Y_b[i] - #Y_b unsuccessful retreatment with C after B resulting in resistance, matched in Y_bc #THIS SHOULD BE T_sr NOT T_s, fixed 5/27
    # E_b*kappa*T_sr*Y_b[i] - #added last resort treatment when E_b is active #5/27 I don't think this should be here
    d * Y_b[i] #Y_b natural recovery, matched in S
  
  update(Z_b[]) <- Z_b[i] + (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i] + #infection 
    E_b*w_b*T_m*Z_s[i] - #Z_s unsuccessful treatment with B resulting in resistance, matched in Z_s
    E_c*(1 - w_c)*T_m*Z_b[i] - #Z_b successful treatment with C, matched in S
    E_a*(1 - w_a)*T_m*Z_b[i] - #Z_b successful treatment with A, matched in S
    E_c*w_c*T_m*Z_b[i] - #Z_b unsuccessful treatment with C, matched in Z_bc
    E_a*w_a*T_m*Z_b[i] - #Z_b unsuccessful treatment with A, matched in Z_ab
    d * Z_b[i] #Z_b natural recovery, matched in S
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + sigma*f_c*sum(BI_YZc[,i])*S[i] + #infection 
    E_c*w_c*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with C resulting in resistance, matched in Y_s
    E_b*(1 - w_b)*T_s*Y_c[i] - #Y_c successful retreatment with B, matched in S
    E_a*(1 - w_a)*T_s*Y_c[i] - #Y_c successful treatment with A, matched in S
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_c[i] - #Y_c successful retreatment iwth A after C, matched in S
    E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_c[i] - #Y_c successful retreatment with B after C, matched in S
    E_a*w_a*T_s*Y_c[i] - #Y_c unsuccessful treatment with A resulting in resistance, matched in Y_ac
    E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] - #Y_c unsuccessful retreatment with A after C resulting in resistance, matched in Y_ac
    E_b*w_b*T_s*Y_c[i] - #Y_c unsuccessful treatment with B resulting in resistance, matched in Y_bc
    E_c*L_b_c*w_b*kappa*T_sr*Y_c[i] - #Y_c unsuccessful retreatment with B after C resulting in resistance, matched in Y_bc #SHOULD BE T_sr NOT T_s, fixed 5/27
    v_c*kappa*T_sr*Y_c[i] - #Y_c last resort retreatment, matched in S
    d * Y_c[i] #Y_c natural recovery, matched in S
  
  update(Z_c[]) <- Z_c[i] + (1 - sigma)*f_c*sum(BI_YZc[,i])*S[i] + #infection
    E_c*w_c*T_m*Z_s[i] - #Z_s unsuccessful treatment with C resulting in resistance, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_c[i] - #Z_c successful treatment with B, matched in S
    E_a*(1 - w_a)*T_m*Z_c[i] - #Z_c successful treatment iwth A, matched in S
    E_a*w_a*T_m*Z_c[i] - #Z_c unsuccessful treatment with A resulting in resistance, matched in Z_ac
    E_b*w_b*T_m*Z_c[i] - #Z_c unsuccessful treatment with B resulting in resistance, matched in Z_bc
    d * Z_c[i] #Z_c natural recovery, matched in S
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + sigma*f_ab*sum(BI_YZab[,i])*S[i] + #infection
    E_b*w_b*T_s*Y_a[i] + #Y_a unsuccessful treatment with B resulting in resistance, matched in Y_a
    E_a*L_b_a*w_b*kappa*T_sr*Y_a[i] + #Y_a unsuccessful retreatment with B after A resulting in resistance, matched in Y_a
    E_a*w_a*T_s*Y_b[i] + #Y_b unsuccessful treatment with A resulting in resistance, matched in Y_a
    E_b*L_a_b*w_a*kappa*T_sr*Y_b[i] - #Y_b unsuccessful retreatment with A resulting in resistance, matched in Y_a
    E_c*(1 - w_c)*T_s*Y_ab[i] - #Y_ab successful treatment with C, matched in S
    E_a*L_c_a*(1 - w_c)*kappa*T_sr*Y_ab[i] - #Y_ab successful retreatment with C after A, matched in S
    E_b*L_c_b*(1 - w_c)*kappa*T_sr*Y_ab[i] - #Y_ab successful retreatment with C after B, matched in S
    E_c*w_c*T_s*Y_ab[i] - #Y_ab unsuccessful treatment with C resulting in resistance, matched in Y_abc
    E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] - #Y_ab unsuccessful retreatment with C after A resulting in resistance, matched in Y_abc
    E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] - #Y_ab unsuccessful retreatment with C after B resulting in resistance, matched in Y_abc
    d * Y_ab[i] #Y_ab natural recovery, matched in S 
  
  update(Z_ab[]) <- Z_ab[i] + (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i] + #infection
    E_b*w_b*T_m*Z_a[i] + #Z_a unsuccessful treatment with B resulting in resistance, matched in Z_a
    E_a*w_a*T_m*Z_b[i] - #Z_b unsuccessful treatment with A, matched in Z_b
    E_c*(1 - w_c)*T_m*Z_ab[i] - #Z_ab successful treatment with C, matched in S
    E_c*w_c*T_m*Z_ab[i] - #Z_ab unsuccessful treatment with C, matched in Z_abc
    d * Z_ab[i] #Z_ab natural recovery, matched in S
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + sigma*f_ac*sum(BI_YZac[,i])*S[i] + #infection
    E_c*w_c*T_s*Y_a[i] + #Y_a unsuccessful treatment with C resulting in resistance, matched in Y_a
    E_a*L_c_a*w_c*kappa*T_sr*Y_a[i] + #Y_a unsuccessful retreatment with C after A resulting in resistance, matched in Y_a
    E_a*w_a*T_s*Y_c[i] + #Y_c unsuccessful treatment with A resulting in resistance, matched in Y_c
    E_c*L_a_c*w_a*kappa*T_sr*Y_c[i] - #Y_c unsuccessful retreatment with A after C resulting in resistance, matched in Y_c
    E_b*(1 - w_b)*T_s*Y_ac[i] - #Y_ac successful treatment with B, matched in S
    E_a*L_b_a*(1 - w_b)*kappa*T_sr*Y_ac[i] - #Y_ac successful retreatment with B after A, matched in S
    E_c*L_b_c*(1 - w_b)*kappa*T_sr*Y_ac[i] - #Y_ac successful retreatment with B after C, matched in S
    E_b*w_b*T_s*Y_ac[i] - #Y_ac unsuccessful treatment with B resulting in resistance, matched in Y_abc
    E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] - #Y_ac unsuccessful retreatment with B after A resulting in resistance, matched in Y_abc
    E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] - #Y_ac unsuccessful retreatment with B after C resulting in resistance, matched in Y_abc
    v_ac*kappa*T_sr*Y_ac[i] - #Y_ac last resort retreatment, matched in S
    d * Y_ac[i] #Y_ac natural recovery, matched in S
  
  update(Z_ac[]) <- Z_ac[i] + (1-sigma)*f_ac*sum(BI_YZac[,i])*S[i] + #infection
    E_c*w_c*T_m*Z_a[i] + #Z_a unsuccessful treatment with C resulting in resistance, matched in Z_a
    E_a*w_a*T_m*Z_c[i] - #Z_c unsuccessful treatment with A resulting in resistance, matched in Z_c
    E_b*(1 - w_b)*T_m*Z_ac[i] - #Z_ac successful treatment with B, matched in S
    E_b*w_b*T_m*Z_ac[i] - #Z_ac unsuccessful treatment with B, matched in Z_abc
    d * Z_ac[i] #Z_ac natural recovery, matched in S
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + sigma*f_bc*sum(BI_YZbc[,i])*S[i] + #infection
    E_c*w_c*T_s*Y_b[i] + #Y_b unsuccessful treatment with C resulting in resistance, matched in Y_b
    E_b*L_c_b*w_c*kappa*T_sr*Y_b[i] + #Y_b unsuccessful retreatment with C after B resulting in resistance, matched in Y_b #THIS SHOULD BE T_sr NOT T_s, fixed 1/27
    E_b*w_b*T_s*Y_c[i] + #Y_c unsuccessful treatment with B resulting in resistance, matched in Y_c
    E_c*L_b_c*w_b*kappa*T_sr*Y_c[i] - #Y_c unsuccessful retreatment with B after C resulting in resistance, matched in Y_c
    E_a*(1 - w_a)*T_s*Y_bc[i] - #Y_bc successful treatment with A, matched in S
    E_b*L_a_b*(1 - w_a)*kappa*T_sr*Y_bc[i] - #Y_bc successful retreatment with A after B, matched in S
    E_c*L_a_c*(1 - w_a)*kappa*T_sr*Y_bc[i] - #Y_bc successful retreatment with A after C, matched in S
    E_a*w_a*T_s*Y_bc[i] - #Y_bc unsuccessful treatment with A resulting in resistance, matched in Y_abc
    E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] - #Y_bc unsuccessful retreatment with A after B resulting in resistance, matched in Y_abc
    E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] - #Y_bc unsuccessful retreatment wtih A after C resulting in resistance, matched in Y_abc
    v_bc*kappa*T_sr*Y_bc[i] - #Y_bc last resort retreatment, matched in S
    d * Y_bc[i] #Y_bc natural recovery, matched in S
  
  update(Z_bc[]) <- Z_bc[i] + (1-sigma)*f_bc*sum(BI_YZbc[,i])*S[i] +
    E_c*w_c*T_m*Z_b[i] + #Z_b unsuccessful treatment with C, matched in Z_b
    E_b*w_b*T_m*Z_c[i] - #Z_c unsuccessful treatment with B resulting in resistance, matched in Z_c
    E_a*(1 - w_a)*T_m*Z_bc[i]  - #Z_bc successful treatment with A, matched in S
    E_a*w_a*T_m*Z_bc[i] - #Z_bc unsuccessful treatment with A resulting in resistance, matched in Z_abc
    d * Z_bc[i] #Z_bc natural recovery, matched in S
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + sigma*f_abc*sum(BI_YZabc[,i])*S[i] + 
    E_a*L_c_a*w_c*kappa*T_sr*Y_ab[i] + #Y_ab unsuccessful retreatment with C after A resulting in resistance, matched in Y_ab
    E_c*w_c*T_s*Y_ab[i] + #Y_ab unsuccessful treatment with C resulting in resistance, matched in Y_ab
    E_b*L_c_b*w_c*kappa*T_sr*Y_ab[i] + #Y_ab unsuccessful retreatment with C after B resulting in resistance, matched in Y_ab
    E_b*w_b*T_s*Y_ac[i] + #Y_ac unsuccessful treatment with B resulting in resistance, matched in Y_ac
    E_a*L_b_a*w_b*kappa*T_sr*Y_ac[i] + #Y_ac unsuccessful retreatment with B after A resulting in resistance, matched in Y_ac
    E_c*L_b_c*w_b*kappa*T_sr*Y_ac[i] + #Y_ac unsuccessful retreatment with B after C resulting in resistance, matched in Y_ac
    E_a*w_a*T_s*Y_bc[i] + #Y_bc unsuccessful treatment with A resulting in resistance, matched in Y_bc
    E_b*L_a_b*w_a*kappa*T_sr*Y_bc[i] + #Y_bc unsuccessful retreatment with A after B resulting in resistance, matched in Y_bc
    E_c*L_a_c*w_a*kappa*T_sr*Y_bc[i] - #Y_bc unsuccessful retreatment wtih A after C resulting in resistance, matched in Y_bc
    v_abc*kappa*T_sr*Y_abc[i] - #Y_abc last resort retreatment, matched in Y_abc
    d * Y_abc[i] #Y_abc natural recovery, matched in S
  
  update(Z_abc[]) <- Z_abc[i] + (1 - sigma)*f_abc*sum(BI_YZabc[,i])*S[i] + 
    E_c*w_c*T_m*Z_ab[i] + #Z_ab unsuccessful treatment with C, matched in Z_ab
    E_b*w_b*T_m*Z_ac[i] + #Z_ac unsuccessful treatment with B, matched in Z_ac
    E_a*w_a*T_m*Z_bc[i] - #Z_bc unsuccessful treatment with A resulting in resistance, matched in Z_bc
    d * Z_abc[i] #Z_abc natural recovery, matched in Z_abc
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)

  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  denom <- sum(N) - sum(S)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  # Binary switches between treatment
  E_a <- (lock_threshold_A < 1)
  E_c <- 1 - (lock_threshold_B < 1)
  E_b <- 1  - E_a - E_c
  
  #Set retreatment terms; for 1-18 sequential model are fixed
  L_b_a <- 0.5
  L_b_c <- 0
  L_c_a <- 0.5
  L_a_b <- 0
  L_a_c <- 0
  L_c_b <- 1
  
  #last resort treatment terms
  v_abc <- 1
  v_bc <- 1 - E_a 
  v_ac <- E_c
  v_c <- E_c
  
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(denom) <- denom
  output(L_b_a) <- L_b_a
  output(L_b_c) <- L_b_c
  output(L_c_a) <- L_c_a
  output(L_a_b) <- L_a_b
  output(L_a_c) <- L_a_c
  output(L_c_b) <- L_c_b 
  output(v_abc) <- v_abc
  output(v_bc) <- v_bc
  output(v_ac) <- v_ac
  output(v_c) <- v_c
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  
  
})


#changing retreatment / last resort treatment mechanisms to simplify, still deterministic, 5/28
sequential_simplified_discrete_5_28 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] -sum(BI_S[,i]) * S[i] + #leaving S and going into infectious compartment
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] + #leaving Y_s after successful treatment; this term is in Y_s
    E_b*(1 - w_b)*T_s*Y_a[i] + #Y_a successful treatment with B, matched in Y_a
    E_c*(1 - w_c)*T_s*Y_a[i] + #Y_a successful treatment with C, matched in Y_a
    E_a*kappa*T_sr*Y_a[i] + #Y_a successful retreatment, matched in  Y_a
    E_c*(1 - w_c)*T_s*Y_b[i] + #Y_b successful initial treatment with C, matched in Y_b
    E_a*(1 - w_a)*T_s*Y_b[i] + #Y_b successful initial treatment with A, matched in Y_b
    E_b*kappa*T_sr*Y_b[i] + #Y_b successful retreatment, matched in Y_b
    E_b*(1 - w_b)*T_s*Y_c[i] + #Y_c successful treatment with B, matched in Y_c
    E_a*(1 - w_a)*T_s*Y_c[i] + #Y_c successful treatment iwth A, matched in Y_c
    E_c*kappa*T_sr*Y_c[i] + #Y_c successful retreatment with A after C, matched in Y_c
    E_c*(1 - w_c)*T_s*Y_ab[i] + #Y_ab successful treatment with C, matched in Y_ab
    (1 - E_c)*kappa*T_sr*Y_ab[i] + #Y_ab successful retreatment, matched in Y_ab
    E_b*(1 - w_b)*T_s*Y_ac[i] + #Y_ac successful treatment with B, matched in Y_ac
    (1 - E_b)*kappa*T_sr*Y_ac[i] + #Y_ac successful retreatment matched in Y_ac
    E_a*(1 - w_a)*T_s*Y_bc[i] + #Y_bc successful treatment with A, matched in Y_bc
    (1 - E_a)*kappa*T_sr*Y_bc[i] + #Y_bc successful retreatment, matched in Y_bc
    kappa*T_sr*Y_abc[i] + #Y_abc last resort retreatment, matched in Y_abc
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] + #Successful treatment with A, B, or C from Z_s, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_a[i] + #Z_a successful treatment with B, matched in Z_a
    E_c*(1 - w_c)*T_m*Z_a[i] + #Z_a successful treatment with C, matched in Z_a
    E_a*(1 - w_a)*T_m*Z_b[i] + #Z_b successful treatment with A, matched in Z_b
    E_c*(1 - w_c)*T_m*Z_b[i] + #Z_b successful treatment with C, matched in Z_b
    E_b*(1 - w_b)*T_m*Z_c[i] + #Z_c successful treatment with B, matched in Z_c
    E_a*(1 - w_a)*T_m*Z_c[i] + #Z_c successful treatment with A, matched in Z_c
    E_c*(1 - w_c)*T_m*Z_ab[i] + #Z_ab successful treatment with C, matched in Z_ab
    E_b*(1 - w_b)*T_m*Z_ac[i] + #Z_ac successful treatment with B, matched in Z_ac
    E_a*(1 - w_a)*T_m*Z_bc[i] + #Z_bc successful treatment with A, matched in Z_bc
    d * Y_s[i] + #Y_s natural recovery, matched in Y_s
    d * Z_s[i] + #Z_s natural recovery, matched in Z_s
    d * Y_a[i] + #Y_a natural recovery, matched in Y_a
    d * Z_a[i] + #Z_a natural recovery, matched in Z_a
    d * Y_b[i] + #Y_b natural recovery, matched in Y_b
    d * Z_b[i] + #Z_b natural recovery, matched in Z_b
    d * Y_c[i] + #Y_c natural recovery, matched in Y_c
    d * Z_c[i] + #Z_c natural recovery, matched in Z_c
    d * Y_ab[i] + #Y_ab natural recovery, matched in Y_ab
    d * Z_ab[i] + #Z_ab natural recovery, matched in Z_ab
    d * Y_ac[i] + #Y_ac natural recovery, matched in Y_ac
    d * Z_ac[i] + #Z_ac natural recovery, matched in Z_ac
    d * Y_bc[i] + #Y_bc natural recovery, matched in Y_bc
    d * Z_bc[i] + #Z_bc natural recovery, matched in Z_bc
    d * Y_abc[i] + #Y_abc natural recovery, matched in Y_abc
    d * Z_abc[i] #Z_abc natural recovery, matched in Z_abc
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + sigma*sum(BI_YZs[,i])*S[i] - #infection 
    E_a*w_a*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with A resulting in resistance, matched in Y_a
    E_b*w_b*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with B resulting in resistance, matched in Y_b
    E_c*w_c*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with C resulting in resistance, matched in Y_c
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] - #successfully treated on the first go, return to S, is matched
    d * Y_s[i] #Y_s natural recovery, matched in S
  
  update(Z_s[]) <- Z_s[i] + (1 - sigma)*sum(BI_YZs[,i])*S[i] - #infection
    E_a*w_a*T_m*Z_s[i] - #Z_s unsuccessful treatment  with A resulting in resistance, matched in Z_a
    E_b*w_b*T_m*Z_s[i] - #Z_s unsuccessful treatment with B resulting in resistance, matched in Z_b
    E_c*w_c*T_m*Z_s[i] - #Z_s unsuccessful treatment with C resulting in resistance, matched in Z_c
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] - #Successful treatment with A, B, or C from Z_s, matched in Z_s
    d * Z_s[i] #Z_s natural recovery, matched in S
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + sigma*sum(BI_YZa[,i])*S[i] +
    E_a*w_a*T_s*Y_s[i] - #Y_s unsuccessful initial treamtment with A resulting in resistance, matched in Y_s
    E_b*(1 - w_b)*T_s*Y_a[i] - #successful treatment with B, goes back to S and is matched
    E_c*(1 - w_c)*T_s*Y_a[i] - #Y_a successful treatment with C, matched in S
    E_b*w_b*T_s*Y_a[i] - #Y_a unsuccessful treatment with B resulting in resistance, matched in Y_ab
    E_c*w_c*T_s*Y_a[i] - #Y_a unsuccessful treatment with C resulting in resistance, matched in Y_ac
    E_a*kappa*T_sr*Y_a[i] - #Y_a retreatment with last resort (always successful by default), matched in S
    d * Y_a[i] #Y_a natural recovery, matched in S
  
  update(Z_a[]) <- Z_a[i] + (1 - sigma)*sum(BI_YZa[,i])*S[i] + #infection
    E_a*w_a*T_m*Z_s[i] - #Z_s unsuccessful treatment with A resulting in resistance, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_a[i] -#Z_a successful treatment with B, matched in S
    E_c*(1 - w_c)*T_m*Z_a[i] - #Z_a successful treatment with C, matched in S
    E_b*w_b*T_m*Z_a[i] - #Z_a unsuccessful treatment with B resulting in resistance, matched in Z_ab
    E_c*w_c*T_m*Z_a[i] - #Z_a unsuccessful treatment with C resulting in resistance, matched in Z_ac
    d * Z_a[i] #Z_a natural recovery, matched in S
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + sigma*sum(BI_YZb[,i])*S[i] +
    E_b*w_b*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with B resulting in resistance, matched in Y_s
    E_c*(1 - w_c)*T_s*Y_b[i] - #Y_b successful initial treatment with C, matched in S
    E_a*(1 - w_a)*T_s*Y_b[i] - #Y_b successful initial treatment with A, matched in S
    E_b*kappa*T_sr*Y_b[i] - #Y_b retreatment with last resort (always successful by default), matched in S
    E_a*w_a*T_s*Y_b[i] - #Y_b unsuccessful treatment with A resulting in resistance, matched in Y_ab
    E_c*w_c*T_s*Y_b[i] - #Y_b unsuccessful treatment with C resulting in resistance, matched in Y_bc
    d * Y_b[i] #Y_b natural recovery, matched in S
  
  update(Z_b[]) <- Z_b[i] + (1 - sigma)*sum(BI_YZb[,i])*S[i] + #infection 
    E_b*w_b*T_m*Z_s[i] - #Z_s unsuccessful treatment with B resulting in resistance, matched in Z_s
    E_c*(1 - w_c)*T_m*Z_b[i] - #Z_b successful treatment with C, matched in S
    E_a*(1 - w_a)*T_m*Z_b[i] - #Z_b successful treatment with A, matched in S
    E_c*w_c*T_m*Z_b[i] - #Z_b unsuccessful treatment with C, matched in Z_bc
    E_a*w_a*T_m*Z_b[i] - #Z_b unsuccessful treatment with A, matched in Z_ab
    d * Z_b[i] #Z_b natural recovery, matched in S
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + sigma*sum(BI_YZc[,i])*S[i] + #infection 
    E_c*w_c*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with C resulting in resistance, matched in Y_s
    E_b*(1 - w_b)*T_s*Y_c[i] - #Y_c successful retreatment with B, matched in S
    E_a*(1 - w_a)*T_s*Y_c[i] - #Y_c successful treatment with A, matched in S
    E_c*kappa*T_sr*Y_c[i] - #Y_c retreatment with last resort (always successful by default), matched in S
    E_a*w_a*T_s*Y_c[i] - #Y_c unsuccessful treatment with A resulting in resistance, matched in Y_ac
    E_b*w_b*T_s*Y_c[i] - #Y_c unsuccessful treatment with B resulting in resistance, matched in Y_bc
    d * Y_c[i] #Y_c natural recovery, matched in S
  
  update(Z_c[]) <- Z_c[i] + (1 - sigma)*sum(BI_YZc[,i])*S[i] + #infection
    E_c*w_c*T_m*Z_s[i] - #Z_s unsuccessful treatment with C resulting in resistance, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_c[i] - #Z_c successful treatment with B, matched in S
    E_a*(1 - w_a)*T_m*Z_c[i] - #Z_c successful treatment iwth A, matched in S
    E_a*w_a*T_m*Z_c[i] - #Z_c unsuccessful treatment with A resulting in resistance, matched in Z_ac
    E_b*w_b*T_m*Z_c[i] - #Z_c unsuccessful treatment with B resulting in resistance, matched in Z_bc
    d * Z_c[i] #Z_c natural recovery, matched in S
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + sigma*sum(BI_YZab[,i])*S[i] + #infection
    E_b*w_b*T_s*Y_a[i] + #Y_a unsuccessful treatment with B resulting in resistance, matched in Y_a
    E_a*w_a*T_s*Y_b[i] - #Y_b unsuccessful treatment with A resulting in resistance, matched in Y_a
    E_c*(1 - w_c)*T_s*Y_ab[i] - #Y_ab successful treatment with C, matched in S
    (1 - E_c)*kappa*T_sr*Y_ab[i] - #Y_ab retreatment with last resort (always successful by default), matched in S
    E_c*w_c*T_s*Y_ab[i] - #Y_ab unsuccessful treatment with C resulting in resistance, matched in Y_abc
    d * Y_ab[i] #Y_ab natural recovery, matched in S 
  
  update(Z_ab[]) <- Z_ab[i] + (1-sigma)*sum(BI_YZab[,i])*S[i] + #infection
    E_b*w_b*T_m*Z_a[i] + #Z_a unsuccessful treatment with B resulting in resistance, matched in Z_a
    E_a*w_a*T_m*Z_b[i] - #Z_b unsuccessful treatment with A, matched in Z_b
    E_c*(1 - w_c)*T_m*Z_ab[i] - #Z_ab successful treatment with C, matched in S
    E_c*w_c*T_m*Z_ab[i] - #Z_ab unsuccessful treatment with C, matched in Z_abc
    d * Z_ab[i] #Z_ab natural recovery, matched in S
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + sigma*sum(BI_YZac[,i])*S[i] + #infection
    E_c*w_c*T_s*Y_a[i] + #Y_a unsuccessful treatment with C resulting in resistance, matched in Y_a
    E_a*w_a*T_s*Y_c[i] - #Y_c unsuccessful treatment with A resulting in resistance, matched in Y_c
    E_b*(1 - w_b)*T_s*Y_ac[i] - #Y_ac successful treatment with B, matched in S
    (1 - E_b)*kappa*T_sr*Y_ac[i] - #Y_ac retreatment with last resort (always successful by default), matched in S
    E_b*w_b*T_s*Y_ac[i] - #Y_ac unsuccessful treatment with B resulting in resistance, matched in Y_abc
    d * Y_ac[i] #Y_ac natural recovery, matched in S
  
  update(Z_ac[]) <- Z_ac[i] + (1-sigma)*sum(BI_YZac[,i])*S[i] + #infection
    E_c*w_c*T_m*Z_a[i] + #Z_a unsuccessful treatment with C resulting in resistance, matched in Z_a
    E_a*w_a*T_m*Z_c[i] - #Z_c unsuccessful treatment with A resulting in resistance, matched in Z_c
    E_b*(1 - w_b)*T_m*Z_ac[i] - #Z_ac successful treatment with B, matched in S
    E_b*w_b*T_m*Z_ac[i] - #Z_ac unsuccessful treatment with B, matched in Z_abc
    d * Z_ac[i] #Z_ac natural recovery, matched in S
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + sigma*sum(BI_YZbc[,i])*S[i] + #infection
    E_c*w_c*T_s*Y_b[i] + #Y_b unsuccessful treatment with C resulting in resistance, matched in Y_b
    E_b*w_b*T_s*Y_c[i] - #Y_c unsuccessful treatment with B resulting in resistance, matched in Y_c
    E_a*(1 - w_a)*T_s*Y_bc[i] - #Y_bc successful treatment with A, matched in S
    (1 - E_a)*kappa*T_sr*Y_bc[i] - #Y_bc retreatment with last resort (always successful by default), matched in S
    E_a*w_a*T_s*Y_bc[i] - #Y_bc unsuccessful treatment with A resulting in resistance, matched in Y_abc
    d * Y_bc[i] #Y_bc natural recovery, matched in S
  
  update(Z_bc[]) <- Z_bc[i] + (1-sigma)*sum(BI_YZbc[,i])*S[i] +
    E_c*w_c*T_m*Z_b[i] + #Z_b unsuccessful treatment with C, matched in Z_b
    E_b*w_b*T_m*Z_c[i] - #Z_c unsuccessful treatment with B resulting in resistance, matched in Z_c
    E_a*(1 - w_a)*T_m*Z_bc[i]  - #Z_bc successful treatment with A, matched in S
    E_a*w_a*T_m*Z_bc[i] - #Z_bc unsuccessful treatment with A resulting in resistance, matched in Z_abc
    d * Z_bc[i] #Z_bc natural recovery, matched in S
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + sigma*sum(BI_YZabc[,i])*S[i] + 
    E_c*w_c*T_s*Y_ab[i] + #Y_ab unsuccessful treatment with C resulting in resistance, matched in Y_ab
    E_b*w_b*T_s*Y_ac[i] + #Y_ac unsuccessful treatment with B resulting in resistance, matched in Y_ac
    E_a*w_a*T_s*Y_bc[i] - #Y_bc unsuccessful treatment with A resulting in resistance, matched in Y_bc
    kappa*T_sr*Y_abc[i] - #Y_abc last resort retreatment (successful by default), matched in S
    d * Y_abc[i] #Y_abc natural recovery, matched in S
  
  update(Z_abc[]) <- Z_abc[i] + (1 - sigma)*sum(BI_YZabc[,i])*S[i] + 
    E_c*w_c*T_m*Z_ab[i] + #Z_ab unsuccessful treatment with C, matched in Z_ab
    E_b*w_b*T_m*Z_ac[i] + #Z_ac unsuccessful treatment with B, matched in Z_ac
    E_a*w_a*T_m*Z_bc[i] - #Z_bc unsuccessful treatment with A resulting in resistance, matched in Z_bc
    d * Z_abc[i] #Z_abc natural recovery, matched in Z_abc
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #Binary switches between treatment
  E_a <- (lock_threshold_A < 1)
  E_c <- 1 - (lock_threshold_B < 1)
  E_b <- 1  - E_a - E_c
  
  sum_tx_probs = E_a + E_c + E_b 
  
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev
  output(sum_tx_probs) <- sum_tx_probs
  
  
})

#changing retreatment / last resort treamtnet mechanisms to simplify for EA model, still deterministic, 5/29
equal_allocation_simplified_discrete_5_28 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] -sum(BI_S[,i]) * S[i] + #leaving S and going into infectious compartment
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] + #leaving Y_s after successful treatment; this term is in Y_s
    E_b*(1 - w_b)*T_s*Y_a[i] + #Y_a successful treatment with B, matched in Y_a
    E_c*(1 - w_c)*T_s*Y_a[i] + #Y_a successful treatment with C, matched in Y_a
    E_a*kappa*T_sr*Y_a[i] + #Y_a successful retreatment, matched in  Y_a
    E_c*(1 - w_c)*T_s*Y_b[i] + #Y_b successful initial treatment with C, matched in Y_b
    E_a*(1 - w_a)*T_s*Y_b[i] + #Y_b successful initial treatment with A, matched in Y_b
    E_b*kappa*T_sr*Y_b[i] + #Y_b successful retreatment, matched in Y_b
    E_b*(1 - w_b)*T_s*Y_c[i] + #Y_c successful treatment with B, matched in Y_c
    E_a*(1 - w_a)*T_s*Y_c[i] + #Y_c successful treatment iwth A, matched in Y_c
    E_c*kappa*T_sr*Y_c[i] + #Y_c successful retreatment with A after C, matched in Y_c
    E_c*(1 - w_c)*T_s*Y_ab[i] + #Y_ab successful treatment with C, matched in Y_ab
    (1 - E_c)*kappa*T_sr*Y_ab[i] + #Y_ab successful retreatment, matched in Y_ab
    E_b*(1 - w_b)*T_s*Y_ac[i] + #Y_ac successful treatment with B, matched in Y_ac
    (1 - E_b)*kappa*T_sr*Y_ac[i] + #Y_ac successful retreatment matched in Y_ac
    E_a*(1 - w_a)*T_s*Y_bc[i] + #Y_bc successful treatment with A, matched in Y_bc
    (1 - E_a)*kappa*T_sr*Y_bc[i] + #Y_bc successful retreatment, matched in Y_bc
    kappa*T_sr*Y_abc[i] + #Y_abc last resort retreatment, matched in Y_abc
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] + #Successful treatment with A, B, or C from Z_s, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_a[i] + #Z_a successful treatment with B, matched in Z_a
    E_c*(1 - w_c)*T_m*Z_a[i] + #Z_a successful treatment with C, matched in Z_a
    E_a*(1 - w_a)*T_m*Z_b[i] + #Z_b successful treatment with A, matched in Z_b
    E_c*(1 - w_c)*T_m*Z_b[i] + #Z_b successful treatment with C, matched in Z_b
    E_b*(1 - w_b)*T_m*Z_c[i] + #Z_c successful treatment with B, matched in Z_c
    E_a*(1 - w_a)*T_m*Z_c[i] + #Z_c successful treatment with A, matched in Z_c
    E_c*(1 - w_c)*T_m*Z_ab[i] + #Z_ab successful treatment with C, matched in Z_ab
    E_b*(1 - w_b)*T_m*Z_ac[i] + #Z_ac successful treatment with B, matched in Z_ac
    E_a*(1 - w_a)*T_m*Z_bc[i] + #Z_bc successful treatment with A, matched in Z_bc
    d * Y_s[i] + #Y_s natural recovery, matched in Y_s
    d * Z_s[i] + #Z_s natural recovery, matched in Z_s
    d * Y_a[i] + #Y_a natural recovery, matched in Y_a
    d * Z_a[i] + #Z_a natural recovery, matched in Z_a
    d * Y_b[i] + #Y_b natural recovery, matched in Y_b
    d * Z_b[i] + #Z_b natural recovery, matched in Z_b
    d * Y_c[i] + #Y_c natural recovery, matched in Y_c
    d * Z_c[i] + #Z_c natural recovery, matched in Z_c
    d * Y_ab[i] + #Y_ab natural recovery, matched in Y_ab
    d * Z_ab[i] + #Z_ab natural recovery, matched in Z_ab
    d * Y_ac[i] + #Y_ac natural recovery, matched in Y_ac
    d * Z_ac[i] + #Z_ac natural recovery, matched in Z_ac
    d * Y_bc[i] + #Y_bc natural recovery, matched in Y_bc
    d * Z_bc[i] + #Z_bc natural recovery, matched in Z_bc
    d * Y_abc[i] + #Y_abc natural recovery, matched in Y_abc
    d * Z_abc[i] #Z_abc natural recovery, matched in Z_abc
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + sigma*sum(BI_YZs[,i])*S[i] - #infection 
    E_a*w_a*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with A resulting in resistance, matched in Y_a
    E_b*w_b*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with B resulting in resistance, matched in Y_b
    E_c*w_c*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with C resulting in resistance, matched in Y_c
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_s*Y_s[i] - #successfully treated on the first go, return to S, is matched
    d * Y_s[i] #Y_s natural recovery, matched in S
  
  update(Z_s[]) <- Z_s[i] + (1 - sigma)*sum(BI_YZs[,i])*S[i] - #infection
    E_a*w_a*T_m*Z_s[i] - #Z_s unsuccessful treatment  with A resulting in resistance, matched in Z_a
    E_b*w_b*T_m*Z_s[i] - #Z_s unsuccessful treatment with B resulting in resistance, matched in Z_b
    E_c*w_c*T_m*Z_s[i] - #Z_s unsuccessful treatment with C resulting in resistance, matched in Z_c
    (1 - E_a*w_a - E_b*w_b - E_c*w_c)*T_m*Z_s[i] - #Successful treatment with A, B, or C from Z_s, matched in Z_s
    d * Z_s[i] #Z_s natural recovery, matched in S
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + sigma*sum(BI_YZa[,i])*S[i] +
    E_a*w_a*T_s*Y_s[i] - #Y_s unsuccessful initial treamtment with A resulting in resistance, matched in Y_s
    E_b*(1 - w_b)*T_s*Y_a[i] - #successful treatment with B, goes back to S and is matched
    E_c*(1 - w_c)*T_s*Y_a[i] - #Y_a successful treatment with C, matched in S
    E_b*w_b*T_s*Y_a[i] - #Y_a unsuccessful treatment with B resulting in resistance, matched in Y_ab
    E_c*w_c*T_s*Y_a[i] - #Y_a unsuccessful treatment with C resulting in resistance, matched in Y_ac
    E_a*kappa*T_sr*Y_a[i] - #Y_a retreatment with last resort (always successful by default), matched in S
    d * Y_a[i] #Y_a natural recovery, matched in S
  
  update(Z_a[]) <- Z_a[i] + (1 - sigma)*sum(BI_YZa[,i])*S[i] + #infection
    E_a*w_a*T_m*Z_s[i] - #Z_s unsuccessful treatment with A resulting in resistance, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_a[i] -#Z_a successful treatment with B, matched in S
    E_c*(1 - w_c)*T_m*Z_a[i] - #Z_a successful treatment with C, matched in S
    E_b*w_b*T_m*Z_a[i] - #Z_a unsuccessful treatment with B resulting in resistance, matched in Z_ab
    E_c*w_c*T_m*Z_a[i] - #Z_a unsuccessful treatment with C resulting in resistance, matched in Z_ac
    d * Z_a[i] #Z_a natural recovery, matched in S
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + sigma*sum(BI_YZb[,i])*S[i] +
    E_b*w_b*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with B resulting in resistance, matched in Y_s
    E_c*(1 - w_c)*T_s*Y_b[i] - #Y_b successful initial treatment with C, matched in S
    E_a*(1 - w_a)*T_s*Y_b[i] - #Y_b successful initial treatment with A, matched in S
    E_b*kappa*T_sr*Y_b[i] - #Y_b retreatment with last resort (always successful by default), matched in S
    E_a*w_a*T_s*Y_b[i] - #Y_b unsuccessful treatment with A resulting in resistance, matched in Y_ab
    E_c*w_c*T_s*Y_b[i] - #Y_b unsuccessful treatment with C resulting in resistance, matched in Y_bc
    d * Y_b[i] #Y_b natural recovery, matched in S
  
  update(Z_b[]) <- Z_b[i] + (1 - sigma)*sum(BI_YZb[,i])*S[i] + #infection 
    E_b*w_b*T_m*Z_s[i] - #Z_s unsuccessful treatment with B resulting in resistance, matched in Z_s
    E_c*(1 - w_c)*T_m*Z_b[i] - #Z_b successful treatment with C, matched in S
    E_a*(1 - w_a)*T_m*Z_b[i] - #Z_b successful treatment with A, matched in S
    E_c*w_c*T_m*Z_b[i] - #Z_b unsuccessful treatment with C, matched in Z_bc
    E_a*w_a*T_m*Z_b[i] - #Z_b unsuccessful treatment with A, matched in Z_ab
    d * Z_b[i] #Z_b natural recovery, matched in S
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + sigma*sum(BI_YZc[,i])*S[i] + #infection 
    E_c*w_c*T_s*Y_s[i] - #Y_s unsuccessful initial treatment with C resulting in resistance, matched in Y_s
    E_b*(1 - w_b)*T_s*Y_c[i] - #Y_c successful retreatment with B, matched in S
    E_a*(1 - w_a)*T_s*Y_c[i] - #Y_c successful treatment with A, matched in S
    E_c*kappa*T_sr*Y_c[i] - #Y_c retreatment with last resort (always successful by default), matched in S
    E_a*w_a*T_s*Y_c[i] - #Y_c unsuccessful treatment with A resulting in resistance, matched in Y_ac
    E_b*w_b*T_s*Y_c[i] - #Y_c unsuccessful treatment with B resulting in resistance, matched in Y_bc
    d * Y_c[i] #Y_c natural recovery, matched in S
  
  update(Z_c[]) <- Z_c[i] + (1 - sigma)*sum(BI_YZc[,i])*S[i] + #infection
    E_c*w_c*T_m*Z_s[i] - #Z_s unsuccessful treatment with C resulting in resistance, matched in Z_s
    E_b*(1 - w_b)*T_m*Z_c[i] - #Z_c successful treatment with B, matched in S
    E_a*(1 - w_a)*T_m*Z_c[i] - #Z_c successful treatment iwth A, matched in S
    E_a*w_a*T_m*Z_c[i] - #Z_c unsuccessful treatment with A resulting in resistance, matched in Z_ac
    E_b*w_b*T_m*Z_c[i] - #Z_c unsuccessful treatment with B resulting in resistance, matched in Z_bc
    d * Z_c[i] #Z_c natural recovery, matched in S
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + sigma*sum(BI_YZab[,i])*S[i] + #infection
    E_b*w_b*T_s*Y_a[i] + #Y_a unsuccessful treatment with B resulting in resistance, matched in Y_a
    E_a*w_a*T_s*Y_b[i] - #Y_b unsuccessful treatment with A resulting in resistance, matched in Y_a
    E_c*(1 - w_c)*T_s*Y_ab[i] - #Y_ab successful treatment with C, matched in S
    (1 - E_c)*kappa*T_sr*Y_ab[i] - #Y_ab retreatment with last resort (always successful by default), matched in S
    E_c*w_c*T_s*Y_ab[i] - #Y_ab unsuccessful treatment with C resulting in resistance, matched in Y_abc
    d * Y_ab[i] #Y_ab natural recovery, matched in S 
  
  update(Z_ab[]) <- Z_ab[i] + (1-sigma)*sum(BI_YZab[,i])*S[i] + #infection
    E_b*w_b*T_m*Z_a[i] + #Z_a unsuccessful treatment with B resulting in resistance, matched in Z_a
    E_a*w_a*T_m*Z_b[i] - #Z_b unsuccessful treatment with A, matched in Z_b
    E_c*(1 - w_c)*T_m*Z_ab[i] - #Z_ab successful treatment with C, matched in S
    E_c*w_c*T_m*Z_ab[i] - #Z_ab unsuccessful treatment with C, matched in Z_abc
    d * Z_ab[i] #Z_ab natural recovery, matched in S
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + sigma*sum(BI_YZac[,i])*S[i] + #infection
    E_c*w_c*T_s*Y_a[i] + #Y_a unsuccessful treatment with C resulting in resistance, matched in Y_a
    E_a*w_a*T_s*Y_c[i] - #Y_c unsuccessful treatment with A resulting in resistance, matched in Y_c
    E_b*(1 - w_b)*T_s*Y_ac[i] - #Y_ac successful treatment with B, matched in S
    (1 - E_b)*kappa*T_sr*Y_ac[i] - #Y_ac retreatment with last resort (always successful by default), matched in S
    E_b*w_b*T_s*Y_ac[i] - #Y_ac unsuccessful treatment with B resulting in resistance, matched in Y_abc
    d * Y_ac[i] #Y_ac natural recovery, matched in S
  
  update(Z_ac[]) <- Z_ac[i] + (1-sigma)*sum(BI_YZac[,i])*S[i] + #infection
    E_c*w_c*T_m*Z_a[i] + #Z_a unsuccessful treatment with C resulting in resistance, matched in Z_a
    E_a*w_a*T_m*Z_c[i] - #Z_c unsuccessful treatment with A resulting in resistance, matched in Z_c
    E_b*(1 - w_b)*T_m*Z_ac[i] - #Z_ac successful treatment with B, matched in S
    E_b*w_b*T_m*Z_ac[i] - #Z_ac unsuccessful treatment with B, matched in Z_abc
    d * Z_ac[i] #Z_ac natural recovery, matched in S
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + sigma*sum(BI_YZbc[,i])*S[i] + #infection
    E_c*w_c*T_s*Y_b[i] + #Y_b unsuccessful treatment with C resulting in resistance, matched in Y_b
    E_b*w_b*T_s*Y_c[i] - #Y_c unsuccessful treatment with B resulting in resistance, matched in Y_c
    E_a*(1 - w_a)*T_s*Y_bc[i] - #Y_bc successful treatment with A, matched in S
    (1 - E_a)*kappa*T_sr*Y_bc[i] - #Y_bc retreatment with last resort (always successful by default), matched in S
    E_a*w_a*T_s*Y_bc[i] - #Y_bc unsuccessful treatment with A resulting in resistance, matched in Y_abc
    d * Y_bc[i] #Y_bc natural recovery, matched in S
  
  update(Z_bc[]) <- Z_bc[i] + (1-sigma)*sum(BI_YZbc[,i])*S[i] +
    E_c*w_c*T_m*Z_b[i] + #Z_b unsuccessful treatment with C, matched in Z_b
    E_b*w_b*T_m*Z_c[i] - #Z_c unsuccessful treatment with B resulting in resistance, matched in Z_c
    E_a*(1 - w_a)*T_m*Z_bc[i]  - #Z_bc successful treatment with A, matched in S
    E_a*w_a*T_m*Z_bc[i] - #Z_bc unsuccessful treatment with A resulting in resistance, matched in Z_abc
    d * Z_bc[i] #Z_bc natural recovery, matched in S
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + sigma*sum(BI_YZabc[,i])*S[i] + 
    E_c*w_c*T_s*Y_ab[i] + #Y_ab unsuccessful treatment with C resulting in resistance, matched in Y_ab
    E_b*w_b*T_s*Y_ac[i] + #Y_ac unsuccessful treatment with B resulting in resistance, matched in Y_ac
    E_a*w_a*T_s*Y_bc[i] - #Y_bc unsuccessful treatment with A resulting in resistance, matched in Y_bc
    kappa*T_sr*Y_abc[i] - #Y_abc last resort retreatment (successful by default), matched in S
    d * Y_abc[i] #Y_abc natural recovery, matched in S
  
  update(Z_abc[]) <- Z_abc[i] + (1 - sigma)*sum(BI_YZabc[,i])*S[i] + 
    E_c*w_c*T_m*Z_ab[i] + #Z_ab unsuccessful treatment with C, matched in Z_ab
    E_b*w_b*T_m*Z_ac[i] + #Z_ac unsuccessful treatment with B, matched in Z_ac
    E_a*w_a*T_m*Z_bc[i] - #Z_bc unsuccessful treatment with A resulting in resistance, matched in Z_bc
    d * Z_abc[i] #Z_abc natural recovery, matched in Z_abc
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #make indicator variables for which ordering we're in 
  all_off <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_A <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_B <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  only_C <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1)
  A_B_on <- (lock_threshold_A >= 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  update(C_remain_lock) <- C_remain_lock + A_B_on
  A_C_on <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1) 
  update(B_remain_lock) <- B_remain_lock + A_C_on
  B_C_on <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C >= 1) 
  update(A_remain_lock) <- A_remain_lock + B_C_on
  
  E_a <- 1/3*all_off + 0.5*(only_B + only_C) + max((B_C_on), (A_remain_lock >= 1))
  E_b <- 1/3*all_off + 0.5*(only_A + only_C) + max((A_C_on),(B_remain_lock >= 1))
  E_c <- 1/3*all_off + 0.5*(only_A + only_B) + max((A_B_on),(C_remain_lock >= 1))
  
  sum_tx_probs = E_a + E_c + E_b 
  
  
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  initial(A_remain_lock) <- 0
  initial(B_remain_lock) <- 0
  initial(C_remain_lock) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev
  output(sum_tx_probs) <- sum_tx_probs
  
  
})

#make the models stochastic, 5/29
sequential_stochastic_6_2 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] - n_SOut[i] + n_YsS[i] + n_YaS[i] + n_YbS[i] + n_YcS[i] + n_YabS[i] + n_YacS[i] + n_YbcS[i] +
    n_YabcOut[i] + n_ZsS[i] + n_ZaS[i] + n_ZbS[i] + n_ZcS[i] + n_ZabS[i] + n_ZacS[i] + n_ZbcS[i] + n_ZabcOut[i] 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + n_SYs[i] - n_YsOut[i]
  update(Z_s[]) <- Z_s[i] + n_SZs[i] - n_ZsOut[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + n_SYa[i] + n_YsYa[i] - n_YaOut[i]
  update(Z_a[]) <- Z_a[i] + n_SZa[i] + n_ZsZa[i] - n_ZaOut[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + n_SYb[i] + n_YsYb[i] - n_YbOut[i]
  update(Z_b[]) <- Z_b[i] + n_SZb[i] + n_ZsZb[i] - n_ZbOut[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + n_SYc[i] + n_YsYc[i] - n_YcOut[i] 
  update(Z_c[]) <- Z_c[i] + n_SZc[i] + n_ZsZc[i] - n_ZcOut[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + n_SYab[i] + n_YaYab[i] + n_YbYab[i] - n_YabOut[i]
  update(Z_ab[]) <- Z_ab[i] + n_SZab[i] + n_ZaZab[i] + n_ZbZab[i] - n_ZabOut[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + n_SYac[i] + n_YaYac[i] + n_YcYac[i] - n_YacOut[i]
  update(Z_ac[]) <- Z_ac[i] + n_SZac[i] + n_ZaZac[i] + n_ZcZac[i] - n_ZacOut[i] 
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + n_SYbc[i] + n_YbYbc[i] + n_YcYbc[i] - n_YbcOut[i]
  update(Z_bc[]) <- Z_bc[i] + n_SZbc[i] + n_ZbZbc[i] + n_ZcZbc[i] - n_ZbcOut[i]

  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + n_SYabc[i] + n_YabYabc[i] + n_YacYabc[i] + n_YbcYabc[i] - n_YabcOut[i]
  update(Z_abc[]) <- Z_abc[i] + n_SZabc[i] + n_ZabZabc[i] + n_ZacZabc[i] + n_ZbcZabc[i] - n_ZabcOut[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  ##STOCHASTIC calculations
  #calculate the probability of leaving S from the rate
  p_SOut[] <- 1 - exp(-sum_all_BI_terms[i])
  
  #calculate the number of people leaving S
  n_SOut[] <- rbinom(S[i], p_SOut[i])
  
  #divide probabilities leaving S 
  p_SYZa[] <- sum_BI_YZa[i]/sum_all_BI_terms[i]
  p_SYZb[] <- sum_BI_YZb[i]/sum_all_BI_terms[i]
  p_SYZc[] <- sum_BI_YZc[i]/sum_all_BI_terms[i]
  p_SYZab[] <- sum_BI_YZab[i]/sum_all_BI_terms[i]
  p_SYZac[] <- sum_BI_YZac[i]/sum_all_BI_terms[i]
  p_SYZbc[] <- sum_BI_YZbc[i]/sum_all_BI_terms[i]
  p_SYZabc[] <- sum_BI_YZabc[i]/sum_all_BI_terms[i] #fixed error here 6/3
  
  #divide the people leaving S into different flows
  n_SYZa[] <- rbinom(n_SOut[i], p_SYZa[i])
  n_SYZb[] <- rbinom(n_SOut[i] - n_SYZa[i], (p_SYZb[i])/(1 - p_SYZa[i] + 1e-8)) #i) subtract those already subtracted from n_SOut ii) renormalize the probability of exit
  n_SYZc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i], (p_SYZc[i])/(1 - p_SYZa[i] - p_SYZb[i] + 1e-8))
  n_SYZab[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i], (p_SYZab[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] + 1e-8))
  n_SYZac[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i], (p_SYZac[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - n_SYZab[i] + 1e-8))
  n_SYZbc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i], (p_SYZbc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - n_SYZab[i] - n_SYZac[i] + 1e-8)) #FIX THIS 6/5
  n_SYZabc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i], (p_SYZabc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] + 1e-8))
  n_SYZs[] <- n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] - n_SYZabc[i]
  
  #split into symptomatic vs asymptomatics
  n_SYs[] <- rbinom(n_SYZs[i], sigma)
  n_SZs[] <- n_SYZs[i] - n_SYs[i] #fixed issue here
  n_SYa[] <- rbinom(n_SYZa[i], sigma)
  n_SZa[] <- n_SYZa[i] - n_SYa[i]
  n_SYb[] <- rbinom(n_SYZb[i], sigma)
  n_SZb[] <- n_SYZb[i] - n_SYb[i]
  n_SYc[] <- rbinom(n_SYZc[i], sigma)
  n_SZc[] <- n_SYZc[i] - n_SYc[i]
  n_SYab[] <- rbinom(n_SYZab[i], sigma)
  n_SZab[] <- n_SYZab[i] - n_SYab[i]
  n_SYac[] <- rbinom(n_SYZac[i], sigma)
  n_SZac[] <- n_SYZac[i] - n_SYac[i]
  n_SYbc[] <- rbinom(n_SYZbc[i], sigma)
  n_SZbc[] <- n_SYZbc[i] - n_SYbc[i]
  n_SYabc[] <- rbinom(n_SYZabc[i], sigma)
  n_SZabc[] <- n_SYZabc[i] - n_SYabc[i]
  
  #calculate the probability of leaving Ys and Zs by converting exit rates
  p_YsOut <- 1 - exp(-(T_s + d))
  p_ZsOut <- 1 - exp(-(T_m + d))
  
  #calculate the number of people leaving Ys and Zs based on the probabilities
  n_YsOut[] <- rbinom(Y_s[i], p_YsOut)
  n_ZsOut[] <- rbinom(Z_s[i], p_ZsOut)
  
  #split the people recovering naturally (nr) vs the people leaving due to treatment (tr)
  p_Ys_tr <- T_s / (T_s + d)
  n_Ys_tr[] <- rbinom(n_YsOut[i], p_Ys_tr)
  n_Ys_nr[] <- n_YsOut[i] - n_Ys_tr[i]
  p_Zs_tr <- T_m / (T_m + d)
  n_Zs_tr[] <- rbinom(n_ZsOut[i], p_Zs_tr)
  n_Zs_nr[] <- n_ZsOut[i] - n_Zs_tr[i]
  
  #of the people being treated, split the flows into recovery vs. resistance
  p_YZs_tr_fail <- E_a*w_a + E_b*w_b + E_c*w_c
  n_Ys_tr_fail[] <- rbinom(n_Ys_tr[i], p_YZs_tr_fail)
  n_Ys_tr_success[] <- n_Ys_tr[i] - n_Ys_tr_fail[i]
  n_Zs_tr_fail[] <- rbinom(n_Zs_tr[i], p_YZs_tr_fail)
  n_Zs_tr_success[] <- n_Zs_tr[i] - n_Zs_tr_fail[i]
  
  #combine those recovering from successful treatment and natural recovery as leaving Ys and Zs and going back to S
  n_YsS[] <- n_Ys_tr_success[i] + n_Ys_nr[i]
  n_ZsS[] <- n_Zs_tr_success[i] + n_Zs_nr[i]
  
  #split those who are treated and develop resistance into the three possible flows
  #first take the individually normalized probabilities of each flow
  p_YZs_YZa <- E_a*w_a/p_YZs_tr_fail
  p_YZs_YZb <- E_b*w_b/p_YZs_tr_fail
  p_YZs_YZc <- E_c*w_c/p_YZs_tr_fail
  
  #do stick-breaking with Y
  n_YsYa[] <- rbinom(n_Ys_tr_fail[i], p_YZs_YZa)
  n_YsYb[] <- rbinom(n_Ys_tr_fail[i] - n_YsYa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_YsYc[] <- n_Ys_tr_fail[i] - n_YsYa[i] - n_YsYb[i]
  
  #do stick-breaking with Z
  n_ZsZa[] <- rbinom(n_Zs_tr_fail[i], p_YZs_YZa)
  n_ZsZb[] <- rbinom(n_Zs_tr_fail[i] - n_ZsZa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_ZsZc[] <- n_Zs_tr_fail[i] - n_ZsZa[i] - n_ZsZb[i]
  
  ##now onto single resistance compartments. Probabilities of exit will be different for Y and Z compartments due to retreatment
  #probabilities for exit flows from Ya, Yb, Yc, Za, Zb, and Zc
  p_YaOut <- 1 - exp(-((E_b + E_c)*T_s + E_a*kappa*T_sr + d))
  p_YbOut <- 1 - exp(-((E_a + E_c)*T_s + E_b*kappa*T_sr + d))
  p_YcOut <- 1 - exp(-((E_a + E_b)*T_s + E_c*kappa*T_sr + d))
  p_ZaOut <- 1 - exp(-((E_b + E_c)*T_m + d))
  p_ZbOut <- 1 - exp(-((E_a + E_c)*T_m + d))
  p_ZcOut <- 1 - exp(-((E_a + E_b)*T_m + d))
  
  #calculate those leaving the single resistance compartments
  n_YaOut[] <- rbinom(Y_a[i], p_YaOut)
  n_YbOut[] <- rbinom(Y_b[i], p_YbOut)
  n_YcOut[] <- rbinom(Y_c[i], p_YcOut)
  n_ZaOut[] <- rbinom(Z_a[i], p_ZaOut)
  n_ZbOut[] <- rbinom(Z_b[i], p_ZbOut)
  n_ZcOut[] <- rbinom(Z_c[i], p_ZcOut)
  
  #calculate the probabilities of leaving single resistance compartments due to natural recovery
  p_Ya_nr <- d/((E_b + E_c)*T_s + E_a*kappa*T_sr + d)
  p_Yb_nr <- d/((E_a + E_c)*T_s + E_b*kappa*T_sr + d)
  p_Yc_nr <- d/((E_a + E_b)*T_s + E_c*kappa*T_sr + d)
  p_Za_nr <- d/((E_b + E_c)*T_m + d)
  p_Zb_nr <- d/((E_a + E_c)*T_m + d)
  p_Zc_nr <- d/((E_a + E_b)*T_m + d)
  
  #split the people leaving single resistance into natural recovery vs treatment
  n_Ya_nr[] <- rbinom(n_YaOut[i], p_Ya_nr)
  n_Ya_tr[] <- n_YaOut[i] - n_Ya_nr[i]
  n_Yb_nr[] <- rbinom(n_YbOut[i], p_Yb_nr)
  n_Yb_tr[] <- n_YbOut[i] - n_Yb_nr[i]
  n_Yc_nr[] <- rbinom(n_YcOut[i], p_Yc_nr)
  n_Yc_tr[] <- n_YcOut[i] - n_Yc_nr[i]
  n_Za_nr[] <- rbinom(n_ZaOut[i], p_Za_nr)
  n_Za_tr[] <- n_ZaOut[i] - n_Za_nr[i]
  n_Zb_nr[] <- rbinom(n_ZbOut[i], p_Zb_nr)
  n_Zb_tr[] <- n_ZbOut[i] - n_Zb_nr[i]
  n_Zc_nr[] <- rbinom(n_ZcOut[i], p_Zc_nr)
  n_Zc_tr[] <- n_ZcOut[i] - n_Zc_nr[i]
  
  #split those leaving single resistance after treatment into those who fail treatment vs. those who recover
  #note that here, treatment failure means people who develop resistance on treatment, not those treated with the drug they are resistance bc they stay in the resistant compartment
  p_Ya_tr_fail <- (E_b*w_b*T_s + E_c*w_c*T_s)/((E_b + E_c)*T_s + E_a*kappa*T_sr)
  p_Yb_tr_fail <- (E_a*w_a*T_s + E_c*w_c*T_s)/((E_a + E_c)*T_s + E_b*kappa*T_sr)
  p_Yc_tr_fail <- (E_a*w_a*T_s + E_b*w_b*T_s)/((E_a + E_b)*T_s + E_c*kappa*T_sr)
  p_Za_tr_fail <- (E_b*w_b + E_c*w_c)/ ((E_b + E_c) + 1e-8) #add very small term to denominator so this probability is always defined
  p_Zb_tr_fail <- (E_a*w_a + E_c*w_c)/ ((E_a + E_c) + 1e-8)
  p_Zc_tr_fail <- (E_a*w_a + E_b*w_b)/ ((E_a + E_b) + 1e-8)
  
  n_Ya_tr_fail[] <- rbinom(n_Ya_tr[i], p_Ya_tr_fail)
  n_Yb_tr_fail[] <- rbinom(n_Yb_tr[i], p_Yb_tr_fail)
  n_Yc_tr_fail[] <- rbinom(n_Yc_tr[i], p_Yc_tr_fail)
  n_Za_tr_fail[] <- rbinom(n_Za_tr[i], p_Za_tr_fail)
  n_Zb_tr_fail[] <- rbinom(n_Zb_tr[i], p_Zb_tr_fail)
  n_Zc_tr_fail[] <- rbinom(n_Zc_tr[i], p_Zc_tr_fail)
  
  n_Ya_tr_success[] <- n_Ya_tr[i] - n_Ya_tr_fail[i]
  n_Yb_tr_success[] <- n_Yb_tr[i] - n_Yb_tr_fail[i]
  n_Yc_tr_success[] <- n_Yc_tr[i] - n_Yc_tr_fail[i]
  n_Za_tr_success[] <- n_Za_tr[i] - n_Za_tr_fail[i]
  n_Zb_tr_success[] <- n_Zb_tr[i] - n_Zb_tr_fail[i]
  n_Zc_tr_success[] <- n_Zc_tr[i] - n_Zc_tr_fail[i]
  
  #combine treatment success and natural recover terms to get single resistance recovery flows
  n_YaS[] <- n_Ya_tr_success[i] + n_Ya_nr[i]
  n_YbS[] <- n_Yb_tr_success[i] + n_Yb_nr[i]
  n_YcS[] <- n_Yc_tr_success[i] + n_Yc_nr[i]
  n_ZaS[] <- n_Za_tr_success[i] + n_Za_nr[i]
  n_ZbS[] <- n_Zb_tr_success[i] + n_Zb_nr[i]
  n_ZcS[] <- n_Zc_tr_success[i] + n_Zc_nr[i]
  
  #now split the treatment failures into the two different resistance compartment flows
  p_YZaYZab <- E_b*w_b / (E_b*w_b + E_c*w_c + 1e-8) #add very small term to denominator so this probability is always defined
  p_YZbYZab <- E_a*w_a / (E_a*w_a + E_c*w_c + 1e-8)
  p_YZcYZac <- E_a*w_a / (E_a*w_a + E_b*w_b + 1e-8)
  
  n_YaYab[] <- rbinom(n_Ya_tr_fail[i], p_YZaYZab)
  n_YaYac[] <- n_Ya_tr_fail[i] - n_YaYab[i]
  n_YbYab[] <- rbinom(n_Yb_tr_fail[i], p_YZbYZab)
  n_YbYbc[] <- n_Yb_tr_fail[i] - n_YbYab[i]
  n_YcYac[] <- rbinom(n_Yc_tr_fail[i], p_YZcYZac)
  n_YcYbc[] <- n_Yc_tr_fail[i] - n_YcYac[i]
  
  n_ZaZab[] <- rbinom(n_Za_tr_fail[i], p_YZaYZab)
  n_ZaZac[] <- n_Za_tr_fail[i] - n_ZaZab[i]
  n_ZbZab[] <- rbinom(n_Zb_tr_fail[i], p_YZbYZab)
  n_ZbZbc[] <- n_Zb_tr_fail[i] - n_ZbZab[i]
  n_ZcZac[] <- rbinom(n_Zc_tr_fail[i], p_YZcYZac)
  n_ZcZbc[] <- n_Zc_tr_fail[i] - n_ZcZac[i]
  
  ##now move to the double resistance compartments
  #first calculate total outflow -- calculate exit probabilities and then draw
  p_YabOut <- 1 - exp(-(E_c*T_s + (1 - E_c)*kappa*T_sr + d))
  p_YacOut <- 1 - exp(-(E_b*T_s + (1 - E_b)*kappa*T_sr + d))
  p_YbcOut <- 1 - exp(-(E_a*T_s + (1 - E_a)*kappa*T_sr + d))
  p_ZabOut <- 1 - exp(-(E_c*T_m + d))
  p_ZacOut <- 1 - exp(-(E_b*T_m + d))
  p_ZbcOut <- 1 - exp(-(E_a*T_m + d))
  
  n_YabOut[] <- rbinom(Y_ab[i], p_YabOut)
  n_YacOut[] <- rbinom(Y_ac[i], p_YacOut)
  n_YbcOut[] <- rbinom(Y_bc[i], p_YbcOut)
  n_ZabOut[] <- rbinom(Z_ab[i], p_ZabOut)
  n_ZacOut[] <- rbinom(Z_ac[i], p_ZacOut)
  n_ZbcOut[] <- rbinom(Z_bc[i], p_ZbcOut)

  #because these are already dual resistance, the only ones leaving and not recovering are developing resistance to the same drug, so separate those and the rest recover
  p_Yab_tr_fail <- (E_c*w_c*T_s)/(E_c*T_s + (1 - E_c)*kappa*T_sr + d)
  p_Yac_tr_fail <- (E_b*w_b*T_s)/(E_b*T_s + (1 - E_b)*kappa*T_sr + d)
  p_Ybc_tr_fail <- (E_a*w_a*T_s)/(E_a*T_s + (1 - E_a)*kappa*T_sr + d)
  p_Zab_tr_fail <- (E_c*w_c*T_m)/(E_c*T_m + d)
  p_Zac_tr_fail <- (E_b*w_b*T_m)/(E_b*T_m + d)
  p_Zbc_tr_fail <- (E_a*w_a*T_m)/(E_a*T_m + d)
  
  #split out resistant flows to Yabc vs recoveries to S
  n_YabYabc[] <- rbinom(n_YabOut[i], p_Yab_tr_fail)
  n_YabS[] <- n_YabOut[i] - n_YabYabc[i]
  n_YacYabc[] <- rbinom(n_YacOut[i], p_Yac_tr_fail)
  n_YacS[] <- n_YacOut[i] - n_YacYabc[i]
  n_YbcYabc[] <- rbinom(n_YbcOut[i], p_Ybc_tr_fail)
  n_YbcS[] <- n_YbcOut[i] - n_YbcYabc[i]
  
  n_ZabZabc[] <- rbinom(n_ZabOut[i], p_Zab_tr_fail)
  n_ZabS[] <- n_ZabOut[i] - n_ZabZabc[i]
  n_ZacZabc[] <- rbinom(n_ZacOut[i], p_Zac_tr_fail)
  n_ZacS[] <- n_ZacOut[i] - n_ZacZabc[i]
  n_ZbcZabc[] <- rbinom(n_ZbcOut[i], p_Zbc_tr_fail)
  n_ZbcS[] <- n_ZbcOut[i] - n_ZbcZabc[i]
  
  #recovery from triple resistance states
  p_YabcOut <- 1 - exp(-(kappa*T_sr + d))
  p_ZabcOut <- 1 - exp(-d)
  n_YabcOut[] <- rbinom(Y_abc[i], p_YabcOut)
  n_ZabcOut[] <- rbinom(Z_abc[i], p_ZabcOut)
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #Binary switches between treatment
  #THIS PART IS UNIQUE TO THE MODEL STRATEGY
  E_a <- (lock_threshold_A < 1)
  E_c <- 1 - (lock_threshold_B < 1)
  E_b <- 1  - E_a - E_c
  
  sum_tx_probs = E_a + E_c + E_b 

  #sum array variables for book-keeping / accounting outputs
  sum_n_SOut <- sum(n_SOut)
  sum_n_SYZa <- sum(n_SYZa)
  sum_n_SYZb <- sum(n_SYZb)
  sum_n_SYZc <- sum(n_SYZc)
  sum_n_SYZab <- sum(n_SYZab)
  sum_n_SYZac <- sum(n_SYZac)
  sum_n_SYZbc <- sum(n_SYZbc)
  sum_n_SYZabc <- sum(n_SYZabc)
  sum_n_SYZs <- sum(n_SYZs)
  sum_n_SYs <- sum(n_SYs)
  sum_n_SZs <- sum(n_SZs)
  sum_n_SYa <- sum(n_SYa)
  sum_n_SZa <- sum(n_SZa)
  sum_n_SYb <- sum(n_SYb)
  sum_n_SZb <- sum(n_SZb)
  sum_n_SYc <- sum(n_SYc)
  sum_n_SZc <- sum(n_SZc)
  sum_n_SYab <- sum(n_SYab)
  sum_n_SZab <- sum(n_SZab)
  sum_n_SYac <- sum(n_SYac)
  sum_n_SZac <- sum(n_SZac)
  sum_n_SYbc <- sum(n_SYbc)
  sum_n_SZbc <- sum(n_SZbc)
  sum_n_SYabc <- sum(n_SYabc)
  sum_n_SZabc <- sum(n_SZabc)
  sum_n_YsOut <- sum(n_YsOut)
  sum_n_ZsOut <- sum(n_ZsOut)
  sum_n_Ys_tr <- sum(n_Ys_tr)
  sum_n_Ys_nr <- sum(n_Ys_nr)
  sum_n_Zs_tr <- sum(n_Zs_tr)
  sum_n_Zs_nr <- sum(n_Zs_nr)
  sum_n_Ys_tr_fail <- sum(n_Ys_tr_fail)
  sum_n_Ys_tr_success <- sum(n_Ys_tr_success)
  sum_n_Zs_tr_fail <- sum(n_Zs_tr_fail)
  sum_n_Zs_tr_success <- sum(n_Zs_tr_success)
  sum_n_YsS <- sum(n_YsS)
  sum_n_ZsS <- sum(n_ZsS)
  sum_n_YsYa <- sum(n_YsYa)
  sum_n_YsYb <- sum(n_YsYb)
  sum_n_YsYc <- sum(n_YsYc)
  sum_n_ZsZa <- sum(n_ZsZa)
  sum_n_ZsZb <- sum(n_ZsZb)
  sum_n_ZsZc <- sum(n_ZsZc)
  sum_n_YaOut <- sum(n_YaOut)
  sum_n_YbOut <- sum(n_YbOut)
  sum_n_YcOut <- sum(n_YcOut)
  sum_n_ZaOut <- sum(n_ZaOut)
  sum_n_ZbOut <- sum(n_ZbOut)
  sum_n_ZcOut <- sum(n_ZcOut)
  sum_n_Ya_nr <- sum(n_Ya_nr)
  sum_n_Ya_tr <- sum(n_Ya_tr)
  sum_n_Yb_nr <- sum(n_Yb_nr)
  sum_n_Yb_tr <- sum(n_Yb_tr)
  sum_n_Yc_nr <- sum(n_Yc_nr)
  sum_n_Yc_tr <- sum(n_Yc_tr)
  sum_n_Za_nr <- sum(n_Za_nr)
  sum_n_Za_tr <- sum(n_Za_tr)
  sum_n_Zb_nr <- sum(n_Zb_nr)
  sum_n_Zb_tr <- sum(n_Zb_tr)
  sum_n_Zc_nr <- sum(n_Zc_nr)
  sum_n_Zc_tr <- sum(n_Zc_tr)
  sum_n_Ya_tr_fail <- sum(n_Ya_tr_fail)
  sum_n_Yb_tr_fail <- sum(n_Yb_tr_fail)
  sum_n_Yc_tr_fail <- sum(n_Yc_tr_fail)
  sum_n_Za_tr_fail <- sum(n_Za_tr_fail)
  sum_n_Zb_tr_fail <- sum(n_Zb_tr_fail)
  sum_n_Zc_tr_fail <- sum(n_Zc_tr_fail)
  sum_n_Ya_tr_success <- sum(n_Ya_tr_success)
  sum_n_Yb_tr_success <- sum(n_Yb_tr_success)
  sum_n_Yc_tr_success <- sum(n_Yc_tr_success)
  sum_n_Za_tr_success <- sum(n_Za_tr_success)
  sum_n_Zb_tr_success <- sum(n_Zb_tr_success)
  sum_n_Zc_tr_success <- sum(n_Zc_tr_success)
  sum_n_YaS <- sum(n_YaS)
  sum_n_YbS <- sum(n_YbS)
  sum_n_YcS <- sum(n_YcS)
  sum_n_ZaS <- sum(n_ZaS)
  sum_n_ZbS <- sum(n_ZbS)
  sum_n_ZcS <- sum(n_ZcS)
  sum_n_YaYab <- sum(n_YaYab)
  sum_n_YaYac <- sum(n_YaYac)
  sum_n_YbYab <- sum(n_YbYab)
  sum_n_YbYbc <- sum(n_YbYbc)
  sum_n_YcYac <- sum(n_YcYac)
  sum_n_YcYbc <- sum(n_YcYbc)
  sum_n_ZaZab <- sum(n_ZaZab)
  sum_n_ZaZac <- sum(n_ZaZac)
  sum_n_ZbZab <- sum(n_ZbZab)
  sum_n_ZbZbc <- sum(n_ZbZbc)
  sum_n_ZcZac <- sum(n_ZcZac)
  sum_n_ZcZbc <- sum(n_ZcZbc)
  sum_n_YabOut <- sum(n_YabOut)
  sum_n_YacOut <- sum(n_YacOut)
  sum_n_YbcOut <- sum(n_YbcOut)
  sum_n_ZabOut <- sum(n_ZabOut)
  sum_n_ZacOut <- sum(n_ZacOut)
  sum_n_ZbcOut <- sum(n_ZbcOut)
  sum_n_YabYabc <- sum(n_YabYabc)
  sum_n_YabS <- sum(n_YabS)
  sum_n_YacYabc <- sum(n_YacYabc)
  sum_n_YacS <- sum(n_YacS)
  sum_n_YbcYabc <- sum(n_YbcYabc)
  sum_n_YbcS<- sum(n_YbcS)
  sum_n_ZabZabc <- sum(n_ZabZabc)
  sum_n_ZabS <- sum(n_ZabS)
  sum_n_ZacZabc <- sum(n_ZacZabc)
  sum_n_ZacS <- sum(n_ZacS)
  sum_n_ZbcZabc <- sum(n_ZbcZabc)
  sum_n_ZbcS<- sum(n_ZbcS)
  sum_n_YabcOut <- sum(n_YabcOut)
  sum_n_ZabcOut <- sum(n_ZabcOut)
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  dim(n_SOut) <- N_risk
  dim(p_SOut) <- N_risk
  dim(p_SYZa) <- N_risk
  dim(p_SYZb) <- N_risk
  dim(p_SYZc) <- N_risk
  dim(p_SYZab) <- N_risk
  dim(p_SYZac) <- N_risk  
  dim(p_SYZbc) <- N_risk
  dim(p_SYZabc) <- N_risk
  dim(n_SYZa) <- N_risk
  dim(n_SYZb) <- N_risk
  dim(n_SYZc) <- N_risk
  dim(n_SYZab) <- N_risk
  dim(n_SYZac) <- N_risk
  dim(n_SYZbc) <- N_risk
  dim(n_SYZabc) <- N_risk
  dim(n_SYZs) <- N_risk
  dim(n_SYs) <- N_risk
  dim(n_SZs) <- N_risk
  dim(n_SYa) <- N_risk
  dim(n_SZa) <- N_risk
  dim(n_SYb) <- N_risk
  dim(n_SZb) <- N_risk
  dim(n_SYc) <- N_risk
  dim(n_SZc) <- N_risk
  dim(n_SYab) <- N_risk
  dim(n_SZab) <- N_risk
  dim(n_SYac) <- N_risk
  dim(n_SZac) <- N_risk
  dim(n_SYbc) <- N_risk
  dim(n_SZbc) <- N_risk
  dim(n_SYabc) <- N_risk
  dim(n_SZabc) <- N_risk
  dim(n_YsOut) <- N_risk
  dim(n_ZsOut) <- N_risk
  dim(n_Ys_tr) <- N_risk
  dim(n_Ys_nr) <- N_risk
  dim(n_Zs_tr) <- N_risk
  dim(n_Zs_nr) <- N_risk
  dim(n_Ys_tr_fail) <- N_risk
  dim(n_Ys_tr_success) <- N_risk
  dim(n_Zs_tr_fail) <- N_risk
  dim(n_Zs_tr_success) <- N_risk
  dim(n_YsS) <- N_risk
  dim(n_ZsS) <- N_risk
  dim(n_YsYa) <- N_risk
  dim(n_YsYb) <- N_risk
  dim(n_YsYc) <- N_risk
  dim(n_ZsZa) <- N_risk
  dim(n_ZsZb) <- N_risk
  dim(n_ZsZc) <- N_risk
  dim(n_YaOut) <- N_risk
  dim(n_YbOut) <- N_risk
  dim(n_YcOut) <- N_risk
  dim(n_ZaOut) <- N_risk
  dim(n_ZbOut) <- N_risk
  dim(n_ZcOut) <- N_risk
  dim(n_Ya_nr) <- N_risk
  dim(n_Ya_tr) <- N_risk
  dim(n_Yb_nr) <- N_risk
  dim(n_Yb_tr) <- N_risk
  dim(n_Yc_nr) <- N_risk
  dim(n_Yc_tr) <- N_risk
  dim(n_Za_nr) <- N_risk
  dim(n_Za_tr) <- N_risk
  dim(n_Zb_nr) <- N_risk
  dim(n_Zb_tr) <- N_risk
  dim(n_Zc_nr) <- N_risk
  dim(n_Zc_tr) <- N_risk
  dim(n_Ya_tr_fail) <- N_risk
  dim(n_Yb_tr_fail) <- N_risk
  dim(n_Yc_tr_fail) <- N_risk
  dim(n_Za_tr_fail) <- N_risk
  dim(n_Zb_tr_fail) <- N_risk
  dim(n_Zc_tr_fail) <- N_risk
  dim(n_Ya_tr_success) <- N_risk
  dim(n_Yb_tr_success) <- N_risk
  dim(n_Yc_tr_success) <- N_risk
  dim(n_Za_tr_success) <- N_risk
  dim(n_Zb_tr_success) <- N_risk
  dim(n_Zc_tr_success) <- N_risk
  dim(n_YaS) <- N_risk
  dim(n_YbS) <- N_risk
  dim(n_YcS) <- N_risk
  dim(n_ZaS) <- N_risk
  dim(n_ZbS) <- N_risk
  dim(n_ZcS) <- N_risk
  dim(n_YaYab) <- N_risk
  dim(n_YaYac) <- N_risk
  dim(n_YbYab) <- N_risk
  dim(n_YbYbc) <- N_risk
  dim(n_YcYac) <- N_risk
  dim(n_YcYbc) <- N_risk
  dim(n_ZaZab) <- N_risk
  dim(n_ZaZac) <- N_risk
  dim(n_ZbZab) <- N_risk
  dim(n_ZbZbc) <- N_risk
  dim(n_ZcZac) <- N_risk
  dim(n_ZcZbc) <- N_risk
  dim(n_YabOut) <- N_risk
  dim(n_YacOut) <- N_risk
  dim(n_YbcOut) <- N_risk
  dim(n_ZabOut) <- N_risk
  dim(n_ZacOut) <- N_risk
  dim(n_ZbcOut) <- N_risk
  dim(n_YabYabc) <- N_risk
  dim(n_YabS) <- N_risk
  dim(n_YacYabc) <- N_risk
  dim(n_YacS) <- N_risk
  dim(n_YbcYabc) <- N_risk
  dim(n_YbcS) <- N_risk
  dim(n_ZabZabc) <- N_risk
  dim(n_ZabS) <- N_risk
  dim(n_ZacZabc) <- N_risk
  dim(n_ZacS) <- N_risk
  dim(n_ZbcZabc) <- N_risk
  dim(n_ZbcS) <- N_risk
  dim(n_YabcOut) <- N_risk
  dim(n_ZabcOut) <- N_risk

  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev
  output(sum_tx_probs) <- sum_tx_probs
  output(p_SOut) <- p_SOut
  output(sum_n_SOut) <- sum_n_SOut
  output(p_SYZa) <- p_SYZa
  output(p_SYZb) <- p_SYZb
  output(p_SYZc) <- p_SYZc
  output(p_SYZab) <- p_SYZab
  output(p_SYZac) <- p_SYZac
  output(p_SYZbc) <- p_SYZbc
  output(p_SYZabc) <- p_SYZabc
  output(sum_n_SYZa) <- sum_n_SYZa
  output(sum_n_SYZb) <- sum_n_SYZb
  output(sum_n_SYZc) <- sum_n_SYZc
  output(sum_n_SYZab) <- sum_n_SYZab
  output(sum_n_SYZac) <- sum_n_SYZac
  output(sum_n_SYZbc) <- sum_n_SYZbc
  output(sum_n_SYZabc) <- sum_n_SYZabc
  output(sum_n_SYZs) <- sum_n_SYZs
  output(sum_n_SYs) <- sum_n_SYs
  output(sum_n_SZs) <- sum_n_SZs
  output(sum_n_SYa) <- sum_n_SYa
  output(sum_n_SZa) <- sum_n_SZa
  output(sum_n_SYb) <- sum_n_SYb
  output(sum_n_SZb) <- sum_n_SZb
  output(sum_n_SYc) <- sum_n_SYc
  output(sum_n_SZc) <- sum_n_SZc
  output(sum_n_SYab) <- sum_n_SYab
  output(sum_n_SZab) <- sum_n_SZab
  output(sum_n_SYac) <- sum_n_SYac
  output(sum_n_SZac) <- sum_n_SZac
  output(sum_n_SYbc) <- sum_n_SYbc
  output(sum_n_SZbc) <- sum_n_SZbc
  output(sum_n_SYabc) <- sum_n_SYabc
  output(sum_n_SZabc) <- sum_n_SZabc
  output(p_YsOut) <- p_YsOut
  output(p_ZsOut) <- p_ZsOut
  output(p_Ys_tr) <- p_Ys_tr
  output(p_Zs_tr) <- p_Zs_tr
  output(p_YZs_tr_fail) <- p_YZs_tr_fail
  output(p_YZs_YZa) <- p_YZs_YZa
  output(p_YZs_YZb) <- p_YZs_YZb
  output(p_YZs_YZc) <- p_YZs_YZc
  output(p_YaOut) <- p_YaOut
  output(p_YbOut) <- p_YbOut
  output(p_YcOut) <- p_YcOut
  output(p_ZaOut) <- p_ZaOut
  output(p_ZbOut) <- p_ZbOut
  output(p_ZcOut) <- p_ZcOut
  output(p_Ya_nr) <- p_Ya_nr
  output(p_Yb_nr) <- p_Yb_nr
  output(p_Yc_nr) <- p_Yc_nr
  output(p_Za_nr) <- p_Za_nr
  output(p_Zb_nr) <- p_Zb_nr
  output(p_Zc_nr) <- p_Zc_nr
  output(p_Ya_tr_fail) <- p_Ya_tr_fail
  output(p_Yb_tr_fail) <- p_Yb_tr_fail
  output(p_Yc_tr_fail) <- p_Yc_tr_fail
  output(p_Za_tr_fail) <- p_Za_tr_fail
  output(p_Zb_tr_fail) <- p_Zb_tr_fail
  output(p_Zc_tr_fail) <- p_Zc_tr_fail
  output(p_YZaYZab) <- p_YZaYZab
  output(p_YZbYZab) <- p_YZbYZab
  output(p_YZcYZac) <- p_YZcYZac
  output(p_YabOut) <- p_YabOut 
  output(p_YacOut) <- p_YacOut 
  output(p_YbcOut) <- p_YbcOut 
  output(p_ZabOut) <- p_ZabOut 
  output(p_ZacOut) <- p_ZacOut 
  output(p_ZbcOut) <- p_ZbcOut
  output(p_Yab_tr_fail) <- p_Yab_tr_fail
  output(p_Yac_tr_fail) <- p_Yac_tr_fail
  output(p_Ybc_tr_fail) <- p_Ybc_tr_fail
  output(p_Zab_tr_fail) <- p_Zab_tr_fail
  output(p_Zac_tr_fail) <- p_Zac_tr_fail
  output(p_Zbc_tr_fail) <- p_Zbc_tr_fail
  output(sum_n_YsOut) <- sum_n_YsOut
  output(sum_n_ZsOut) <- sum_n_ZsOut
  output(sum_n_Ys_tr) <- sum_n_Ys_tr
  output(sum_n_Ys_nr) <- sum_n_Ys_nr
  output(sum_n_Zs_tr) <- sum_n_Zs_tr
  output(sum_n_Zs_nr) <- sum_n_Zs_nr
  output(sum_n_Ys_tr_fail) <- sum_n_Ys_tr_fail
  output(sum_n_Ys_tr_success) <- sum_n_Ys_tr_success
  output(sum_n_Zs_tr_fail) <- sum_n_Zs_tr_fail
  output(sum_n_Zs_tr_success) <- sum_n_Zs_tr_success
  output(sum_n_YsS) <- sum_n_YsS
  output(sum_n_ZsS) <- sum_n_ZsS
  output(sum_n_YsYa) <- sum_n_YsYa
  output(sum_n_YsYb) <- sum_n_YsYb
  output(sum_n_YsYc) <- sum_n_YsYc
  output(sum_n_ZsZa) <- sum_n_ZsZa
  output(sum_n_ZsZb) <- sum_n_ZsZb
  output(sum_n_ZsZc) <- sum_n_ZsZc
  output(sum_n_YaOut) <- sum_n_YaOut
  output(sum_n_YbOut) <- sum_n_YbOut
  output(sum_n_YcOut) <- sum_n_YcOut
  output(sum_n_ZaOut) <- sum_n_ZaOut
  output(sum_n_ZbOut) <- sum_n_ZbOut
  output(sum_n_ZcOut) <- sum_n_ZcOut
  output(sum_n_Ya_nr) <- sum_n_Ya_nr
  output(sum_n_Ya_tr) <- sum_n_Ya_tr
  output(sum_n_Yb_nr) <- sum_n_Yb_nr
  output(sum_n_Yb_tr) <- sum_n_Yb_tr
  output(sum_n_Yc_nr) <- sum_n_Yc_nr
  output(sum_n_Yc_tr) <- sum_n_Yc_tr
  output(sum_n_Za_nr) <- sum_n_Za_nr
  output(sum_n_Za_tr) <- sum_n_Za_tr
  output(sum_n_Zb_nr) <- sum_n_Zb_nr
  output(sum_n_Zb_tr) <- sum_n_Zb_tr
  output(sum_n_Zc_nr) <- sum_n_Zc_nr
  output(sum_n_Zc_tr) <- sum_n_Zc_tr
  output(sum_n_Ya_tr_fail) <- sum_n_Ya_tr_fail
  output(sum_n_Yb_tr_fail) <- sum_n_Yb_tr_fail
  output(sum_n_Yc_tr_fail) <- sum_n_Yc_tr_fail
  output(sum_n_Za_tr_fail) <- sum_n_Za_tr_fail
  output(sum_n_Zb_tr_fail) <- sum_n_Zb_tr_fail
  output(sum_n_Zc_tr_fail) <- sum_n_Zc_tr_fail
  output(sum_n_Ya_tr_success) <- sum_n_Ya_tr_success
  output(sum_n_Yb_tr_success) <- sum_n_Yb_tr_success
  output(sum_n_Yc_tr_success) <- sum_n_Yc_tr_success
  output(sum_n_Za_tr_success) <- sum_n_Za_tr_success
  output(sum_n_Zb_tr_success) <- sum_n_Zb_tr_success
  output(sum_n_Zc_tr_success) <- sum_n_Zc_tr_success
  output(sum_n_YaS) <- sum_n_YaS
  output(sum_n_YbS) <- sum_n_YbS
  output(sum_n_YcS) <- sum_n_YcS
  output(sum_n_ZaS) <- sum_n_ZaS
  output(sum_n_ZbS) <- sum_n_ZbS
  output(sum_n_ZcS) <- sum_n_ZcS
  output(sum_n_YaYab) <- sum_n_YaYab
  output(sum_n_YaYac) <- sum_n_YaYac
  output(sum_n_YbYab) <- sum_n_YbYab
  output(sum_n_YbYbc) <- sum_n_YbYbc
  output(sum_n_YcYac) <- sum_n_YcYac 
  output(sum_n_YcYbc) <- sum_n_YcYbc 
  output(sum_n_ZaZab) <- sum_n_ZaZab
  output(sum_n_ZaZac) <- sum_n_ZaZac
  output(sum_n_ZbZab) <- sum_n_ZbZab
  output(sum_n_ZbZbc) <- sum_n_ZbZbc
  output(sum_n_ZcZac) <- sum_n_ZcZac 
  output(sum_n_ZcZbc) <- sum_n_ZcZbc 
  output(sum_n_YabOut) <- sum_n_YabOut 
  output(sum_n_YacOut) <- sum_n_YacOut 
  output(sum_n_YbcOut) <- sum_n_YbcOut
  output(sum_n_ZabOut) <- sum_n_ZabOut 
  output(sum_n_ZacOut) <- sum_n_ZacOut 
  output(sum_n_ZbcOut) <- sum_n_ZbcOut
  output(sum_n_YabYabc) <- sum_n_YabYabc
  output(sum_n_YabS) <- sum_n_YabS
  output(sum_n_YacYabc) <- sum_n_YacYabc
  output(sum_n_YacS) <- sum_n_YacS
  output(sum_n_YbcYabc) <- sum_n_YbcYabc
  output(sum_n_YbcS) <- sum_n_YbcS
  output(sum_n_ZabZabc) <- sum_n_ZabZabc
  output(sum_n_ZabS) <- sum_n_ZabS
  output(sum_n_ZacZabc) <- sum_n_ZacZabc
  output(sum_n_ZacS) <- sum_n_ZacS
  output(sum_n_ZbcZabc) <- sum_n_ZbcZabc
  output(sum_n_ZbcS) <- sum_n_ZbcS
  output(sum_n_ZabcOut) <- sum_n_ZabcOut 
  output(p_YabcOut) <- p_YabcOut
  output(p_ZabcOut) <- p_ZabcOut
  output(sum_n_YabcOut) <- sum_n_YabcOut
})

equal_allocation_stochastic_6_3 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] - n_SOut[i] + n_YsS[i] + n_YaS[i] + n_YbS[i] + n_YcS[i] + n_YabS[i] + n_YacS[i] + n_YbcS[i] +
    n_YabcOut[i] + n_ZsS[i] + n_ZaS[i] + n_ZbS[i] + n_ZcS[i] + n_ZabS[i] + n_ZacS[i] + n_ZbcS[i] + n_ZabcOut[i] 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + n_SYs[i] - n_YsOut[i]
  update(Z_s[]) <- Z_s[i] + n_SZs[i] - n_ZsOut[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + n_SYa[i] + n_YsYa[i] - n_YaOut[i]
  update(Z_a[]) <- Z_a[i] + n_SZa[i] + n_ZsZa[i] - n_ZaOut[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + n_SYb[i] + n_YsYb[i] - n_YbOut[i]
  update(Z_b[]) <- Z_b[i] + n_SZb[i] + n_ZsZb[i] - n_ZbOut[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + n_SYc[i] + n_YsYc[i] - n_YcOut[i] 
  update(Z_c[]) <- Z_c[i] + n_SZc[i] + n_ZsZc[i] - n_ZcOut[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + n_SYab[i] + n_YaYab[i] + n_YbYab[i] - n_YabOut[i]
  update(Z_ab[]) <- Z_ab[i] + n_SZab[i] + n_ZaZab[i] + n_ZbZab[i] - n_ZabOut[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + n_SYac[i] + n_YaYac[i] + n_YcYac[i] - n_YacOut[i]
  update(Z_ac[]) <- Z_ac[i] + n_SZac[i] + n_ZaZac[i] + n_ZcZac[i] - n_ZacOut[i] 
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + n_SYbc[i] + n_YbYbc[i] + n_YcYbc[i] - n_YbcOut[i]
  update(Z_bc[]) <- Z_bc[i] + n_SZbc[i] + n_ZbZbc[i] + n_ZcZbc[i] - n_ZbcOut[i]
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + n_SYabc[i] + n_YabYabc[i] + n_YacYabc[i] + n_YbcYabc[i] - n_YabcOut[i]
  update(Z_abc[]) <- Z_abc[i] + n_SZabc[i] + n_ZabZabc[i] + n_ZacZabc[i] + n_ZbcZabc[i] - n_ZabcOut[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  ##STOCHASTIC calculations
  #calculate the probability of leaving S from the rate
  p_SOut[] <- 1 - exp(-sum_all_BI_terms[i])
  
  #calculate the number of people leaving S
  n_SOut[] <- rbinom(S[i], p_SOut[i])
  
  #divide probabilities leaving S 
  p_SYZa[] <- sum_BI_YZa[i]/sum_all_BI_terms[i]
  p_SYZb[] <- sum_BI_YZb[i]/sum_all_BI_terms[i]
  p_SYZc[] <- sum_BI_YZc[i]/sum_all_BI_terms[i]
  p_SYZab[] <- sum_BI_YZab[i]/sum_all_BI_terms[i]
  p_SYZac[] <- sum_BI_YZac[i]/sum_all_BI_terms[i]
  p_SYZbc[] <- sum_BI_YZbc[i]/sum_all_BI_terms[i]
  p_SYZabc[] <- sum_BI_YZabc[i]/sum_all_BI_terms[i] #fixed error here 6/3
  
  #divide the people leaving S into different flows
  n_SYZa[] <- rbinom(n_SOut[i], p_SYZa[i])
  n_SYZb[] <- rbinom(n_SOut[i] - n_SYZa[i], (p_SYZb[i])/(1 - p_SYZa[i] + 1e-8)) #i) subtract those already subtracted from n_SOut ii) renormalize the probability of exit
  n_SYZc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i], (p_SYZc[i])/(1 - p_SYZa[i] - p_SYZb[i] + 1e-8))
  n_SYZab[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i], (p_SYZab[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] + 1e-8))
  n_SYZac[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i], (p_SYZac[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - n_SYZab[i] + 1e-8))
  n_SYZbc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i], (p_SYZbc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - n_SYZab[i] - n_SYZac[i] + 1e-8))
  n_SYZabc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i], (p_SYZabc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] + 1e-8))
  n_SYZs[] <- n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] - n_SYZabc[i]
  
  #split into symptomatic vs asymptomatics
  n_SYs[] <- rbinom(n_SYZs[i], sigma)
  n_SZs[] <- n_SYZs[i] - n_SYs[i] #fixed issue here
  n_SYa[] <- rbinom(n_SYZa[i], sigma)
  n_SZa[] <- n_SYZa[i] - n_SYa[i]
  n_SYb[] <- rbinom(n_SYZb[i], sigma)
  n_SZb[] <- n_SYZb[i] - n_SYb[i]
  n_SYc[] <- rbinom(n_SYZc[i], sigma)
  n_SZc[] <- n_SYZc[i] - n_SYc[i]
  n_SYab[] <- rbinom(n_SYZab[i], sigma)
  n_SZab[] <- n_SYZab[i] - n_SYab[i]
  n_SYac[] <- rbinom(n_SYZac[i], sigma)
  n_SZac[] <- n_SYZac[i] - n_SYac[i]
  n_SYbc[] <- rbinom(n_SYZbc[i], sigma)
  n_SZbc[] <- n_SYZbc[i] - n_SYbc[i]
  n_SYabc[] <- rbinom(n_SYZabc[i], sigma)
  n_SZabc[] <- n_SYZabc[i] - n_SYabc[i]
  
  #calculate the probability of leaving Ys and Zs by converting exit rates
  p_YsOut <- 1 - exp(-(T_s + d))
  p_ZsOut <- 1 - exp(-(T_m + d))
  
  #calculate the number of people leaving Ys and Zs based on the probabilities
  n_YsOut[] <- rbinom(Y_s[i], p_YsOut)
  n_ZsOut[] <- rbinom(Z_s[i], p_ZsOut)
  
  #split the people recovering naturally (nr) vs the people leaving due to treatment (tr)
  p_Ys_tr <- T_s / (T_s + d)
  n_Ys_tr[] <- rbinom(n_YsOut[i], p_Ys_tr)
  n_Ys_nr[] <- n_YsOut[i] - n_Ys_tr[i]
  p_Zs_tr <- T_m / (T_m + d)
  n_Zs_tr[] <- rbinom(n_ZsOut[i], p_Zs_tr)
  n_Zs_nr[] <- n_ZsOut[i] - n_Zs_tr[i]
  
  #of the people being treated, split the flows into recovery vs. resistance
  p_YZs_tr_fail <- E_a*w_a + E_b*w_b + E_c*w_c
  n_Ys_tr_fail[] <- rbinom(n_Ys_tr[i], p_YZs_tr_fail)
  n_Ys_tr_success[] <- n_Ys_tr[i] - n_Ys_tr_fail[i]
  n_Zs_tr_fail[] <- rbinom(n_Zs_tr[i], p_YZs_tr_fail)
  n_Zs_tr_success[] <- n_Zs_tr[i] - n_Zs_tr_fail[i]
  
  #combine those recovering from successful treatment and natural recovery as leaving Ys and Zs and going back to S
  n_YsS[] <- n_Ys_tr_success[i] + n_Ys_nr[i]
  n_ZsS[] <- n_Zs_tr_success[i] + n_Zs_nr[i]
  
  #split those who are treated and develop resistance into the three possible flows
  #first take the individually normalized probabilities of each flow
  p_YZs_YZa <- E_a*w_a/p_YZs_tr_fail
  p_YZs_YZb <- E_b*w_b/p_YZs_tr_fail
  p_YZs_YZc <- E_c*w_c/p_YZs_tr_fail
  
  #do stick-breaking with Y
  n_YsYa[] <- rbinom(n_Ys_tr_fail[i], p_YZs_YZa)
  n_YsYb[] <- rbinom(n_Ys_tr_fail[i] - n_YsYa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_YsYc[] <- n_Ys_tr_fail[i] - n_YsYa[i] - n_YsYb[i]
  
  #do stick-breaking with Z
  n_ZsZa[] <- rbinom(n_Zs_tr_fail[i], p_YZs_YZa)
  n_ZsZb[] <- rbinom(n_Zs_tr_fail[i] - n_ZsZa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_ZsZc[] <- n_Zs_tr_fail[i] - n_ZsZa[i] - n_ZsZb[i]
  
  ##now onto single resistance compartments. Probabilities of exit will be different for Y and Z compartments due to retreatment
  #probabilities for exit flows from Ya, Yb, Yc, Za, Zb, and Zc
  p_YaOut <- 1 - exp(-((E_b + E_c)*T_s + E_a*kappa*T_sr + d))
  p_YbOut <- 1 - exp(-((E_a + E_c)*T_s + E_b*kappa*T_sr + d))
  p_YcOut <- 1 - exp(-((E_a + E_b)*T_s + E_c*kappa*T_sr + d))
  p_ZaOut <- 1 - exp(-((E_b + E_c)*T_m + d))
  p_ZbOut <- 1 - exp(-((E_a + E_c)*T_m + d))
  p_ZcOut <- 1 - exp(-((E_a + E_b)*T_m + d))
  
  #calculate those leaving the single resistance compartments
  n_YaOut[] <- rbinom(Y_a[i], p_YaOut)
  n_YbOut[] <- rbinom(Y_b[i], p_YbOut)
  n_YcOut[] <- rbinom(Y_c[i], p_YcOut)
  n_ZaOut[] <- rbinom(Z_a[i], p_ZaOut)
  n_ZbOut[] <- rbinom(Z_b[i], p_ZbOut)
  n_ZcOut[] <- rbinom(Z_c[i], p_ZcOut)
  
  #calculate the probabilities of leaving single resistance compartments due to natural recovery
  p_Ya_nr <- d/((E_b + E_c)*T_s + E_a*kappa*T_sr + d)
  p_Yb_nr <- d/((E_a + E_c)*T_s + E_b*kappa*T_sr + d)
  p_Yc_nr <- d/((E_a + E_b)*T_s + E_c*kappa*T_sr + d)
  p_Za_nr <- d/((E_b + E_c)*T_m + d)
  p_Zb_nr <- d/((E_a + E_c)*T_m + d)
  p_Zc_nr <- d/((E_a + E_b)*T_m + d)
  
  #split the people leaving single resistance into natural recovery vs treatment
  n_Ya_nr[] <- rbinom(n_YaOut[i], p_Ya_nr)
  n_Ya_tr[] <- n_YaOut[i] - n_Ya_nr[i]
  n_Yb_nr[] <- rbinom(n_YbOut[i], p_Yb_nr)
  n_Yb_tr[] <- n_YbOut[i] - n_Yb_nr[i]
  n_Yc_nr[] <- rbinom(n_YcOut[i], p_Yc_nr)
  n_Yc_tr[] <- n_YcOut[i] - n_Yc_nr[i]
  n_Za_nr[] <- rbinom(n_ZaOut[i], p_Za_nr)
  n_Za_tr[] <- n_ZaOut[i] - n_Za_nr[i]
  n_Zb_nr[] <- rbinom(n_ZbOut[i], p_Zb_nr)
  n_Zb_tr[] <- n_ZbOut[i] - n_Zb_nr[i]
  n_Zc_nr[] <- rbinom(n_ZcOut[i], p_Zc_nr)
  n_Zc_tr[] <- n_ZcOut[i] - n_Zc_nr[i]
  
  #split those leaving single resistance after treatment into those who fail treatment vs. those who recover
  #note that here, treatment failure means people who develop resistance on treatment, not those treated with the drug they are resistance bc they stay in the resistant compartment
  p_Ya_tr_fail <- (E_b*w_b*T_s + E_c*w_c*T_s)/((E_b + E_c)*T_s + E_a*kappa*T_sr)
  p_Yb_tr_fail <- (E_a*w_a*T_s + E_c*w_c*T_s)/((E_a + E_c)*T_s + E_b*kappa*T_sr)
  p_Yc_tr_fail <- (E_a*w_a*T_s + E_b*w_b*T_s)/((E_a + E_b)*T_s + E_c*kappa*T_sr)
  p_Za_tr_fail <- (E_b*w_b + E_c*w_c)/ ((E_b + E_c) + 1e-8) #add very small term to denominator so this probability is always defined
  p_Zb_tr_fail <- (E_a*w_a + E_c*w_c)/ ((E_a + E_c) + 1e-8)
  p_Zc_tr_fail <- (E_a*w_a + E_b*w_b)/ ((E_a + E_b) + 1e-8)
  
  n_Ya_tr_fail[] <- rbinom(n_Ya_tr[i], p_Ya_tr_fail)
  n_Yb_tr_fail[] <- rbinom(n_Yb_tr[i], p_Yb_tr_fail)
  n_Yc_tr_fail[] <- rbinom(n_Yc_tr[i], p_Yc_tr_fail)
  n_Za_tr_fail[] <- rbinom(n_Za_tr[i], p_Za_tr_fail)
  n_Zb_tr_fail[] <- rbinom(n_Zb_tr[i], p_Zb_tr_fail)
  n_Zc_tr_fail[] <- rbinom(n_Zc_tr[i], p_Zc_tr_fail)
  
  n_Ya_tr_success[] <- n_Ya_tr[i] - n_Ya_tr_fail[i]
  n_Yb_tr_success[] <- n_Yb_tr[i] - n_Yb_tr_fail[i]
  n_Yc_tr_success[] <- n_Yc_tr[i] - n_Yc_tr_fail[i]
  n_Za_tr_success[] <- n_Za_tr[i] - n_Za_tr_fail[i]
  n_Zb_tr_success[] <- n_Zb_tr[i] - n_Zb_tr_fail[i]
  n_Zc_tr_success[] <- n_Zc_tr[i] - n_Zc_tr_fail[i]
  
  #combine treatment success and natural recover terms to get single resistance recovery flows
  n_YaS[] <- n_Ya_tr_success[i] + n_Ya_nr[i]
  n_YbS[] <- n_Yb_tr_success[i] + n_Yb_nr[i]
  n_YcS[] <- n_Yc_tr_success[i] + n_Yc_nr[i]
  n_ZaS[] <- n_Za_tr_success[i] + n_Za_nr[i]
  n_ZbS[] <- n_Zb_tr_success[i] + n_Zb_nr[i]
  n_ZcS[] <- n_Zc_tr_success[i] + n_Zc_nr[i]
  
  #now split the treatment failures into the two different resistance compartment flows
  p_YZaYZab <- E_b*w_b / (E_b*w_b + E_c*w_c + 1e-8) #add very small term to denominator so this probability is always defined
  p_YZbYZab <- E_a*w_a / (E_a*w_a + E_c*w_c + 1e-8)
  p_YZcYZac <- E_a*w_a / (E_a*w_a + E_b*w_b + 1e-8)
  
  n_YaYab[] <- rbinom(n_Ya_tr_fail[i], p_YZaYZab)
  n_YaYac[] <- n_Ya_tr_fail[i] - n_YaYab[i]
  n_YbYab[] <- rbinom(n_Yb_tr_fail[i], p_YZbYZab)
  n_YbYbc[] <- n_Yb_tr_fail[i] - n_YbYab[i]
  n_YcYac[] <- rbinom(n_Yc_tr_fail[i], p_YZcYZac)
  n_YcYbc[] <- n_Yc_tr_fail[i] - n_YcYac[i]
  
  n_ZaZab[] <- rbinom(n_Za_tr_fail[i], p_YZaYZab)
  n_ZaZac[] <- n_Za_tr_fail[i] - n_ZaZab[i]
  n_ZbZab[] <- rbinom(n_Zb_tr_fail[i], p_YZbYZab)
  n_ZbZbc[] <- n_Zb_tr_fail[i] - n_ZbZab[i]
  n_ZcZac[] <- rbinom(n_Zc_tr_fail[i], p_YZcYZac)
  n_ZcZbc[] <- n_Zc_tr_fail[i] - n_ZcZac[i]
  
  ##now move to the double resistance compartments
  #first calculate total outflow -- calculate exit probabilities and then draw
  p_YabOut <- 1 - exp(-(E_c*T_s + (1 - E_c)*kappa*T_sr + d))
  p_YacOut <- 1 - exp(-(E_b*T_s + (1 - E_b)*kappa*T_sr + d))
  p_YbcOut <- 1 - exp(-(E_a*T_s + (1 - E_a)*kappa*T_sr + d))
  p_ZabOut <- 1 - exp(-(E_c*T_m + d))
  p_ZacOut <- 1 - exp(-(E_b*T_m + d))
  p_ZbcOut <- 1 - exp(-(E_a*T_m + d))
  
  n_YabOut[] <- rbinom(Y_ab[i], p_YabOut)
  n_YacOut[] <- rbinom(Y_ac[i], p_YacOut)
  n_YbcOut[] <- rbinom(Y_bc[i], p_YbcOut)
  n_ZabOut[] <- rbinom(Z_ab[i], p_ZabOut)
  n_ZacOut[] <- rbinom(Z_ac[i], p_ZacOut)
  n_ZbcOut[] <- rbinom(Z_bc[i], p_ZbcOut)
  
  #because these are already dual resistance, the only ones leaving and not recovering are developing resistance to the same drug, so separate those and the rest recover
  p_Yab_tr_fail <- (E_c*w_c*T_s)/(E_c*T_s + (1 - E_c)*kappa*T_sr + d)
  p_Yac_tr_fail <- (E_b*w_b*T_s)/(E_b*T_s + (1 - E_b)*kappa*T_sr + d)
  p_Ybc_tr_fail <- (E_a*w_a*T_s)/(E_a*T_s + (1 - E_a)*kappa*T_sr + d)
  p_Zab_tr_fail <- (E_c*w_c*T_m)/(E_c*T_m + d)
  p_Zac_tr_fail <- (E_b*w_b*T_m)/(E_b*T_m + d)
  p_Zbc_tr_fail <- (E_a*w_a*T_m)/(E_a*T_m + d)
  
  #split out resistant flows to Yabc vs recoveries to S
  n_YabYabc[] <- rbinom(n_YabOut[i], p_Yab_tr_fail)
  n_YabS[] <- n_YabOut[i] - n_YabYabc[i]
  n_YacYabc[] <- rbinom(n_YacOut[i], p_Yac_tr_fail)
  n_YacS[] <- n_YacOut[i] - n_YacYabc[i]
  n_YbcYabc[] <- rbinom(n_YbcOut[i], p_Ybc_tr_fail)
  n_YbcS[] <- n_YbcOut[i] - n_YbcYabc[i]
  
  n_ZabZabc[] <- rbinom(n_ZabOut[i], p_Zab_tr_fail)
  n_ZabS[] <- n_ZabOut[i] - n_ZabZabc[i]
  n_ZacZabc[] <- rbinom(n_ZacOut[i], p_Zac_tr_fail)
  n_ZacS[] <- n_ZacOut[i] - n_ZacZabc[i]
  n_ZbcZabc[] <- rbinom(n_ZbcOut[i], p_Zbc_tr_fail)
  n_ZbcS[] <- n_ZbcOut[i] - n_ZbcZabc[i]
  
  #recovery from triple resistance states
  p_YabcOut <- 1 - exp(-(kappa*T_sr + d))
  p_ZabcOut <- 1 - exp(-d)
  n_YabcOut[] <- rbinom(Y_abc[i], p_YabcOut)
  n_ZabcOut[] <- rbinom(Z_abc[i], p_ZabcOut)
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #make indicator variables for which ordering we're in 
  all_off <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_A <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_B <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  only_C <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1)
  A_B_on <- (lock_threshold_A >= 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  update(C_remain_lock) <- C_remain_lock + A_B_on
  A_C_on <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1) 
  update(B_remain_lock) <- B_remain_lock + A_C_on
  B_C_on <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C >= 1) 
  update(A_remain_lock) <- A_remain_lock + B_C_on
  
  E_a <- 1/3*all_off + 0.5*(only_B + only_C) + max((B_C_on), (A_remain_lock >= 1))
  E_b <- 1/3*all_off + 0.5*(only_A + only_C) + max((A_C_on),(B_remain_lock >= 1))
  E_c <- 1/3*all_off + 0.5*(only_A + only_B) + max((A_B_on),(C_remain_lock >= 1))
  
  sum_tx_probs = E_a + E_c + E_b 
  
  #sum array variables for book-keeping / accounting outputs
  sum_n_SOut <- sum(n_SOut)
  sum_n_SYZa <- sum(n_SYZa)
  sum_n_SYZb <- sum(n_SYZb)
  sum_n_SYZc <- sum(n_SYZc)
  sum_n_SYZab <- sum(n_SYZab)
  sum_n_SYZac <- sum(n_SYZac)
  sum_n_SYZbc <- sum(n_SYZbc)
  sum_n_SYZabc <- sum(n_SYZabc)
  sum_n_SYZs <- sum(n_SYZs)
  sum_n_SYs <- sum(n_SYs)
  sum_n_SZs <- sum(n_SZs)
  sum_n_SYa <- sum(n_SYa)
  sum_n_SZa <- sum(n_SZa)
  sum_n_SYb <- sum(n_SYb)
  sum_n_SZb <- sum(n_SZb)
  sum_n_SYc <- sum(n_SYc)
  sum_n_SZc <- sum(n_SZc)
  sum_n_SYab <- sum(n_SYab)
  sum_n_SZab <- sum(n_SZab)
  sum_n_SYac <- sum(n_SYac)
  sum_n_SZac <- sum(n_SZac)
  sum_n_SYbc <- sum(n_SYbc)
  sum_n_SZbc <- sum(n_SZbc)
  sum_n_SYabc <- sum(n_SYabc)
  sum_n_SZabc <- sum(n_SZabc)
  sum_n_YsOut <- sum(n_YsOut)
  sum_n_ZsOut <- sum(n_ZsOut)
  sum_n_Ys_tr <- sum(n_Ys_tr)
  sum_n_Ys_nr <- sum(n_Ys_nr)
  sum_n_Zs_tr <- sum(n_Zs_tr)
  sum_n_Zs_nr <- sum(n_Zs_nr)
  sum_n_Ys_tr_fail <- sum(n_Ys_tr_fail)
  sum_n_Ys_tr_success <- sum(n_Ys_tr_success)
  sum_n_Zs_tr_fail <- sum(n_Zs_tr_fail)
  sum_n_Zs_tr_success <- sum(n_Zs_tr_success)
  sum_n_YsS <- sum(n_YsS)
  sum_n_ZsS <- sum(n_ZsS)
  sum_n_YsYa <- sum(n_YsYa)
  sum_n_YsYb <- sum(n_YsYb)
  sum_n_YsYc <- sum(n_YsYc)
  sum_n_ZsZa <- sum(n_ZsZa)
  sum_n_ZsZb <- sum(n_ZsZb)
  sum_n_ZsZc <- sum(n_ZsZc)
  sum_n_YaOut <- sum(n_YaOut)
  sum_n_YbOut <- sum(n_YbOut)
  sum_n_YcOut <- sum(n_YcOut)
  sum_n_ZaOut <- sum(n_ZaOut)
  sum_n_ZbOut <- sum(n_ZbOut)
  sum_n_ZcOut <- sum(n_ZcOut)
  sum_n_Ya_nr <- sum(n_Ya_nr)
  sum_n_Ya_tr <- sum(n_Ya_tr)
  sum_n_Yb_nr <- sum(n_Yb_nr)
  sum_n_Yb_tr <- sum(n_Yb_tr)
  sum_n_Yc_nr <- sum(n_Yc_nr)
  sum_n_Yc_tr <- sum(n_Yc_tr)
  sum_n_Za_nr <- sum(n_Za_nr)
  sum_n_Za_tr <- sum(n_Za_tr)
  sum_n_Zb_nr <- sum(n_Zb_nr)
  sum_n_Zb_tr <- sum(n_Zb_tr)
  sum_n_Zc_nr <- sum(n_Zc_nr)
  sum_n_Zc_tr <- sum(n_Zc_tr)
  sum_n_Ya_tr_fail <- sum(n_Ya_tr_fail)
  sum_n_Yb_tr_fail <- sum(n_Yb_tr_fail)
  sum_n_Yc_tr_fail <- sum(n_Yc_tr_fail)
  sum_n_Za_tr_fail <- sum(n_Za_tr_fail)
  sum_n_Zb_tr_fail <- sum(n_Zb_tr_fail)
  sum_n_Zc_tr_fail <- sum(n_Zc_tr_fail)
  sum_n_Ya_tr_success <- sum(n_Ya_tr_success)
  sum_n_Yb_tr_success <- sum(n_Yb_tr_success)
  sum_n_Yc_tr_success <- sum(n_Yc_tr_success)
  sum_n_Za_tr_success <- sum(n_Za_tr_success)
  sum_n_Zb_tr_success <- sum(n_Zb_tr_success)
  sum_n_Zc_tr_success <- sum(n_Zc_tr_success)
  sum_n_YaS <- sum(n_YaS)
  sum_n_YbS <- sum(n_YbS)
  sum_n_YcS <- sum(n_YcS)
  sum_n_ZaS <- sum(n_ZaS)
  sum_n_ZbS <- sum(n_ZbS)
  sum_n_ZcS <- sum(n_ZcS)
  sum_n_YaYab <- sum(n_YaYab)
  sum_n_YaYac <- sum(n_YaYac)
  sum_n_YbYab <- sum(n_YbYab)
  sum_n_YbYbc <- sum(n_YbYbc)
  sum_n_YcYac <- sum(n_YcYac)
  sum_n_YcYbc <- sum(n_YcYbc)
  sum_n_ZaZab <- sum(n_ZaZab)
  sum_n_ZaZac <- sum(n_ZaZac)
  sum_n_ZbZab <- sum(n_ZbZab)
  sum_n_ZbZbc <- sum(n_ZbZbc)
  sum_n_ZcZac <- sum(n_ZcZac)
  sum_n_ZcZbc <- sum(n_ZcZbc)
  sum_n_YabOut <- sum(n_YabOut)
  sum_n_YacOut <- sum(n_YacOut)
  sum_n_YbcOut <- sum(n_YbcOut)
  sum_n_ZabOut <- sum(n_ZabOut)
  sum_n_ZacOut <- sum(n_ZacOut)
  sum_n_ZbcOut <- sum(n_ZbcOut)
  sum_n_YabYabc <- sum(n_YabYabc)
  sum_n_YabS <- sum(n_YabS)
  sum_n_YacYabc <- sum(n_YacYabc)
  sum_n_YacS <- sum(n_YacS)
  sum_n_YbcYabc <- sum(n_YbcYabc)
  sum_n_YbcS<- sum(n_YbcS)
  sum_n_ZabZabc <- sum(n_ZabZabc)
  sum_n_ZabS <- sum(n_ZabS)
  sum_n_ZacZabc <- sum(n_ZacZabc)
  sum_n_ZacS <- sum(n_ZacS)
  sum_n_ZbcZabc <- sum(n_ZbcZabc)
  sum_n_ZbcS<- sum(n_ZbcS)
  sum_n_YabcOut <- sum(n_YabcOut)
  sum_n_ZabcOut <- sum(n_ZabcOut)
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  initial(A_remain_lock) <- 0
  initial(B_remain_lock) <- 0
  initial(C_remain_lock) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  dim(n_SOut) <- N_risk
  dim(p_SOut) <- N_risk
  dim(p_SYZa) <- N_risk
  dim(p_SYZb) <- N_risk
  dim(p_SYZc) <- N_risk
  dim(p_SYZab) <- N_risk
  dim(p_SYZac) <- N_risk  
  dim(p_SYZbc) <- N_risk
  dim(p_SYZabc) <- N_risk
  dim(n_SYZa) <- N_risk
  dim(n_SYZb) <- N_risk
  dim(n_SYZc) <- N_risk
  dim(n_SYZab) <- N_risk
  dim(n_SYZac) <- N_risk
  dim(n_SYZbc) <- N_risk
  dim(n_SYZabc) <- N_risk
  dim(n_SYZs) <- N_risk
  dim(n_SYs) <- N_risk
  dim(n_SZs) <- N_risk
  dim(n_SYa) <- N_risk
  dim(n_SZa) <- N_risk
  dim(n_SYb) <- N_risk
  dim(n_SZb) <- N_risk
  dim(n_SYc) <- N_risk
  dim(n_SZc) <- N_risk
  dim(n_SYab) <- N_risk
  dim(n_SZab) <- N_risk
  dim(n_SYac) <- N_risk
  dim(n_SZac) <- N_risk
  dim(n_SYbc) <- N_risk
  dim(n_SZbc) <- N_risk
  dim(n_SYabc) <- N_risk
  dim(n_SZabc) <- N_risk
  dim(n_YsOut) <- N_risk
  dim(n_ZsOut) <- N_risk
  dim(n_Ys_tr) <- N_risk
  dim(n_Ys_nr) <- N_risk
  dim(n_Zs_tr) <- N_risk
  dim(n_Zs_nr) <- N_risk
  dim(n_Ys_tr_fail) <- N_risk
  dim(n_Ys_tr_success) <- N_risk
  dim(n_Zs_tr_fail) <- N_risk
  dim(n_Zs_tr_success) <- N_risk
  dim(n_YsS) <- N_risk
  dim(n_ZsS) <- N_risk
  dim(n_YsYa) <- N_risk
  dim(n_YsYb) <- N_risk
  dim(n_YsYc) <- N_risk
  dim(n_ZsZa) <- N_risk
  dim(n_ZsZb) <- N_risk
  dim(n_ZsZc) <- N_risk
  dim(n_YaOut) <- N_risk
  dim(n_YbOut) <- N_risk
  dim(n_YcOut) <- N_risk
  dim(n_ZaOut) <- N_risk
  dim(n_ZbOut) <- N_risk
  dim(n_ZcOut) <- N_risk
  dim(n_Ya_nr) <- N_risk
  dim(n_Ya_tr) <- N_risk
  dim(n_Yb_nr) <- N_risk
  dim(n_Yb_tr) <- N_risk
  dim(n_Yc_nr) <- N_risk
  dim(n_Yc_tr) <- N_risk
  dim(n_Za_nr) <- N_risk
  dim(n_Za_tr) <- N_risk
  dim(n_Zb_nr) <- N_risk
  dim(n_Zb_tr) <- N_risk
  dim(n_Zc_nr) <- N_risk
  dim(n_Zc_tr) <- N_risk
  dim(n_Ya_tr_fail) <- N_risk
  dim(n_Yb_tr_fail) <- N_risk
  dim(n_Yc_tr_fail) <- N_risk
  dim(n_Za_tr_fail) <- N_risk
  dim(n_Zb_tr_fail) <- N_risk
  dim(n_Zc_tr_fail) <- N_risk
  dim(n_Ya_tr_success) <- N_risk
  dim(n_Yb_tr_success) <- N_risk
  dim(n_Yc_tr_success) <- N_risk
  dim(n_Za_tr_success) <- N_risk
  dim(n_Zb_tr_success) <- N_risk
  dim(n_Zc_tr_success) <- N_risk
  dim(n_YaS) <- N_risk
  dim(n_YbS) <- N_risk
  dim(n_YcS) <- N_risk
  dim(n_ZaS) <- N_risk
  dim(n_ZbS) <- N_risk
  dim(n_ZcS) <- N_risk
  dim(n_YaYab) <- N_risk
  dim(n_YaYac) <- N_risk
  dim(n_YbYab) <- N_risk
  dim(n_YbYbc) <- N_risk
  dim(n_YcYac) <- N_risk
  dim(n_YcYbc) <- N_risk
  dim(n_ZaZab) <- N_risk
  dim(n_ZaZac) <- N_risk
  dim(n_ZbZab) <- N_risk
  dim(n_ZbZbc) <- N_risk
  dim(n_ZcZac) <- N_risk
  dim(n_ZcZbc) <- N_risk
  dim(n_YabOut) <- N_risk
  dim(n_YacOut) <- N_risk
  dim(n_YbcOut) <- N_risk
  dim(n_ZabOut) <- N_risk
  dim(n_ZacOut) <- N_risk
  dim(n_ZbcOut) <- N_risk
  dim(n_YabYabc) <- N_risk
  dim(n_YabS) <- N_risk
  dim(n_YacYabc) <- N_risk
  dim(n_YacS) <- N_risk
  dim(n_YbcYabc) <- N_risk
  dim(n_YbcS) <- N_risk
  dim(n_ZabZabc) <- N_risk
  dim(n_ZabS) <- N_risk
  dim(n_ZacZabc) <- N_risk
  dim(n_ZacS) <- N_risk
  dim(n_ZbcZabc) <- N_risk
  dim(n_ZbcS) <- N_risk
  dim(n_YabcOut) <- N_risk
  dim(n_ZabcOut) <- N_risk
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev
  output(sum_tx_probs) <- sum_tx_probs
  output(p_SOut) <- p_SOut
  output(sum_n_SOut) <- sum_n_SOut
  output(p_SYZa) <- p_SYZa
  output(p_SYZb) <- p_SYZb
  output(p_SYZc) <- p_SYZc
  output(p_SYZab) <- p_SYZab
  output(p_SYZac) <- p_SYZac
  output(p_SYZbc) <- p_SYZbc
  output(p_SYZabc) <- p_SYZabc
  output(sum_n_SYZa) <- sum_n_SYZa
  output(sum_n_SYZb) <- sum_n_SYZb
  output(sum_n_SYZc) <- sum_n_SYZc
  output(sum_n_SYZab) <- sum_n_SYZab
  output(sum_n_SYZac) <- sum_n_SYZac
  output(sum_n_SYZbc) <- sum_n_SYZbc
  output(sum_n_SYZabc) <- sum_n_SYZabc
  output(sum_n_SYZs) <- sum_n_SYZs
  output(sum_n_SYs) <- sum_n_SYs
  output(sum_n_SZs) <- sum_n_SZs
  output(sum_n_SYa) <- sum_n_SYa
  output(sum_n_SZa) <- sum_n_SZa
  output(sum_n_SYb) <- sum_n_SYb
  output(sum_n_SZb) <- sum_n_SZb
  output(sum_n_SYc) <- sum_n_SYc
  output(sum_n_SZc) <- sum_n_SZc
  output(sum_n_SYab) <- sum_n_SYab
  output(sum_n_SZab) <- sum_n_SZab
  output(sum_n_SYac) <- sum_n_SYac
  output(sum_n_SZac) <- sum_n_SZac
  output(sum_n_SYbc) <- sum_n_SYbc
  output(sum_n_SZbc) <- sum_n_SZbc
  output(sum_n_SYabc) <- sum_n_SYabc
  output(sum_n_SZabc) <- sum_n_SZabc
  output(p_YsOut) <- p_YsOut
  output(p_ZsOut) <- p_ZsOut
  output(p_Ys_tr) <- p_Ys_tr
  output(p_Zs_tr) <- p_Zs_tr
  output(p_YZs_tr_fail) <- p_YZs_tr_fail
  output(p_YZs_YZa) <- p_YZs_YZa
  output(p_YZs_YZb) <- p_YZs_YZb
  output(p_YZs_YZc) <- p_YZs_YZc
  output(p_YaOut) <- p_YaOut
  output(p_YbOut) <- p_YbOut
  output(p_YcOut) <- p_YcOut
  output(p_ZaOut) <- p_ZaOut
  output(p_ZbOut) <- p_ZbOut
  output(p_ZcOut) <- p_ZcOut
  output(p_Ya_nr) <- p_Ya_nr
  output(p_Yb_nr) <- p_Yb_nr
  output(p_Yc_nr) <- p_Yc_nr
  output(p_Za_nr) <- p_Za_nr
  output(p_Zb_nr) <- p_Zb_nr
  output(p_Zc_nr) <- p_Zc_nr
  output(p_Ya_tr_fail) <- p_Ya_tr_fail
  output(p_Yb_tr_fail) <- p_Yb_tr_fail
  output(p_Yc_tr_fail) <- p_Yc_tr_fail
  output(p_Za_tr_fail) <- p_Za_tr_fail
  output(p_Zb_tr_fail) <- p_Zb_tr_fail
  output(p_Zc_tr_fail) <- p_Zc_tr_fail
  output(p_YZaYZab) <- p_YZaYZab
  output(p_YZbYZab) <- p_YZbYZab
  output(p_YZcYZac) <- p_YZcYZac
  output(p_YabOut) <- p_YabOut 
  output(p_YacOut) <- p_YacOut 
  output(p_YbcOut) <- p_YbcOut 
  output(p_ZabOut) <- p_ZabOut 
  output(p_ZacOut) <- p_ZacOut 
  output(p_ZbcOut) <- p_ZbcOut
  output(p_Yab_tr_fail) <- p_Yab_tr_fail
  output(p_Yac_tr_fail) <- p_Yac_tr_fail
  output(p_Ybc_tr_fail) <- p_Ybc_tr_fail
  output(p_Zab_tr_fail) <- p_Zab_tr_fail
  output(p_Zac_tr_fail) <- p_Zac_tr_fail
  output(p_Zbc_tr_fail) <- p_Zbc_tr_fail
  output(sum_n_YsOut) <- sum_n_YsOut
  output(sum_n_ZsOut) <- sum_n_ZsOut
  output(sum_n_Ys_tr) <- sum_n_Ys_tr
  output(sum_n_Ys_nr) <- sum_n_Ys_nr
  output(sum_n_Zs_tr) <- sum_n_Zs_tr
  output(sum_n_Zs_nr) <- sum_n_Zs_nr
  output(sum_n_Ys_tr_fail) <- sum_n_Ys_tr_fail
  output(sum_n_Ys_tr_success) <- sum_n_Ys_tr_success
  output(sum_n_Zs_tr_fail) <- sum_n_Zs_tr_fail
  output(sum_n_Zs_tr_success) <- sum_n_Zs_tr_success
  output(sum_n_YsS) <- sum_n_YsS
  output(sum_n_ZsS) <- sum_n_ZsS
  output(sum_n_YsYa) <- sum_n_YsYa
  output(sum_n_YsYb) <- sum_n_YsYb
  output(sum_n_YsYc) <- sum_n_YsYc
  output(sum_n_ZsZa) <- sum_n_ZsZa
  output(sum_n_ZsZb) <- sum_n_ZsZb
  output(sum_n_ZsZc) <- sum_n_ZsZc
  output(sum_n_YaOut) <- sum_n_YaOut
  output(sum_n_YbOut) <- sum_n_YbOut
  output(sum_n_YcOut) <- sum_n_YcOut
  output(sum_n_ZaOut) <- sum_n_ZaOut
  output(sum_n_ZbOut) <- sum_n_ZbOut
  output(sum_n_ZcOut) <- sum_n_ZcOut
  output(sum_n_Ya_nr) <- sum_n_Ya_nr
  output(sum_n_Ya_tr) <- sum_n_Ya_tr
  output(sum_n_Yb_nr) <- sum_n_Yb_nr
  output(sum_n_Yb_tr) <- sum_n_Yb_tr
  output(sum_n_Yc_nr) <- sum_n_Yc_nr
  output(sum_n_Yc_tr) <- sum_n_Yc_tr
  output(sum_n_Za_nr) <- sum_n_Za_nr
  output(sum_n_Za_tr) <- sum_n_Za_tr
  output(sum_n_Zb_nr) <- sum_n_Zb_nr
  output(sum_n_Zb_tr) <- sum_n_Zb_tr
  output(sum_n_Zc_nr) <- sum_n_Zc_nr
  output(sum_n_Zc_tr) <- sum_n_Zc_tr
  output(sum_n_Ya_tr_fail) <- sum_n_Ya_tr_fail
  output(sum_n_Yb_tr_fail) <- sum_n_Yb_tr_fail
  output(sum_n_Yc_tr_fail) <- sum_n_Yc_tr_fail
  output(sum_n_Za_tr_fail) <- sum_n_Za_tr_fail
  output(sum_n_Zb_tr_fail) <- sum_n_Zb_tr_fail
  output(sum_n_Zc_tr_fail) <- sum_n_Zc_tr_fail
  output(sum_n_Ya_tr_success) <- sum_n_Ya_tr_success
  output(sum_n_Yb_tr_success) <- sum_n_Yb_tr_success
  output(sum_n_Yc_tr_success) <- sum_n_Yc_tr_success
  output(sum_n_Za_tr_success) <- sum_n_Za_tr_success
  output(sum_n_Zb_tr_success) <- sum_n_Zb_tr_success
  output(sum_n_Zc_tr_success) <- sum_n_Zc_tr_success
  output(sum_n_YaS) <- sum_n_YaS
  output(sum_n_YbS) <- sum_n_YbS
  output(sum_n_YcS) <- sum_n_YcS
  output(sum_n_ZaS) <- sum_n_ZaS
  output(sum_n_ZbS) <- sum_n_ZbS
  output(sum_n_ZcS) <- sum_n_ZcS
  output(sum_n_YaYab) <- sum_n_YaYab
  output(sum_n_YaYac) <- sum_n_YaYac
  output(sum_n_YbYab) <- sum_n_YbYab
  output(sum_n_YbYbc) <- sum_n_YbYbc
  output(sum_n_YcYac) <- sum_n_YcYac 
  output(sum_n_YcYbc) <- sum_n_YcYbc 
  output(sum_n_ZaZab) <- sum_n_ZaZab
  output(sum_n_ZaZac) <- sum_n_ZaZac
  output(sum_n_ZbZab) <- sum_n_ZbZab
  output(sum_n_ZbZbc) <- sum_n_ZbZbc
  output(sum_n_ZcZac) <- sum_n_ZcZac 
  output(sum_n_ZcZbc) <- sum_n_ZcZbc 
  output(sum_n_YabOut) <- sum_n_YabOut 
  output(sum_n_YacOut) <- sum_n_YacOut 
  output(sum_n_YbcOut) <- sum_n_YbcOut
  output(sum_n_ZabOut) <- sum_n_ZabOut 
  output(sum_n_ZacOut) <- sum_n_ZacOut 
  output(sum_n_ZbcOut) <- sum_n_ZbcOut
  output(sum_n_YabYabc) <- sum_n_YabYabc
  output(sum_n_YabS) <- sum_n_YabS
  output(sum_n_YacYabc) <- sum_n_YacYabc
  output(sum_n_YacS) <- sum_n_YacS
  output(sum_n_YbcYabc) <- sum_n_YbcYabc
  output(sum_n_YbcS) <- sum_n_YbcS
  output(sum_n_ZabZabc) <- sum_n_ZabZabc
  output(sum_n_ZabS) <- sum_n_ZabS
  output(sum_n_ZacZabc) <- sum_n_ZacZabc
  output(sum_n_ZacS) <- sum_n_ZacS
  output(sum_n_ZbcZabc) <- sum_n_ZbcZabc
  output(sum_n_ZbcS) <- sum_n_ZbcS
  output(sum_n_ZabcOut) <- sum_n_ZabcOut 
  output(p_YabcOut) <- p_YabcOut
  output(p_ZabcOut) <- p_ZabcOut
  output(sum_n_YabcOut) <- sum_n_YabcOut
})

#fixing stickbreaking error
sequential_stochastic_6_5 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] - n_SOut[i] + n_YsS[i] + n_YaS[i] + n_YbS[i] + n_YcS[i] + n_YabS[i] + n_YacS[i] + n_YbcS[i] +
    n_YabcOut[i] + n_ZsS[i] + n_ZaS[i] + n_ZbS[i] + n_ZcS[i] + n_ZabS[i] + n_ZacS[i] + n_ZbcS[i] + n_ZabcOut[i] 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + n_SYs[i] - n_YsOut[i]
  update(Z_s[]) <- Z_s[i] + n_SZs[i] - n_ZsOut[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + n_SYa[i] + n_YsYa[i] - n_YaOut[i]
  update(Z_a[]) <- Z_a[i] + n_SZa[i] + n_ZsZa[i] - n_ZaOut[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + n_SYb[i] + n_YsYb[i] - n_YbOut[i]
  update(Z_b[]) <- Z_b[i] + n_SZb[i] + n_ZsZb[i] - n_ZbOut[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + n_SYc[i] + n_YsYc[i] - n_YcOut[i] 
  update(Z_c[]) <- Z_c[i] + n_SZc[i] + n_ZsZc[i] - n_ZcOut[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + n_SYab[i] + n_YaYab[i] + n_YbYab[i] - n_YabOut[i]
  update(Z_ab[]) <- Z_ab[i] + n_SZab[i] + n_ZaZab[i] + n_ZbZab[i] - n_ZabOut[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + n_SYac[i] + n_YaYac[i] + n_YcYac[i] - n_YacOut[i]
  update(Z_ac[]) <- Z_ac[i] + n_SZac[i] + n_ZaZac[i] + n_ZcZac[i] - n_ZacOut[i] 
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + n_SYbc[i] + n_YbYbc[i] + n_YcYbc[i] - n_YbcOut[i]
  update(Z_bc[]) <- Z_bc[i] + n_SZbc[i] + n_ZbZbc[i] + n_ZcZbc[i] - n_ZbcOut[i]
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + n_SYabc[i] + n_YabYabc[i] + n_YacYabc[i] + n_YbcYabc[i] - n_YabcOut[i]
  update(Z_abc[]) <- Z_abc[i] + n_SZabc[i] + n_ZabZabc[i] + n_ZacZabc[i] + n_ZbcZabc[i] - n_ZabcOut[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  ##STOCHASTIC calculations
  #calculate the probability of leaving S from the rate
  p_SOut[] <- 1 - exp(-sum_all_BI_terms[i])
  
  #calculate the number of people leaving S
  n_SOut[] <- rbinom(S[i], p_SOut[i])
  
  #divide probabilities leaving S 
  p_SYZa[] <- sum_BI_YZa[i]/sum_all_BI_terms[i]
  p_SYZb[] <- sum_BI_YZb[i]/sum_all_BI_terms[i]
  p_SYZc[] <- sum_BI_YZc[i]/sum_all_BI_terms[i]
  p_SYZab[] <- sum_BI_YZab[i]/sum_all_BI_terms[i]
  p_SYZac[] <- sum_BI_YZac[i]/sum_all_BI_terms[i]
  p_SYZbc[] <- sum_BI_YZbc[i]/sum_all_BI_terms[i]
  p_SYZabc[] <- sum_BI_YZabc[i]/sum_all_BI_terms[i] #fixed error here 6/3
  
  #divide the people leaving S into different flows
  n_SYZa[] <- rbinom(n_SOut[i], p_SYZa[i])
  n_SYZb[] <- rbinom(n_SOut[i] - n_SYZa[i], (p_SYZb[i])/(1 - p_SYZa[i] + 1e-8)) #i) subtract those already subtracted from n_SOut ii) renormalize the probability of exit
  n_SYZc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i], (p_SYZc[i])/(1 - p_SYZa[i] - p_SYZb[i] + 1e-8))
  n_SYZab[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i], (p_SYZab[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] + 1e-8))
  n_SYZac[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i], (p_SYZac[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] + 1e-8)) #fixed n to p issue here and remaining 2 eqs 6/5
  n_SYZbc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i], (p_SYZbc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] + 1e-8)) 
  n_SYZabc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i], (p_SYZabc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] - p_SYZbc[i] + 1e-8))
  n_SYZs[] <- n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] - n_SYZabc[i]
  
  #split into symptomatic vs asymptomatics
  n_SYs[] <- rbinom(n_SYZs[i], sigma)
  n_SZs[] <- n_SYZs[i] - n_SYs[i] #fixed issue here
  n_SYa[] <- rbinom(n_SYZa[i], sigma)
  n_SZa[] <- n_SYZa[i] - n_SYa[i]
  n_SYb[] <- rbinom(n_SYZb[i], sigma)
  n_SZb[] <- n_SYZb[i] - n_SYb[i]
  n_SYc[] <- rbinom(n_SYZc[i], sigma)
  n_SZc[] <- n_SYZc[i] - n_SYc[i]
  n_SYab[] <- rbinom(n_SYZab[i], sigma)
  n_SZab[] <- n_SYZab[i] - n_SYab[i]
  n_SYac[] <- rbinom(n_SYZac[i], sigma)
  n_SZac[] <- n_SYZac[i] - n_SYac[i]
  n_SYbc[] <- rbinom(n_SYZbc[i], sigma)
  n_SZbc[] <- n_SYZbc[i] - n_SYbc[i]
  n_SYabc[] <- rbinom(n_SYZabc[i], sigma)
  n_SZabc[] <- n_SYZabc[i] - n_SYabc[i]
  
  #calculate the probability of leaving Ys and Zs by converting exit rates
  p_YsOut <- 1 - exp(-(T_s + d))
  p_ZsOut <- 1 - exp(-(T_m + d))
  
  #calculate the number of people leaving Ys and Zs based on the probabilities
  n_YsOut[] <- rbinom(Y_s[i], p_YsOut)
  n_ZsOut[] <- rbinom(Z_s[i], p_ZsOut)
  
  #split the people recovering naturally (nr) vs the people leaving due to treatment (tr)
  p_Ys_tr <- T_s / (T_s + d)
  n_Ys_tr[] <- rbinom(n_YsOut[i], p_Ys_tr)
  n_Ys_nr[] <- n_YsOut[i] - n_Ys_tr[i]
  p_Zs_tr <- T_m / (T_m + d)
  n_Zs_tr[] <- rbinom(n_ZsOut[i], p_Zs_tr)
  n_Zs_nr[] <- n_ZsOut[i] - n_Zs_tr[i]
  
  #of the people being treated, split the flows into recovery vs. resistance
  p_YZs_tr_fail <- E_a*w_a + E_b*w_b + E_c*w_c
  n_Ys_tr_fail[] <- rbinom(n_Ys_tr[i], p_YZs_tr_fail)
  n_Ys_tr_success[] <- n_Ys_tr[i] - n_Ys_tr_fail[i]
  n_Zs_tr_fail[] <- rbinom(n_Zs_tr[i], p_YZs_tr_fail)
  n_Zs_tr_success[] <- n_Zs_tr[i] - n_Zs_tr_fail[i]
  
  #combine those recovering from successful treatment and natural recovery as leaving Ys and Zs and going back to S
  n_YsS[] <- n_Ys_tr_success[i] + n_Ys_nr[i]
  n_ZsS[] <- n_Zs_tr_success[i] + n_Zs_nr[i]
  
  #split those who are treated and develop resistance into the three possible flows
  #first take the individually normalized probabilities of each flow
  p_YZs_YZa <- E_a*w_a/p_YZs_tr_fail
  p_YZs_YZb <- E_b*w_b/p_YZs_tr_fail
  p_YZs_YZc <- E_c*w_c/p_YZs_tr_fail
  
  #do stick-breaking with Y
  n_YsYa[] <- rbinom(n_Ys_tr_fail[i], p_YZs_YZa)
  n_YsYb[] <- rbinom(n_Ys_tr_fail[i] - n_YsYa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_YsYc[] <- n_Ys_tr_fail[i] - n_YsYa[i] - n_YsYb[i]
  
  #do stick-breaking with Z
  n_ZsZa[] <- rbinom(n_Zs_tr_fail[i], p_YZs_YZa)
  n_ZsZb[] <- rbinom(n_Zs_tr_fail[i] - n_ZsZa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_ZsZc[] <- n_Zs_tr_fail[i] - n_ZsZa[i] - n_ZsZb[i]
  
  ##now onto single resistance compartments. Probabilities of exit will be different for Y and Z compartments due to retreatment
  #probabilities for exit flows from Ya, Yb, Yc, Za, Zb, and Zc
  p_YaOut <- 1 - exp(-((E_b + E_c)*T_s + E_a*kappa*T_sr + d))
  p_YbOut <- 1 - exp(-((E_a + E_c)*T_s + E_b*kappa*T_sr + d))
  p_YcOut <- 1 - exp(-((E_a + E_b)*T_s + E_c*kappa*T_sr + d))
  p_ZaOut <- 1 - exp(-((E_b + E_c)*T_m + d))
  p_ZbOut <- 1 - exp(-((E_a + E_c)*T_m + d))
  p_ZcOut <- 1 - exp(-((E_a + E_b)*T_m + d))
  
  #calculate those leaving the single resistance compartments
  n_YaOut[] <- rbinom(Y_a[i], p_YaOut)
  n_YbOut[] <- rbinom(Y_b[i], p_YbOut)
  n_YcOut[] <- rbinom(Y_c[i], p_YcOut)
  n_ZaOut[] <- rbinom(Z_a[i], p_ZaOut)
  n_ZbOut[] <- rbinom(Z_b[i], p_ZbOut)
  n_ZcOut[] <- rbinom(Z_c[i], p_ZcOut)
  
  #calculate the probabilities of leaving single resistance compartments due to natural recovery
  p_Ya_nr <- d/((E_b + E_c)*T_s + E_a*kappa*T_sr + d)
  p_Yb_nr <- d/((E_a + E_c)*T_s + E_b*kappa*T_sr + d)
  p_Yc_nr <- d/((E_a + E_b)*T_s + E_c*kappa*T_sr + d)
  p_Za_nr <- d/((E_b + E_c)*T_m + d)
  p_Zb_nr <- d/((E_a + E_c)*T_m + d)
  p_Zc_nr <- d/((E_a + E_b)*T_m + d)
  
  #split the people leaving single resistance into natural recovery vs treatment
  n_Ya_nr[] <- rbinom(n_YaOut[i], p_Ya_nr)
  n_Ya_tr[] <- n_YaOut[i] - n_Ya_nr[i]
  n_Yb_nr[] <- rbinom(n_YbOut[i], p_Yb_nr)
  n_Yb_tr[] <- n_YbOut[i] - n_Yb_nr[i]
  n_Yc_nr[] <- rbinom(n_YcOut[i], p_Yc_nr)
  n_Yc_tr[] <- n_YcOut[i] - n_Yc_nr[i]
  n_Za_nr[] <- rbinom(n_ZaOut[i], p_Za_nr)
  n_Za_tr[] <- n_ZaOut[i] - n_Za_nr[i]
  n_Zb_nr[] <- rbinom(n_ZbOut[i], p_Zb_nr)
  n_Zb_tr[] <- n_ZbOut[i] - n_Zb_nr[i]
  n_Zc_nr[] <- rbinom(n_ZcOut[i], p_Zc_nr)
  n_Zc_tr[] <- n_ZcOut[i] - n_Zc_nr[i]
  
  #split those leaving single resistance after treatment into those who fail treatment vs. those who recover
  #note that here, treatment failure means people who develop resistance on treatment, not those treated with the drug they are resistance bc they stay in the resistant compartment
  p_Ya_tr_fail <- (E_b*w_b*T_s + E_c*w_c*T_s)/((E_b + E_c)*T_s + E_a*kappa*T_sr)
  p_Yb_tr_fail <- (E_a*w_a*T_s + E_c*w_c*T_s)/((E_a + E_c)*T_s + E_b*kappa*T_sr)
  p_Yc_tr_fail <- (E_a*w_a*T_s + E_b*w_b*T_s)/((E_a + E_b)*T_s + E_c*kappa*T_sr)
  p_Za_tr_fail <- (E_b*w_b + E_c*w_c)/ ((E_b + E_c) + 1e-8) #add very small term to denominator so this probability is always defined
  p_Zb_tr_fail <- (E_a*w_a + E_c*w_c)/ ((E_a + E_c) + 1e-8)
  p_Zc_tr_fail <- (E_a*w_a + E_b*w_b)/ ((E_a + E_b) + 1e-8)
  
  n_Ya_tr_fail[] <- rbinom(n_Ya_tr[i], p_Ya_tr_fail)
  n_Yb_tr_fail[] <- rbinom(n_Yb_tr[i], p_Yb_tr_fail)
  n_Yc_tr_fail[] <- rbinom(n_Yc_tr[i], p_Yc_tr_fail)
  n_Za_tr_fail[] <- rbinom(n_Za_tr[i], p_Za_tr_fail)
  n_Zb_tr_fail[] <- rbinom(n_Zb_tr[i], p_Zb_tr_fail)
  n_Zc_tr_fail[] <- rbinom(n_Zc_tr[i], p_Zc_tr_fail)
  
  n_Ya_tr_success[] <- n_Ya_tr[i] - n_Ya_tr_fail[i]
  n_Yb_tr_success[] <- n_Yb_tr[i] - n_Yb_tr_fail[i]
  n_Yc_tr_success[] <- n_Yc_tr[i] - n_Yc_tr_fail[i]
  n_Za_tr_success[] <- n_Za_tr[i] - n_Za_tr_fail[i]
  n_Zb_tr_success[] <- n_Zb_tr[i] - n_Zb_tr_fail[i]
  n_Zc_tr_success[] <- n_Zc_tr[i] - n_Zc_tr_fail[i]
  
  #combine treatment success and natural recover terms to get single resistance recovery flows
  n_YaS[] <- n_Ya_tr_success[i] + n_Ya_nr[i]
  n_YbS[] <- n_Yb_tr_success[i] + n_Yb_nr[i]
  n_YcS[] <- n_Yc_tr_success[i] + n_Yc_nr[i]
  n_ZaS[] <- n_Za_tr_success[i] + n_Za_nr[i]
  n_ZbS[] <- n_Zb_tr_success[i] + n_Zb_nr[i]
  n_ZcS[] <- n_Zc_tr_success[i] + n_Zc_nr[i]
  
  #now split the treatment failures into the two different resistance compartment flows
  p_YZaYZab <- E_b*w_b / (E_b*w_b + E_c*w_c + 1e-8) #add very small term to denominator so this probability is always defined
  p_YZbYZab <- E_a*w_a / (E_a*w_a + E_c*w_c + 1e-8)
  p_YZcYZac <- E_a*w_a / (E_a*w_a + E_b*w_b + 1e-8)
  
  n_YaYab[] <- rbinom(n_Ya_tr_fail[i], p_YZaYZab)
  n_YaYac[] <- n_Ya_tr_fail[i] - n_YaYab[i]
  n_YbYab[] <- rbinom(n_Yb_tr_fail[i], p_YZbYZab)
  n_YbYbc[] <- n_Yb_tr_fail[i] - n_YbYab[i]
  n_YcYac[] <- rbinom(n_Yc_tr_fail[i], p_YZcYZac)
  n_YcYbc[] <- n_Yc_tr_fail[i] - n_YcYac[i]
  
  n_ZaZab[] <- rbinom(n_Za_tr_fail[i], p_YZaYZab)
  n_ZaZac[] <- n_Za_tr_fail[i] - n_ZaZab[i]
  n_ZbZab[] <- rbinom(n_Zb_tr_fail[i], p_YZbYZab)
  n_ZbZbc[] <- n_Zb_tr_fail[i] - n_ZbZab[i]
  n_ZcZac[] <- rbinom(n_Zc_tr_fail[i], p_YZcYZac)
  n_ZcZbc[] <- n_Zc_tr_fail[i] - n_ZcZac[i]
  
  ##now move to the double resistance compartments
  #first calculate total outflow -- calculate exit probabilities and then draw
  p_YabOut <- 1 - exp(-(E_c*T_s + (1 - E_c)*kappa*T_sr + d))
  p_YacOut <- 1 - exp(-(E_b*T_s + (1 - E_b)*kappa*T_sr + d))
  p_YbcOut <- 1 - exp(-(E_a*T_s + (1 - E_a)*kappa*T_sr + d))
  p_ZabOut <- 1 - exp(-(E_c*T_m + d))
  p_ZacOut <- 1 - exp(-(E_b*T_m + d))
  p_ZbcOut <- 1 - exp(-(E_a*T_m + d))
  
  n_YabOut[] <- rbinom(Y_ab[i], p_YabOut)
  n_YacOut[] <- rbinom(Y_ac[i], p_YacOut)
  n_YbcOut[] <- rbinom(Y_bc[i], p_YbcOut)
  n_ZabOut[] <- rbinom(Z_ab[i], p_ZabOut)
  n_ZacOut[] <- rbinom(Z_ac[i], p_ZacOut)
  n_ZbcOut[] <- rbinom(Z_bc[i], p_ZbcOut)
  
  #because these are already dual resistance, the only ones leaving and not recovering are developing resistance to the same drug, so separate those and the rest recover
  p_Yab_tr_fail <- (E_c*w_c*T_s)/(E_c*T_s + (1 - E_c)*kappa*T_sr + d)
  p_Yac_tr_fail <- (E_b*w_b*T_s)/(E_b*T_s + (1 - E_b)*kappa*T_sr + d)
  p_Ybc_tr_fail <- (E_a*w_a*T_s)/(E_a*T_s + (1 - E_a)*kappa*T_sr + d)
  p_Zab_tr_fail <- (E_c*w_c*T_m)/(E_c*T_m + d)
  p_Zac_tr_fail <- (E_b*w_b*T_m)/(E_b*T_m + d)
  p_Zbc_tr_fail <- (E_a*w_a*T_m)/(E_a*T_m + d)
  
  #split out resistant flows to Yabc vs recoveries to S
  n_YabYabc[] <- rbinom(n_YabOut[i], p_Yab_tr_fail)
  n_YabS[] <- n_YabOut[i] - n_YabYabc[i]
  n_YacYabc[] <- rbinom(n_YacOut[i], p_Yac_tr_fail)
  n_YacS[] <- n_YacOut[i] - n_YacYabc[i]
  n_YbcYabc[] <- rbinom(n_YbcOut[i], p_Ybc_tr_fail)
  n_YbcS[] <- n_YbcOut[i] - n_YbcYabc[i]
  
  n_ZabZabc[] <- rbinom(n_ZabOut[i], p_Zab_tr_fail)
  n_ZabS[] <- n_ZabOut[i] - n_ZabZabc[i]
  n_ZacZabc[] <- rbinom(n_ZacOut[i], p_Zac_tr_fail)
  n_ZacS[] <- n_ZacOut[i] - n_ZacZabc[i]
  n_ZbcZabc[] <- rbinom(n_ZbcOut[i], p_Zbc_tr_fail)
  n_ZbcS[] <- n_ZbcOut[i] - n_ZbcZabc[i]
  
  #recovery from triple resistance states
  p_YabcOut <- 1 - exp(-(kappa*T_sr + d))
  p_ZabcOut <- 1 - exp(-d)
  n_YabcOut[] <- rbinom(Y_abc[i], p_YabcOut)
  n_ZabcOut[] <- rbinom(Z_abc[i], p_ZabcOut)
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #Binary switches between treatment
  #THIS PART IS UNIQUE TO THE MODEL STRATEGY
  E_a <- (lock_threshold_A < 1)
  E_c <- 1 - (lock_threshold_B < 1)
  E_b <- 1  - E_a - E_c
  
  sum_tx_probs = E_a + E_c + E_b 
  
  #sum array variables for book-keeping / accounting outputs
  sum_n_SOut <- sum(n_SOut)
  sum_n_SYZa <- sum(n_SYZa)
  sum_n_SYZb <- sum(n_SYZb)
  sum_n_SYZc <- sum(n_SYZc)
  sum_n_SYZab <- sum(n_SYZab)
  sum_n_SYZac <- sum(n_SYZac)
  sum_n_SYZbc <- sum(n_SYZbc)
  sum_n_SYZabc <- sum(n_SYZabc)
  sum_n_SYZs <- sum(n_SYZs)
  sum_n_SYs <- sum(n_SYs)
  sum_n_SZs <- sum(n_SZs)
  sum_n_SYa <- sum(n_SYa)
  sum_n_SZa <- sum(n_SZa)
  sum_n_SYb <- sum(n_SYb)
  sum_n_SZb <- sum(n_SZb)
  sum_n_SYc <- sum(n_SYc)
  sum_n_SZc <- sum(n_SZc)
  sum_n_SYab <- sum(n_SYab)
  sum_n_SZab <- sum(n_SZab)
  sum_n_SYac <- sum(n_SYac)
  sum_n_SZac <- sum(n_SZac)
  sum_n_SYbc <- sum(n_SYbc)
  sum_n_SZbc <- sum(n_SZbc)
  sum_n_SYabc <- sum(n_SYabc)
  sum_n_SZabc <- sum(n_SZabc)
  sum_n_YsOut <- sum(n_YsOut)
  sum_n_ZsOut <- sum(n_ZsOut)
  sum_n_Ys_tr <- sum(n_Ys_tr)
  sum_n_Ys_nr <- sum(n_Ys_nr)
  sum_n_Zs_tr <- sum(n_Zs_tr)
  sum_n_Zs_nr <- sum(n_Zs_nr)
  sum_n_Ys_tr_fail <- sum(n_Ys_tr_fail)
  sum_n_Ys_tr_success <- sum(n_Ys_tr_success)
  sum_n_Zs_tr_fail <- sum(n_Zs_tr_fail)
  sum_n_Zs_tr_success <- sum(n_Zs_tr_success)
  sum_n_YsS <- sum(n_YsS)
  sum_n_ZsS <- sum(n_ZsS)
  sum_n_YsYa <- sum(n_YsYa)
  sum_n_YsYb <- sum(n_YsYb)
  sum_n_YsYc <- sum(n_YsYc)
  sum_n_ZsZa <- sum(n_ZsZa)
  sum_n_ZsZb <- sum(n_ZsZb)
  sum_n_ZsZc <- sum(n_ZsZc)
  sum_n_YaOut <- sum(n_YaOut)
  sum_n_YbOut <- sum(n_YbOut)
  sum_n_YcOut <- sum(n_YcOut)
  sum_n_ZaOut <- sum(n_ZaOut)
  sum_n_ZbOut <- sum(n_ZbOut)
  sum_n_ZcOut <- sum(n_ZcOut)
  sum_n_Ya_nr <- sum(n_Ya_nr)
  sum_n_Ya_tr <- sum(n_Ya_tr)
  sum_n_Yb_nr <- sum(n_Yb_nr)
  sum_n_Yb_tr <- sum(n_Yb_tr)
  sum_n_Yc_nr <- sum(n_Yc_nr)
  sum_n_Yc_tr <- sum(n_Yc_tr)
  sum_n_Za_nr <- sum(n_Za_nr)
  sum_n_Za_tr <- sum(n_Za_tr)
  sum_n_Zb_nr <- sum(n_Zb_nr)
  sum_n_Zb_tr <- sum(n_Zb_tr)
  sum_n_Zc_nr <- sum(n_Zc_nr)
  sum_n_Zc_tr <- sum(n_Zc_tr)
  sum_n_Ya_tr_fail <- sum(n_Ya_tr_fail)
  sum_n_Yb_tr_fail <- sum(n_Yb_tr_fail)
  sum_n_Yc_tr_fail <- sum(n_Yc_tr_fail)
  sum_n_Za_tr_fail <- sum(n_Za_tr_fail)
  sum_n_Zb_tr_fail <- sum(n_Zb_tr_fail)
  sum_n_Zc_tr_fail <- sum(n_Zc_tr_fail)
  sum_n_Ya_tr_success <- sum(n_Ya_tr_success)
  sum_n_Yb_tr_success <- sum(n_Yb_tr_success)
  sum_n_Yc_tr_success <- sum(n_Yc_tr_success)
  sum_n_Za_tr_success <- sum(n_Za_tr_success)
  sum_n_Zb_tr_success <- sum(n_Zb_tr_success)
  sum_n_Zc_tr_success <- sum(n_Zc_tr_success)
  sum_n_YaS <- sum(n_YaS)
  sum_n_YbS <- sum(n_YbS)
  sum_n_YcS <- sum(n_YcS)
  sum_n_ZaS <- sum(n_ZaS)
  sum_n_ZbS <- sum(n_ZbS)
  sum_n_ZcS <- sum(n_ZcS)
  sum_n_YaYab <- sum(n_YaYab)
  sum_n_YaYac <- sum(n_YaYac)
  sum_n_YbYab <- sum(n_YbYab)
  sum_n_YbYbc <- sum(n_YbYbc)
  sum_n_YcYac <- sum(n_YcYac)
  sum_n_YcYbc <- sum(n_YcYbc)
  sum_n_ZaZab <- sum(n_ZaZab)
  sum_n_ZaZac <- sum(n_ZaZac)
  sum_n_ZbZab <- sum(n_ZbZab)
  sum_n_ZbZbc <- sum(n_ZbZbc)
  sum_n_ZcZac <- sum(n_ZcZac)
  sum_n_ZcZbc <- sum(n_ZcZbc)
  sum_n_YabOut <- sum(n_YabOut)
  sum_n_YacOut <- sum(n_YacOut)
  sum_n_YbcOut <- sum(n_YbcOut)
  sum_n_ZabOut <- sum(n_ZabOut)
  sum_n_ZacOut <- sum(n_ZacOut)
  sum_n_ZbcOut <- sum(n_ZbcOut)
  sum_n_YabYabc <- sum(n_YabYabc)
  sum_n_YabS <- sum(n_YabS)
  sum_n_YacYabc <- sum(n_YacYabc)
  sum_n_YacS <- sum(n_YacS)
  sum_n_YbcYabc <- sum(n_YbcYabc)
  sum_n_YbcS<- sum(n_YbcS)
  sum_n_ZabZabc <- sum(n_ZabZabc)
  sum_n_ZabS <- sum(n_ZabS)
  sum_n_ZacZabc <- sum(n_ZacZabc)
  sum_n_ZacS <- sum(n_ZacS)
  sum_n_ZbcZabc <- sum(n_ZbcZabc)
  sum_n_ZbcS<- sum(n_ZbcS)
  sum_n_YabcOut <- sum(n_YabcOut)
  sum_n_ZabcOut <- sum(n_ZabcOut)
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  dim(n_SOut) <- N_risk
  dim(p_SOut) <- N_risk
  dim(p_SYZa) <- N_risk
  dim(p_SYZb) <- N_risk
  dim(p_SYZc) <- N_risk
  dim(p_SYZab) <- N_risk
  dim(p_SYZac) <- N_risk  
  dim(p_SYZbc) <- N_risk
  dim(p_SYZabc) <- N_risk
  dim(n_SYZa) <- N_risk
  dim(n_SYZb) <- N_risk
  dim(n_SYZc) <- N_risk
  dim(n_SYZab) <- N_risk
  dim(n_SYZac) <- N_risk
  dim(n_SYZbc) <- N_risk
  dim(n_SYZabc) <- N_risk
  dim(n_SYZs) <- N_risk
  dim(n_SYs) <- N_risk
  dim(n_SZs) <- N_risk
  dim(n_SYa) <- N_risk
  dim(n_SZa) <- N_risk
  dim(n_SYb) <- N_risk
  dim(n_SZb) <- N_risk
  dim(n_SYc) <- N_risk
  dim(n_SZc) <- N_risk
  dim(n_SYab) <- N_risk
  dim(n_SZab) <- N_risk
  dim(n_SYac) <- N_risk
  dim(n_SZac) <- N_risk
  dim(n_SYbc) <- N_risk
  dim(n_SZbc) <- N_risk
  dim(n_SYabc) <- N_risk
  dim(n_SZabc) <- N_risk
  dim(n_YsOut) <- N_risk
  dim(n_ZsOut) <- N_risk
  dim(n_Ys_tr) <- N_risk
  dim(n_Ys_nr) <- N_risk
  dim(n_Zs_tr) <- N_risk
  dim(n_Zs_nr) <- N_risk
  dim(n_Ys_tr_fail) <- N_risk
  dim(n_Ys_tr_success) <- N_risk
  dim(n_Zs_tr_fail) <- N_risk
  dim(n_Zs_tr_success) <- N_risk
  dim(n_YsS) <- N_risk
  dim(n_ZsS) <- N_risk
  dim(n_YsYa) <- N_risk
  dim(n_YsYb) <- N_risk
  dim(n_YsYc) <- N_risk
  dim(n_ZsZa) <- N_risk
  dim(n_ZsZb) <- N_risk
  dim(n_ZsZc) <- N_risk
  dim(n_YaOut) <- N_risk
  dim(n_YbOut) <- N_risk
  dim(n_YcOut) <- N_risk
  dim(n_ZaOut) <- N_risk
  dim(n_ZbOut) <- N_risk
  dim(n_ZcOut) <- N_risk
  dim(n_Ya_nr) <- N_risk
  dim(n_Ya_tr) <- N_risk
  dim(n_Yb_nr) <- N_risk
  dim(n_Yb_tr) <- N_risk
  dim(n_Yc_nr) <- N_risk
  dim(n_Yc_tr) <- N_risk
  dim(n_Za_nr) <- N_risk
  dim(n_Za_tr) <- N_risk
  dim(n_Zb_nr) <- N_risk
  dim(n_Zb_tr) <- N_risk
  dim(n_Zc_nr) <- N_risk
  dim(n_Zc_tr) <- N_risk
  dim(n_Ya_tr_fail) <- N_risk
  dim(n_Yb_tr_fail) <- N_risk
  dim(n_Yc_tr_fail) <- N_risk
  dim(n_Za_tr_fail) <- N_risk
  dim(n_Zb_tr_fail) <- N_risk
  dim(n_Zc_tr_fail) <- N_risk
  dim(n_Ya_tr_success) <- N_risk
  dim(n_Yb_tr_success) <- N_risk
  dim(n_Yc_tr_success) <- N_risk
  dim(n_Za_tr_success) <- N_risk
  dim(n_Zb_tr_success) <- N_risk
  dim(n_Zc_tr_success) <- N_risk
  dim(n_YaS) <- N_risk
  dim(n_YbS) <- N_risk
  dim(n_YcS) <- N_risk
  dim(n_ZaS) <- N_risk
  dim(n_ZbS) <- N_risk
  dim(n_ZcS) <- N_risk
  dim(n_YaYab) <- N_risk
  dim(n_YaYac) <- N_risk
  dim(n_YbYab) <- N_risk
  dim(n_YbYbc) <- N_risk
  dim(n_YcYac) <- N_risk
  dim(n_YcYbc) <- N_risk
  dim(n_ZaZab) <- N_risk
  dim(n_ZaZac) <- N_risk
  dim(n_ZbZab) <- N_risk
  dim(n_ZbZbc) <- N_risk
  dim(n_ZcZac) <- N_risk
  dim(n_ZcZbc) <- N_risk
  dim(n_YabOut) <- N_risk
  dim(n_YacOut) <- N_risk
  dim(n_YbcOut) <- N_risk
  dim(n_ZabOut) <- N_risk
  dim(n_ZacOut) <- N_risk
  dim(n_ZbcOut) <- N_risk
  dim(n_YabYabc) <- N_risk
  dim(n_YabS) <- N_risk
  dim(n_YacYabc) <- N_risk
  dim(n_YacS) <- N_risk
  dim(n_YbcYabc) <- N_risk
  dim(n_YbcS) <- N_risk
  dim(n_ZabZabc) <- N_risk
  dim(n_ZabS) <- N_risk
  dim(n_ZacZabc) <- N_risk
  dim(n_ZacS) <- N_risk
  dim(n_ZbcZabc) <- N_risk
  dim(n_ZbcS) <- N_risk
  dim(n_YabcOut) <- N_risk
  dim(n_ZabcOut) <- N_risk
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev
  output(sum_tx_probs) <- sum_tx_probs
  output(p_SOut) <- p_SOut
  output(sum_n_SOut) <- sum_n_SOut
  output(p_SYZa) <- p_SYZa
  output(p_SYZb) <- p_SYZb
  output(p_SYZc) <- p_SYZc
  output(p_SYZab) <- p_SYZab
  output(p_SYZac) <- p_SYZac
  output(p_SYZbc) <- p_SYZbc
  output(p_SYZabc) <- p_SYZabc
  output(sum_n_SYZa) <- sum_n_SYZa
  output(sum_n_SYZb) <- sum_n_SYZb
  output(sum_n_SYZc) <- sum_n_SYZc
  output(sum_n_SYZab) <- sum_n_SYZab
  output(sum_n_SYZac) <- sum_n_SYZac
  output(sum_n_SYZbc) <- sum_n_SYZbc
  output(sum_n_SYZabc) <- sum_n_SYZabc
  output(sum_n_SYZs) <- sum_n_SYZs
  output(sum_n_SYs) <- sum_n_SYs
  output(sum_n_SZs) <- sum_n_SZs
  output(sum_n_SYa) <- sum_n_SYa
  output(sum_n_SZa) <- sum_n_SZa
  output(sum_n_SYb) <- sum_n_SYb
  output(sum_n_SZb) <- sum_n_SZb
  output(sum_n_SYc) <- sum_n_SYc
  output(sum_n_SZc) <- sum_n_SZc
  output(sum_n_SYab) <- sum_n_SYab
  output(sum_n_SZab) <- sum_n_SZab
  output(sum_n_SYac) <- sum_n_SYac
  output(sum_n_SZac) <- sum_n_SZac
  output(sum_n_SYbc) <- sum_n_SYbc
  output(sum_n_SZbc) <- sum_n_SZbc
  output(sum_n_SYabc) <- sum_n_SYabc
  output(sum_n_SZabc) <- sum_n_SZabc
  output(p_YsOut) <- p_YsOut
  output(p_ZsOut) <- p_ZsOut
  output(p_Ys_tr) <- p_Ys_tr
  output(p_Zs_tr) <- p_Zs_tr
  output(p_YZs_tr_fail) <- p_YZs_tr_fail
  output(p_YZs_YZa) <- p_YZs_YZa
  output(p_YZs_YZb) <- p_YZs_YZb
  output(p_YZs_YZc) <- p_YZs_YZc
  output(p_YaOut) <- p_YaOut
  output(p_YbOut) <- p_YbOut
  output(p_YcOut) <- p_YcOut
  output(p_ZaOut) <- p_ZaOut
  output(p_ZbOut) <- p_ZbOut
  output(p_ZcOut) <- p_ZcOut
  output(p_Ya_nr) <- p_Ya_nr
  output(p_Yb_nr) <- p_Yb_nr
  output(p_Yc_nr) <- p_Yc_nr
  output(p_Za_nr) <- p_Za_nr
  output(p_Zb_nr) <- p_Zb_nr
  output(p_Zc_nr) <- p_Zc_nr
  output(p_Ya_tr_fail) <- p_Ya_tr_fail
  output(p_Yb_tr_fail) <- p_Yb_tr_fail
  output(p_Yc_tr_fail) <- p_Yc_tr_fail
  output(p_Za_tr_fail) <- p_Za_tr_fail
  output(p_Zb_tr_fail) <- p_Zb_tr_fail
  output(p_Zc_tr_fail) <- p_Zc_tr_fail
  output(p_YZaYZab) <- p_YZaYZab
  output(p_YZbYZab) <- p_YZbYZab
  output(p_YZcYZac) <- p_YZcYZac
  output(p_YabOut) <- p_YabOut 
  output(p_YacOut) <- p_YacOut 
  output(p_YbcOut) <- p_YbcOut 
  output(p_ZabOut) <- p_ZabOut 
  output(p_ZacOut) <- p_ZacOut 
  output(p_ZbcOut) <- p_ZbcOut
  output(p_Yab_tr_fail) <- p_Yab_tr_fail
  output(p_Yac_tr_fail) <- p_Yac_tr_fail
  output(p_Ybc_tr_fail) <- p_Ybc_tr_fail
  output(p_Zab_tr_fail) <- p_Zab_tr_fail
  output(p_Zac_tr_fail) <- p_Zac_tr_fail
  output(p_Zbc_tr_fail) <- p_Zbc_tr_fail
  output(sum_n_YsOut) <- sum_n_YsOut
  output(sum_n_ZsOut) <- sum_n_ZsOut
  output(sum_n_Ys_tr) <- sum_n_Ys_tr
  output(sum_n_Ys_nr) <- sum_n_Ys_nr
  output(sum_n_Zs_tr) <- sum_n_Zs_tr
  output(sum_n_Zs_nr) <- sum_n_Zs_nr
  output(sum_n_Ys_tr_fail) <- sum_n_Ys_tr_fail
  output(sum_n_Ys_tr_success) <- sum_n_Ys_tr_success
  output(sum_n_Zs_tr_fail) <- sum_n_Zs_tr_fail
  output(sum_n_Zs_tr_success) <- sum_n_Zs_tr_success
  output(sum_n_YsS) <- sum_n_YsS
  output(sum_n_ZsS) <- sum_n_ZsS
  output(sum_n_YsYa) <- sum_n_YsYa
  output(sum_n_YsYb) <- sum_n_YsYb
  output(sum_n_YsYc) <- sum_n_YsYc
  output(sum_n_ZsZa) <- sum_n_ZsZa
  output(sum_n_ZsZb) <- sum_n_ZsZb
  output(sum_n_ZsZc) <- sum_n_ZsZc
  output(sum_n_YaOut) <- sum_n_YaOut
  output(sum_n_YbOut) <- sum_n_YbOut
  output(sum_n_YcOut) <- sum_n_YcOut
  output(sum_n_ZaOut) <- sum_n_ZaOut
  output(sum_n_ZbOut) <- sum_n_ZbOut
  output(sum_n_ZcOut) <- sum_n_ZcOut
  output(sum_n_Ya_nr) <- sum_n_Ya_nr
  output(sum_n_Ya_tr) <- sum_n_Ya_tr
  output(sum_n_Yb_nr) <- sum_n_Yb_nr
  output(sum_n_Yb_tr) <- sum_n_Yb_tr
  output(sum_n_Yc_nr) <- sum_n_Yc_nr
  output(sum_n_Yc_tr) <- sum_n_Yc_tr
  output(sum_n_Za_nr) <- sum_n_Za_nr
  output(sum_n_Za_tr) <- sum_n_Za_tr
  output(sum_n_Zb_nr) <- sum_n_Zb_nr
  output(sum_n_Zb_tr) <- sum_n_Zb_tr
  output(sum_n_Zc_nr) <- sum_n_Zc_nr
  output(sum_n_Zc_tr) <- sum_n_Zc_tr
  output(sum_n_Ya_tr_fail) <- sum_n_Ya_tr_fail
  output(sum_n_Yb_tr_fail) <- sum_n_Yb_tr_fail
  output(sum_n_Yc_tr_fail) <- sum_n_Yc_tr_fail
  output(sum_n_Za_tr_fail) <- sum_n_Za_tr_fail
  output(sum_n_Zb_tr_fail) <- sum_n_Zb_tr_fail
  output(sum_n_Zc_tr_fail) <- sum_n_Zc_tr_fail
  output(sum_n_Ya_tr_success) <- sum_n_Ya_tr_success
  output(sum_n_Yb_tr_success) <- sum_n_Yb_tr_success
  output(sum_n_Yc_tr_success) <- sum_n_Yc_tr_success
  output(sum_n_Za_tr_success) <- sum_n_Za_tr_success
  output(sum_n_Zb_tr_success) <- sum_n_Zb_tr_success
  output(sum_n_Zc_tr_success) <- sum_n_Zc_tr_success
  output(sum_n_YaS) <- sum_n_YaS
  output(sum_n_YbS) <- sum_n_YbS
  output(sum_n_YcS) <- sum_n_YcS
  output(sum_n_ZaS) <- sum_n_ZaS
  output(sum_n_ZbS) <- sum_n_ZbS
  output(sum_n_ZcS) <- sum_n_ZcS
  output(sum_n_YaYab) <- sum_n_YaYab
  output(sum_n_YaYac) <- sum_n_YaYac
  output(sum_n_YbYab) <- sum_n_YbYab
  output(sum_n_YbYbc) <- sum_n_YbYbc
  output(sum_n_YcYac) <- sum_n_YcYac 
  output(sum_n_YcYbc) <- sum_n_YcYbc 
  output(sum_n_ZaZab) <- sum_n_ZaZab
  output(sum_n_ZaZac) <- sum_n_ZaZac
  output(sum_n_ZbZab) <- sum_n_ZbZab
  output(sum_n_ZbZbc) <- sum_n_ZbZbc
  output(sum_n_ZcZac) <- sum_n_ZcZac 
  output(sum_n_ZcZbc) <- sum_n_ZcZbc 
  output(sum_n_YabOut) <- sum_n_YabOut 
  output(sum_n_YacOut) <- sum_n_YacOut 
  output(sum_n_YbcOut) <- sum_n_YbcOut
  output(sum_n_ZabOut) <- sum_n_ZabOut 
  output(sum_n_ZacOut) <- sum_n_ZacOut 
  output(sum_n_ZbcOut) <- sum_n_ZbcOut
  output(sum_n_YabYabc) <- sum_n_YabYabc
  output(sum_n_YabS) <- sum_n_YabS
  output(sum_n_YacYabc) <- sum_n_YacYabc
  output(sum_n_YacS) <- sum_n_YacS
  output(sum_n_YbcYabc) <- sum_n_YbcYabc
  output(sum_n_YbcS) <- sum_n_YbcS
  output(sum_n_ZabZabc) <- sum_n_ZabZabc
  output(sum_n_ZabS) <- sum_n_ZabS
  output(sum_n_ZacZabc) <- sum_n_ZacZabc
  output(sum_n_ZacS) <- sum_n_ZacS
  output(sum_n_ZbcZabc) <- sum_n_ZbcZabc
  output(sum_n_ZbcS) <- sum_n_ZbcS
  output(sum_n_ZabcOut) <- sum_n_ZabcOut 
  output(p_YabcOut) <- p_YabcOut
  output(p_ZabcOut) <- p_ZabcOut
  output(sum_n_YabcOut) <- sum_n_YabcOut
})

equal_allocation_stochastic_6_5 <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] - n_SOut[i] + n_YsS[i] + n_YaS[i] + n_YbS[i] + n_YcS[i] + n_YabS[i] + n_YacS[i] + n_YbcS[i] +
    n_YabcOut[i] + n_ZsS[i] + n_ZaS[i] + n_ZbS[i] + n_ZcS[i] + n_ZabS[i] + n_ZacS[i] + n_ZbcS[i] + n_ZabcOut[i] 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + n_SYs[i] - n_YsOut[i]
  update(Z_s[]) <- Z_s[i] + n_SZs[i] - n_ZsOut[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + n_SYa[i] + n_YsYa[i] - n_YaOut[i]
  update(Z_a[]) <- Z_a[i] + n_SZa[i] + n_ZsZa[i] - n_ZaOut[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + n_SYb[i] + n_YsYb[i] - n_YbOut[i]
  update(Z_b[]) <- Z_b[i] + n_SZb[i] + n_ZsZb[i] - n_ZbOut[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + n_SYc[i] + n_YsYc[i] - n_YcOut[i] 
  update(Z_c[]) <- Z_c[i] + n_SZc[i] + n_ZsZc[i] - n_ZcOut[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + n_SYab[i] + n_YaYab[i] + n_YbYab[i] - n_YabOut[i]
  update(Z_ab[]) <- Z_ab[i] + n_SZab[i] + n_ZaZab[i] + n_ZbZab[i] - n_ZabOut[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + n_SYac[i] + n_YaYac[i] + n_YcYac[i] - n_YacOut[i]
  update(Z_ac[]) <- Z_ac[i] + n_SZac[i] + n_ZaZac[i] + n_ZcZac[i] - n_ZacOut[i] 
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + n_SYbc[i] + n_YbYbc[i] + n_YcYbc[i] - n_YbcOut[i]
  update(Z_bc[]) <- Z_bc[i] + n_SZbc[i] + n_ZbZbc[i] + n_ZcZbc[i] - n_ZbcOut[i]
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + n_SYabc[i] + n_YabYabc[i] + n_YacYabc[i] + n_YbcYabc[i] - n_YabcOut[i]
  update(Z_abc[]) <- Z_abc[i] + n_SZabc[i] + n_ZabZabc[i] + n_ZacZabc[i] + n_ZbcZabc[i] - n_ZabcOut[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  ##STOCHASTIC calculations
  #calculate the probability of leaving S from the rate
  p_SOut[] <- 1 - exp(-sum_all_BI_terms[i])
  
  #calculate the number of people leaving S
  n_SOut[] <- rbinom(S[i], p_SOut[i])
  
  #divide probabilities leaving S 
  p_SYZa[] <- sum_BI_YZa[i]/sum_all_BI_terms[i]
  p_SYZb[] <- sum_BI_YZb[i]/sum_all_BI_terms[i]
  p_SYZc[] <- sum_BI_YZc[i]/sum_all_BI_terms[i]
  p_SYZab[] <- sum_BI_YZab[i]/sum_all_BI_terms[i]
  p_SYZac[] <- sum_BI_YZac[i]/sum_all_BI_terms[i]
  p_SYZbc[] <- sum_BI_YZbc[i]/sum_all_BI_terms[i]
  p_SYZabc[] <- sum_BI_YZabc[i]/sum_all_BI_terms[i] #fixed error here 6/3
  
  #divide the people leaving S into different flows
  n_SYZa[] <- rbinom(n_SOut[i], p_SYZa[i])
  n_SYZb[] <- rbinom(n_SOut[i] - n_SYZa[i], (p_SYZb[i])/(1 - p_SYZa[i] + 1e-8)) #i) subtract those already subtracted from n_SOut ii) renormalize the probability of exit
  n_SYZc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i], (p_SYZc[i])/(1 - p_SYZa[i] - p_SYZb[i] + 1e-8))
  n_SYZab[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i], (p_SYZab[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] + 1e-8))
  n_SYZac[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i], (p_SYZac[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] + 1e-8)) #fixed n to p issue here and remaining 2 eqs 6/5
  n_SYZbc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i], (p_SYZbc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] + 1e-8)) 
  n_SYZabc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i], (p_SYZabc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] - p_SYZbc[i] + 1e-8))
  n_SYZs[] <- n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] - n_SYZabc[i]
  
  #split into symptomatic vs asymptomatics
  n_SYs[] <- rbinom(n_SYZs[i], sigma)
  n_SZs[] <- n_SYZs[i] - n_SYs[i] #fixed issue here
  n_SYa[] <- rbinom(n_SYZa[i], sigma)
  n_SZa[] <- n_SYZa[i] - n_SYa[i]
  n_SYb[] <- rbinom(n_SYZb[i], sigma)
  n_SZb[] <- n_SYZb[i] - n_SYb[i]
  n_SYc[] <- rbinom(n_SYZc[i], sigma)
  n_SZc[] <- n_SYZc[i] - n_SYc[i]
  n_SYab[] <- rbinom(n_SYZab[i], sigma)
  n_SZab[] <- n_SYZab[i] - n_SYab[i]
  n_SYac[] <- rbinom(n_SYZac[i], sigma)
  n_SZac[] <- n_SYZac[i] - n_SYac[i]
  n_SYbc[] <- rbinom(n_SYZbc[i], sigma)
  n_SZbc[] <- n_SYZbc[i] - n_SYbc[i]
  n_SYabc[] <- rbinom(n_SYZabc[i], sigma)
  n_SZabc[] <- n_SYZabc[i] - n_SYabc[i]
  
  #calculate the probability of leaving Ys and Zs by converting exit rates
  p_YsOut <- 1 - exp(-(T_s + d))
  p_ZsOut <- 1 - exp(-(T_m + d))
  
  #calculate the number of people leaving Ys and Zs based on the probabilities
  n_YsOut[] <- rbinom(Y_s[i], p_YsOut)
  n_ZsOut[] <- rbinom(Z_s[i], p_ZsOut)
  
  #split the people recovering naturally (nr) vs the people leaving due to treatment (tr)
  p_Ys_tr <- T_s / (T_s + d)
  n_Ys_tr[] <- rbinom(n_YsOut[i], p_Ys_tr)
  n_Ys_nr[] <- n_YsOut[i] - n_Ys_tr[i]
  p_Zs_tr <- T_m / (T_m + d)
  n_Zs_tr[] <- rbinom(n_ZsOut[i], p_Zs_tr)
  n_Zs_nr[] <- n_ZsOut[i] - n_Zs_tr[i]
  
  #of the people being treated, split the flows into recovery vs. resistance
  p_YZs_tr_fail <- E_a*w_a + E_b*w_b + E_c*w_c
  n_Ys_tr_fail[] <- rbinom(n_Ys_tr[i], p_YZs_tr_fail)
  n_Ys_tr_success[] <- n_Ys_tr[i] - n_Ys_tr_fail[i]
  n_Zs_tr_fail[] <- rbinom(n_Zs_tr[i], p_YZs_tr_fail)
  n_Zs_tr_success[] <- n_Zs_tr[i] - n_Zs_tr_fail[i]
  
  #combine those recovering from successful treatment and natural recovery as leaving Ys and Zs and going back to S
  n_YsS[] <- n_Ys_tr_success[i] + n_Ys_nr[i]
  n_ZsS[] <- n_Zs_tr_success[i] + n_Zs_nr[i]
  
  #split those who are treated and develop resistance into the three possible flows
  #first take the individually normalized probabilities of each flow
  p_YZs_YZa <- E_a*w_a/p_YZs_tr_fail
  p_YZs_YZb <- E_b*w_b/p_YZs_tr_fail
  p_YZs_YZc <- E_c*w_c/p_YZs_tr_fail
  
  #do stick-breaking with Y
  n_YsYa[] <- rbinom(n_Ys_tr_fail[i], p_YZs_YZa)
  n_YsYb[] <- rbinom(n_Ys_tr_fail[i] - n_YsYa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_YsYc[] <- n_Ys_tr_fail[i] - n_YsYa[i] - n_YsYb[i]
  
  #do stick-breaking with Z
  n_ZsZa[] <- rbinom(n_Zs_tr_fail[i], p_YZs_YZa)
  n_ZsZb[] <- rbinom(n_Zs_tr_fail[i] - n_ZsZa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_ZsZc[] <- n_Zs_tr_fail[i] - n_ZsZa[i] - n_ZsZb[i]
  
  ##now onto single resistance compartments. Probabilities of exit will be different for Y and Z compartments due to retreatment
  #probabilities for exit flows from Ya, Yb, Yc, Za, Zb, and Zc
  p_YaOut <- 1 - exp(-((E_b + E_c)*T_s + E_a*kappa*T_sr + d))
  p_YbOut <- 1 - exp(-((E_a + E_c)*T_s + E_b*kappa*T_sr + d))
  p_YcOut <- 1 - exp(-((E_a + E_b)*T_s + E_c*kappa*T_sr + d))
  p_ZaOut <- 1 - exp(-((E_b + E_c)*T_m + d))
  p_ZbOut <- 1 - exp(-((E_a + E_c)*T_m + d))
  p_ZcOut <- 1 - exp(-((E_a + E_b)*T_m + d))
  
  #calculate those leaving the single resistance compartments
  n_YaOut[] <- rbinom(Y_a[i], p_YaOut)
  n_YbOut[] <- rbinom(Y_b[i], p_YbOut)
  n_YcOut[] <- rbinom(Y_c[i], p_YcOut)
  n_ZaOut[] <- rbinom(Z_a[i], p_ZaOut)
  n_ZbOut[] <- rbinom(Z_b[i], p_ZbOut)
  n_ZcOut[] <- rbinom(Z_c[i], p_ZcOut)
  
  #calculate the probabilities of leaving single resistance compartments due to natural recovery
  p_Ya_nr <- d/((E_b + E_c)*T_s + E_a*kappa*T_sr + d)
  p_Yb_nr <- d/((E_a + E_c)*T_s + E_b*kappa*T_sr + d)
  p_Yc_nr <- d/((E_a + E_b)*T_s + E_c*kappa*T_sr + d)
  p_Za_nr <- d/((E_b + E_c)*T_m + d)
  p_Zb_nr <- d/((E_a + E_c)*T_m + d)
  p_Zc_nr <- d/((E_a + E_b)*T_m + d)
  
  #split the people leaving single resistance into natural recovery vs treatment
  n_Ya_nr[] <- rbinom(n_YaOut[i], p_Ya_nr)
  n_Ya_tr[] <- n_YaOut[i] - n_Ya_nr[i]
  n_Yb_nr[] <- rbinom(n_YbOut[i], p_Yb_nr)
  n_Yb_tr[] <- n_YbOut[i] - n_Yb_nr[i]
  n_Yc_nr[] <- rbinom(n_YcOut[i], p_Yc_nr)
  n_Yc_tr[] <- n_YcOut[i] - n_Yc_nr[i]
  n_Za_nr[] <- rbinom(n_ZaOut[i], p_Za_nr)
  n_Za_tr[] <- n_ZaOut[i] - n_Za_nr[i]
  n_Zb_nr[] <- rbinom(n_ZbOut[i], p_Zb_nr)
  n_Zb_tr[] <- n_ZbOut[i] - n_Zb_nr[i]
  n_Zc_nr[] <- rbinom(n_ZcOut[i], p_Zc_nr)
  n_Zc_tr[] <- n_ZcOut[i] - n_Zc_nr[i]
  
  #split those leaving single resistance after treatment into those who fail treatment vs. those who recover
  #note that here, treatment failure means people who develop resistance on treatment, not those treated with the drug they are resistance bc they stay in the resistant compartment
  p_Ya_tr_fail <- (E_b*w_b*T_s + E_c*w_c*T_s)/((E_b + E_c)*T_s + E_a*kappa*T_sr)
  p_Yb_tr_fail <- (E_a*w_a*T_s + E_c*w_c*T_s)/((E_a + E_c)*T_s + E_b*kappa*T_sr)
  p_Yc_tr_fail <- (E_a*w_a*T_s + E_b*w_b*T_s)/((E_a + E_b)*T_s + E_c*kappa*T_sr)
  p_Za_tr_fail <- (E_b*w_b + E_c*w_c)/ ((E_b + E_c) + 1e-8) #add very small term to denominator so this probability is always defined
  p_Zb_tr_fail <- (E_a*w_a + E_c*w_c)/ ((E_a + E_c) + 1e-8)
  p_Zc_tr_fail <- (E_a*w_a + E_b*w_b)/ ((E_a + E_b) + 1e-8)
  
  n_Ya_tr_fail[] <- rbinom(n_Ya_tr[i], p_Ya_tr_fail)
  n_Yb_tr_fail[] <- rbinom(n_Yb_tr[i], p_Yb_tr_fail)
  n_Yc_tr_fail[] <- rbinom(n_Yc_tr[i], p_Yc_tr_fail)
  n_Za_tr_fail[] <- rbinom(n_Za_tr[i], p_Za_tr_fail)
  n_Zb_tr_fail[] <- rbinom(n_Zb_tr[i], p_Zb_tr_fail)
  n_Zc_tr_fail[] <- rbinom(n_Zc_tr[i], p_Zc_tr_fail)
  
  n_Ya_tr_success[] <- n_Ya_tr[i] - n_Ya_tr_fail[i]
  n_Yb_tr_success[] <- n_Yb_tr[i] - n_Yb_tr_fail[i]
  n_Yc_tr_success[] <- n_Yc_tr[i] - n_Yc_tr_fail[i]
  n_Za_tr_success[] <- n_Za_tr[i] - n_Za_tr_fail[i]
  n_Zb_tr_success[] <- n_Zb_tr[i] - n_Zb_tr_fail[i]
  n_Zc_tr_success[] <- n_Zc_tr[i] - n_Zc_tr_fail[i]
  
  #combine treatment success and natural recover terms to get single resistance recovery flows
  n_YaS[] <- n_Ya_tr_success[i] + n_Ya_nr[i]
  n_YbS[] <- n_Yb_tr_success[i] + n_Yb_nr[i]
  n_YcS[] <- n_Yc_tr_success[i] + n_Yc_nr[i]
  n_ZaS[] <- n_Za_tr_success[i] + n_Za_nr[i]
  n_ZbS[] <- n_Zb_tr_success[i] + n_Zb_nr[i]
  n_ZcS[] <- n_Zc_tr_success[i] + n_Zc_nr[i]
  
  #now split the treatment failures into the two different resistance compartment flows
  p_YZaYZab <- E_b*w_b / (E_b*w_b + E_c*w_c + 1e-8) #add very small term to denominator so this probability is always defined
  p_YZbYZab <- E_a*w_a / (E_a*w_a + E_c*w_c + 1e-8)
  p_YZcYZac <- E_a*w_a / (E_a*w_a + E_b*w_b + 1e-8)
  
  n_YaYab[] <- rbinom(n_Ya_tr_fail[i], p_YZaYZab)
  n_YaYac[] <- n_Ya_tr_fail[i] - n_YaYab[i]
  n_YbYab[] <- rbinom(n_Yb_tr_fail[i], p_YZbYZab)
  n_YbYbc[] <- n_Yb_tr_fail[i] - n_YbYab[i]
  n_YcYac[] <- rbinom(n_Yc_tr_fail[i], p_YZcYZac)
  n_YcYbc[] <- n_Yc_tr_fail[i] - n_YcYac[i]
  
  n_ZaZab[] <- rbinom(n_Za_tr_fail[i], p_YZaYZab)
  n_ZaZac[] <- n_Za_tr_fail[i] - n_ZaZab[i]
  n_ZbZab[] <- rbinom(n_Zb_tr_fail[i], p_YZbYZab)
  n_ZbZbc[] <- n_Zb_tr_fail[i] - n_ZbZab[i]
  n_ZcZac[] <- rbinom(n_Zc_tr_fail[i], p_YZcYZac)
  n_ZcZbc[] <- n_Zc_tr_fail[i] - n_ZcZac[i]
  
  ##now move to the double resistance compartments
  #first calculate total outflow -- calculate exit probabilities and then draw
  p_YabOut <- 1 - exp(-(E_c*T_s + (1 - E_c)*kappa*T_sr + d))
  p_YacOut <- 1 - exp(-(E_b*T_s + (1 - E_b)*kappa*T_sr + d))
  p_YbcOut <- 1 - exp(-(E_a*T_s + (1 - E_a)*kappa*T_sr + d))
  p_ZabOut <- 1 - exp(-(E_c*T_m + d))
  p_ZacOut <- 1 - exp(-(E_b*T_m + d))
  p_ZbcOut <- 1 - exp(-(E_a*T_m + d))
  
  n_YabOut[] <- rbinom(Y_ab[i], p_YabOut)
  n_YacOut[] <- rbinom(Y_ac[i], p_YacOut)
  n_YbcOut[] <- rbinom(Y_bc[i], p_YbcOut)
  n_ZabOut[] <- rbinom(Z_ab[i], p_ZabOut)
  n_ZacOut[] <- rbinom(Z_ac[i], p_ZacOut)
  n_ZbcOut[] <- rbinom(Z_bc[i], p_ZbcOut)
  
  #because these are already dual resistance, the only ones leaving and not recovering are developing resistance to the same drug, so separate those and the rest recover
  p_Yab_tr_fail <- (E_c*w_c*T_s)/(E_c*T_s + (1 - E_c)*kappa*T_sr + d)
  p_Yac_tr_fail <- (E_b*w_b*T_s)/(E_b*T_s + (1 - E_b)*kappa*T_sr + d)
  p_Ybc_tr_fail <- (E_a*w_a*T_s)/(E_a*T_s + (1 - E_a)*kappa*T_sr + d)
  p_Zab_tr_fail <- (E_c*w_c*T_m)/(E_c*T_m + d)
  p_Zac_tr_fail <- (E_b*w_b*T_m)/(E_b*T_m + d)
  p_Zbc_tr_fail <- (E_a*w_a*T_m)/(E_a*T_m + d)
  
  #split out resistant flows to Yabc vs recoveries to S
  n_YabYabc[] <- rbinom(n_YabOut[i], p_Yab_tr_fail)
  n_YabS[] <- n_YabOut[i] - n_YabYabc[i]
  n_YacYabc[] <- rbinom(n_YacOut[i], p_Yac_tr_fail)
  n_YacS[] <- n_YacOut[i] - n_YacYabc[i]
  n_YbcYabc[] <- rbinom(n_YbcOut[i], p_Ybc_tr_fail)
  n_YbcS[] <- n_YbcOut[i] - n_YbcYabc[i]
  
  n_ZabZabc[] <- rbinom(n_ZabOut[i], p_Zab_tr_fail)
  n_ZabS[] <- n_ZabOut[i] - n_ZabZabc[i]
  n_ZacZabc[] <- rbinom(n_ZacOut[i], p_Zac_tr_fail)
  n_ZacS[] <- n_ZacOut[i] - n_ZacZabc[i]
  n_ZbcZabc[] <- rbinom(n_ZbcOut[i], p_Zbc_tr_fail)
  n_ZbcS[] <- n_ZbcOut[i] - n_ZbcZabc[i]
  
  #recovery from triple resistance states
  p_YabcOut <- 1 - exp(-(kappa*T_sr + d))
  p_ZabcOut <- 1 - exp(-d)
  n_YabcOut[] <- rbinom(Y_abc[i], p_YabcOut)
  n_ZabcOut[] <- rbinom(Z_abc[i], p_ZabcOut)
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #make indicator variables for which ordering we're in 
  all_off <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_A <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_B <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  only_C <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1)
  A_B_on <- (lock_threshold_A >= 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  update(C_remain_lock) <- C_remain_lock + A_B_on
  A_C_on <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1) 
  update(B_remain_lock) <- B_remain_lock + A_C_on
  B_C_on <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C >= 1) 
  update(A_remain_lock) <- A_remain_lock + B_C_on
  
  E_a <- 1/3*all_off + 0.5*(only_B + only_C) + max((B_C_on), (A_remain_lock >= 1))
  E_b <- 1/3*all_off + 0.5*(only_A + only_C) + max((A_C_on),(B_remain_lock >= 1))
  E_c <- 1/3*all_off + 0.5*(only_A + only_B) + max((A_B_on),(C_remain_lock >= 1))
  
  sum_tx_probs = E_a + E_c + E_b 
  
  #sum array variables for book-keeping / accounting outputs
  sum_n_SOut <- sum(n_SOut)
  sum_n_SYZa <- sum(n_SYZa)
  sum_n_SYZb <- sum(n_SYZb)
  sum_n_SYZc <- sum(n_SYZc)
  sum_n_SYZab <- sum(n_SYZab)
  sum_n_SYZac <- sum(n_SYZac)
  sum_n_SYZbc <- sum(n_SYZbc)
  sum_n_SYZabc <- sum(n_SYZabc)
  sum_n_SYZs <- sum(n_SYZs)
  sum_n_SYs <- sum(n_SYs)
  sum_n_SZs <- sum(n_SZs)
  sum_n_SYa <- sum(n_SYa)
  sum_n_SZa <- sum(n_SZa)
  sum_n_SYb <- sum(n_SYb)
  sum_n_SZb <- sum(n_SZb)
  sum_n_SYc <- sum(n_SYc)
  sum_n_SZc <- sum(n_SZc)
  sum_n_SYab <- sum(n_SYab)
  sum_n_SZab <- sum(n_SZab)
  sum_n_SYac <- sum(n_SYac)
  sum_n_SZac <- sum(n_SZac)
  sum_n_SYbc <- sum(n_SYbc)
  sum_n_SZbc <- sum(n_SZbc)
  sum_n_SYabc <- sum(n_SYabc)
  sum_n_SZabc <- sum(n_SZabc)
  sum_n_YsOut <- sum(n_YsOut)
  sum_n_ZsOut <- sum(n_ZsOut)
  sum_n_Ys_tr <- sum(n_Ys_tr)
  sum_n_Ys_nr <- sum(n_Ys_nr)
  sum_n_Zs_tr <- sum(n_Zs_tr)
  sum_n_Zs_nr <- sum(n_Zs_nr)
  sum_n_Ys_tr_fail <- sum(n_Ys_tr_fail)
  sum_n_Ys_tr_success <- sum(n_Ys_tr_success)
  sum_n_Zs_tr_fail <- sum(n_Zs_tr_fail)
  sum_n_Zs_tr_success <- sum(n_Zs_tr_success)
  sum_n_YsS <- sum(n_YsS)
  sum_n_ZsS <- sum(n_ZsS)
  sum_n_YsYa <- sum(n_YsYa)
  sum_n_YsYb <- sum(n_YsYb)
  sum_n_YsYc <- sum(n_YsYc)
  sum_n_ZsZa <- sum(n_ZsZa)
  sum_n_ZsZb <- sum(n_ZsZb)
  sum_n_ZsZc <- sum(n_ZsZc)
  sum_n_YaOut <- sum(n_YaOut)
  sum_n_YbOut <- sum(n_YbOut)
  sum_n_YcOut <- sum(n_YcOut)
  sum_n_ZaOut <- sum(n_ZaOut)
  sum_n_ZbOut <- sum(n_ZbOut)
  sum_n_ZcOut <- sum(n_ZcOut)
  sum_n_Ya_nr <- sum(n_Ya_nr)
  sum_n_Ya_tr <- sum(n_Ya_tr)
  sum_n_Yb_nr <- sum(n_Yb_nr)
  sum_n_Yb_tr <- sum(n_Yb_tr)
  sum_n_Yc_nr <- sum(n_Yc_nr)
  sum_n_Yc_tr <- sum(n_Yc_tr)
  sum_n_Za_nr <- sum(n_Za_nr)
  sum_n_Za_tr <- sum(n_Za_tr)
  sum_n_Zb_nr <- sum(n_Zb_nr)
  sum_n_Zb_tr <- sum(n_Zb_tr)
  sum_n_Zc_nr <- sum(n_Zc_nr)
  sum_n_Zc_tr <- sum(n_Zc_tr)
  sum_n_Ya_tr_fail <- sum(n_Ya_tr_fail)
  sum_n_Yb_tr_fail <- sum(n_Yb_tr_fail)
  sum_n_Yc_tr_fail <- sum(n_Yc_tr_fail)
  sum_n_Za_tr_fail <- sum(n_Za_tr_fail)
  sum_n_Zb_tr_fail <- sum(n_Zb_tr_fail)
  sum_n_Zc_tr_fail <- sum(n_Zc_tr_fail)
  sum_n_Ya_tr_success <- sum(n_Ya_tr_success)
  sum_n_Yb_tr_success <- sum(n_Yb_tr_success)
  sum_n_Yc_tr_success <- sum(n_Yc_tr_success)
  sum_n_Za_tr_success <- sum(n_Za_tr_success)
  sum_n_Zb_tr_success <- sum(n_Zb_tr_success)
  sum_n_Zc_tr_success <- sum(n_Zc_tr_success)
  sum_n_YaS <- sum(n_YaS)
  sum_n_YbS <- sum(n_YbS)
  sum_n_YcS <- sum(n_YcS)
  sum_n_ZaS <- sum(n_ZaS)
  sum_n_ZbS <- sum(n_ZbS)
  sum_n_ZcS <- sum(n_ZcS)
  sum_n_YaYab <- sum(n_YaYab)
  sum_n_YaYac <- sum(n_YaYac)
  sum_n_YbYab <- sum(n_YbYab)
  sum_n_YbYbc <- sum(n_YbYbc)
  sum_n_YcYac <- sum(n_YcYac)
  sum_n_YcYbc <- sum(n_YcYbc)
  sum_n_ZaZab <- sum(n_ZaZab)
  sum_n_ZaZac <- sum(n_ZaZac)
  sum_n_ZbZab <- sum(n_ZbZab)
  sum_n_ZbZbc <- sum(n_ZbZbc)
  sum_n_ZcZac <- sum(n_ZcZac)
  sum_n_ZcZbc <- sum(n_ZcZbc)
  sum_n_YabOut <- sum(n_YabOut)
  sum_n_YacOut <- sum(n_YacOut)
  sum_n_YbcOut <- sum(n_YbcOut)
  sum_n_ZabOut <- sum(n_ZabOut)
  sum_n_ZacOut <- sum(n_ZacOut)
  sum_n_ZbcOut <- sum(n_ZbcOut)
  sum_n_YabYabc <- sum(n_YabYabc)
  sum_n_YabS <- sum(n_YabS)
  sum_n_YacYabc <- sum(n_YacYabc)
  sum_n_YacS <- sum(n_YacS)
  sum_n_YbcYabc <- sum(n_YbcYabc)
  sum_n_YbcS<- sum(n_YbcS)
  sum_n_ZabZabc <- sum(n_ZabZabc)
  sum_n_ZabS <- sum(n_ZabS)
  sum_n_ZacZabc <- sum(n_ZacZabc)
  sum_n_ZacS <- sum(n_ZacS)
  sum_n_ZbcZabc <- sum(n_ZbcZabc)
  sum_n_ZbcS<- sum(n_ZbcS)
  sum_n_YabcOut <- sum(n_YabcOut)
  sum_n_ZabcOut <- sum(n_ZabcOut)
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  initial(A_remain_lock) <- 0
  initial(B_remain_lock) <- 0
  initial(C_remain_lock) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  dim(n_SOut) <- N_risk
  dim(p_SOut) <- N_risk
  dim(p_SYZa) <- N_risk
  dim(p_SYZb) <- N_risk
  dim(p_SYZc) <- N_risk
  dim(p_SYZab) <- N_risk
  dim(p_SYZac) <- N_risk  
  dim(p_SYZbc) <- N_risk
  dim(p_SYZabc) <- N_risk
  dim(n_SYZa) <- N_risk
  dim(n_SYZb) <- N_risk
  dim(n_SYZc) <- N_risk
  dim(n_SYZab) <- N_risk
  dim(n_SYZac) <- N_risk
  dim(n_SYZbc) <- N_risk
  dim(n_SYZabc) <- N_risk
  dim(n_SYZs) <- N_risk
  dim(n_SYs) <- N_risk
  dim(n_SZs) <- N_risk
  dim(n_SYa) <- N_risk
  dim(n_SZa) <- N_risk
  dim(n_SYb) <- N_risk
  dim(n_SZb) <- N_risk
  dim(n_SYc) <- N_risk
  dim(n_SZc) <- N_risk
  dim(n_SYab) <- N_risk
  dim(n_SZab) <- N_risk
  dim(n_SYac) <- N_risk
  dim(n_SZac) <- N_risk
  dim(n_SYbc) <- N_risk
  dim(n_SZbc) <- N_risk
  dim(n_SYabc) <- N_risk
  dim(n_SZabc) <- N_risk
  dim(n_YsOut) <- N_risk
  dim(n_ZsOut) <- N_risk
  dim(n_Ys_tr) <- N_risk
  dim(n_Ys_nr) <- N_risk
  dim(n_Zs_tr) <- N_risk
  dim(n_Zs_nr) <- N_risk
  dim(n_Ys_tr_fail) <- N_risk
  dim(n_Ys_tr_success) <- N_risk
  dim(n_Zs_tr_fail) <- N_risk
  dim(n_Zs_tr_success) <- N_risk
  dim(n_YsS) <- N_risk
  dim(n_ZsS) <- N_risk
  dim(n_YsYa) <- N_risk
  dim(n_YsYb) <- N_risk
  dim(n_YsYc) <- N_risk
  dim(n_ZsZa) <- N_risk
  dim(n_ZsZb) <- N_risk
  dim(n_ZsZc) <- N_risk
  dim(n_YaOut) <- N_risk
  dim(n_YbOut) <- N_risk
  dim(n_YcOut) <- N_risk
  dim(n_ZaOut) <- N_risk
  dim(n_ZbOut) <- N_risk
  dim(n_ZcOut) <- N_risk
  dim(n_Ya_nr) <- N_risk
  dim(n_Ya_tr) <- N_risk
  dim(n_Yb_nr) <- N_risk
  dim(n_Yb_tr) <- N_risk
  dim(n_Yc_nr) <- N_risk
  dim(n_Yc_tr) <- N_risk
  dim(n_Za_nr) <- N_risk
  dim(n_Za_tr) <- N_risk
  dim(n_Zb_nr) <- N_risk
  dim(n_Zb_tr) <- N_risk
  dim(n_Zc_nr) <- N_risk
  dim(n_Zc_tr) <- N_risk
  dim(n_Ya_tr_fail) <- N_risk
  dim(n_Yb_tr_fail) <- N_risk
  dim(n_Yc_tr_fail) <- N_risk
  dim(n_Za_tr_fail) <- N_risk
  dim(n_Zb_tr_fail) <- N_risk
  dim(n_Zc_tr_fail) <- N_risk
  dim(n_Ya_tr_success) <- N_risk
  dim(n_Yb_tr_success) <- N_risk
  dim(n_Yc_tr_success) <- N_risk
  dim(n_Za_tr_success) <- N_risk
  dim(n_Zb_tr_success) <- N_risk
  dim(n_Zc_tr_success) <- N_risk
  dim(n_YaS) <- N_risk
  dim(n_YbS) <- N_risk
  dim(n_YcS) <- N_risk
  dim(n_ZaS) <- N_risk
  dim(n_ZbS) <- N_risk
  dim(n_ZcS) <- N_risk
  dim(n_YaYab) <- N_risk
  dim(n_YaYac) <- N_risk
  dim(n_YbYab) <- N_risk
  dim(n_YbYbc) <- N_risk
  dim(n_YcYac) <- N_risk
  dim(n_YcYbc) <- N_risk
  dim(n_ZaZab) <- N_risk
  dim(n_ZaZac) <- N_risk
  dim(n_ZbZab) <- N_risk
  dim(n_ZbZbc) <- N_risk
  dim(n_ZcZac) <- N_risk
  dim(n_ZcZbc) <- N_risk
  dim(n_YabOut) <- N_risk
  dim(n_YacOut) <- N_risk
  dim(n_YbcOut) <- N_risk
  dim(n_ZabOut) <- N_risk
  dim(n_ZacOut) <- N_risk
  dim(n_ZbcOut) <- N_risk
  dim(n_YabYabc) <- N_risk
  dim(n_YabS) <- N_risk
  dim(n_YacYabc) <- N_risk
  dim(n_YacS) <- N_risk
  dim(n_YbcYabc) <- N_risk
  dim(n_YbcS) <- N_risk
  dim(n_ZabZabc) <- N_risk
  dim(n_ZabS) <- N_risk
  dim(n_ZacZabc) <- N_risk
  dim(n_ZacS) <- N_risk
  dim(n_ZbcZabc) <- N_risk
  dim(n_ZbcS) <- N_risk
  dim(n_YabcOut) <- N_risk
  dim(n_ZabcOut) <- N_risk
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  
  #define desired non-state outputs
  output(ninf_Ys) <- ninf_Ys
  output(ninf_Zs) <- ninf_Zs
  output(ninf_Ya) <- ninf_Ya
  output(ninf_Za) <- ninf_Za
  output(ninf_Yb) <- ninf_Yb
  output(ninf_Zb) <- ninf_Zb
  output(ninf_Yab) <- ninf_Yab
  output(ninf_Zab) <- ninf_Zab
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(total_sum_BI_S) <- total_sum_BI_S
  output(total_sum_BI_YZs) <- total_sum_BI_YZs
  output(total_sum_BI_YZa) <- total_sum_BI_YZa
  output(total_sum_BI_YZb) <- total_sum_BI_YZb
  output(total_sum_BI_YZc) <- total_sum_BI_YZc
  output(total_sum_BI_YZab) <- total_sum_BI_YZab
  output(total_sum_BI_YZac) <- total_sum_BI_YZac
  output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
  output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
  output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev
  output(sum_tx_probs) <- sum_tx_probs
  output(p_SOut) <- p_SOut
  output(sum_n_SOut) <- sum_n_SOut
  output(p_SYZa) <- p_SYZa
  output(p_SYZb) <- p_SYZb
  output(p_SYZc) <- p_SYZc
  output(p_SYZab) <- p_SYZab
  output(p_SYZac) <- p_SYZac
  output(p_SYZbc) <- p_SYZbc
  output(p_SYZabc) <- p_SYZabc
  output(sum_n_SYZa) <- sum_n_SYZa
  output(sum_n_SYZb) <- sum_n_SYZb
  output(sum_n_SYZc) <- sum_n_SYZc
  output(sum_n_SYZab) <- sum_n_SYZab
  output(sum_n_SYZac) <- sum_n_SYZac
  output(sum_n_SYZbc) <- sum_n_SYZbc
  output(sum_n_SYZabc) <- sum_n_SYZabc
  output(sum_n_SYZs) <- sum_n_SYZs
  output(sum_n_SYs) <- sum_n_SYs
  output(sum_n_SZs) <- sum_n_SZs
  output(sum_n_SYa) <- sum_n_SYa
  output(sum_n_SZa) <- sum_n_SZa
  output(sum_n_SYb) <- sum_n_SYb
  output(sum_n_SZb) <- sum_n_SZb
  output(sum_n_SYc) <- sum_n_SYc
  output(sum_n_SZc) <- sum_n_SZc
  output(sum_n_SYab) <- sum_n_SYab
  output(sum_n_SZab) <- sum_n_SZab
  output(sum_n_SYac) <- sum_n_SYac
  output(sum_n_SZac) <- sum_n_SZac
  output(sum_n_SYbc) <- sum_n_SYbc
  output(sum_n_SZbc) <- sum_n_SZbc
  output(sum_n_SYabc) <- sum_n_SYabc
  output(sum_n_SZabc) <- sum_n_SZabc
  output(p_YsOut) <- p_YsOut
  output(p_ZsOut) <- p_ZsOut
  output(p_Ys_tr) <- p_Ys_tr
  output(p_Zs_tr) <- p_Zs_tr
  output(p_YZs_tr_fail) <- p_YZs_tr_fail
  output(p_YZs_YZa) <- p_YZs_YZa
  output(p_YZs_YZb) <- p_YZs_YZb
  output(p_YZs_YZc) <- p_YZs_YZc
  output(p_YaOut) <- p_YaOut
  output(p_YbOut) <- p_YbOut
  output(p_YcOut) <- p_YcOut
  output(p_ZaOut) <- p_ZaOut
  output(p_ZbOut) <- p_ZbOut
  output(p_ZcOut) <- p_ZcOut
  output(p_Ya_nr) <- p_Ya_nr
  output(p_Yb_nr) <- p_Yb_nr
  output(p_Yc_nr) <- p_Yc_nr
  output(p_Za_nr) <- p_Za_nr
  output(p_Zb_nr) <- p_Zb_nr
  output(p_Zc_nr) <- p_Zc_nr
  output(p_Ya_tr_fail) <- p_Ya_tr_fail
  output(p_Yb_tr_fail) <- p_Yb_tr_fail
  output(p_Yc_tr_fail) <- p_Yc_tr_fail
  output(p_Za_tr_fail) <- p_Za_tr_fail
  output(p_Zb_tr_fail) <- p_Zb_tr_fail
  output(p_Zc_tr_fail) <- p_Zc_tr_fail
  output(p_YZaYZab) <- p_YZaYZab
  output(p_YZbYZab) <- p_YZbYZab
  output(p_YZcYZac) <- p_YZcYZac
  output(p_YabOut) <- p_YabOut 
  output(p_YacOut) <- p_YacOut 
  output(p_YbcOut) <- p_YbcOut 
  output(p_ZabOut) <- p_ZabOut 
  output(p_ZacOut) <- p_ZacOut 
  output(p_ZbcOut) <- p_ZbcOut
  output(p_Yab_tr_fail) <- p_Yab_tr_fail
  output(p_Yac_tr_fail) <- p_Yac_tr_fail
  output(p_Ybc_tr_fail) <- p_Ybc_tr_fail
  output(p_Zab_tr_fail) <- p_Zab_tr_fail
  output(p_Zac_tr_fail) <- p_Zac_tr_fail
  output(p_Zbc_tr_fail) <- p_Zbc_tr_fail
  output(sum_n_YsOut) <- sum_n_YsOut
  output(sum_n_ZsOut) <- sum_n_ZsOut
  output(sum_n_Ys_tr) <- sum_n_Ys_tr
  output(sum_n_Ys_nr) <- sum_n_Ys_nr
  output(sum_n_Zs_tr) <- sum_n_Zs_tr
  output(sum_n_Zs_nr) <- sum_n_Zs_nr
  output(sum_n_Ys_tr_fail) <- sum_n_Ys_tr_fail
  output(sum_n_Ys_tr_success) <- sum_n_Ys_tr_success
  output(sum_n_Zs_tr_fail) <- sum_n_Zs_tr_fail
  output(sum_n_Zs_tr_success) <- sum_n_Zs_tr_success
  output(sum_n_YsS) <- sum_n_YsS
  output(sum_n_ZsS) <- sum_n_ZsS
  output(sum_n_YsYa) <- sum_n_YsYa
  output(sum_n_YsYb) <- sum_n_YsYb
  output(sum_n_YsYc) <- sum_n_YsYc
  output(sum_n_ZsZa) <- sum_n_ZsZa
  output(sum_n_ZsZb) <- sum_n_ZsZb
  output(sum_n_ZsZc) <- sum_n_ZsZc
  output(sum_n_YaOut) <- sum_n_YaOut
  output(sum_n_YbOut) <- sum_n_YbOut
  output(sum_n_YcOut) <- sum_n_YcOut
  output(sum_n_ZaOut) <- sum_n_ZaOut
  output(sum_n_ZbOut) <- sum_n_ZbOut
  output(sum_n_ZcOut) <- sum_n_ZcOut
  output(sum_n_Ya_nr) <- sum_n_Ya_nr
  output(sum_n_Ya_tr) <- sum_n_Ya_tr
  output(sum_n_Yb_nr) <- sum_n_Yb_nr
  output(sum_n_Yb_tr) <- sum_n_Yb_tr
  output(sum_n_Yc_nr) <- sum_n_Yc_nr
  output(sum_n_Yc_tr) <- sum_n_Yc_tr
  output(sum_n_Za_nr) <- sum_n_Za_nr
  output(sum_n_Za_tr) <- sum_n_Za_tr
  output(sum_n_Zb_nr) <- sum_n_Zb_nr
  output(sum_n_Zb_tr) <- sum_n_Zb_tr
  output(sum_n_Zc_nr) <- sum_n_Zc_nr
  output(sum_n_Zc_tr) <- sum_n_Zc_tr
  output(sum_n_Ya_tr_fail) <- sum_n_Ya_tr_fail
  output(sum_n_Yb_tr_fail) <- sum_n_Yb_tr_fail
  output(sum_n_Yc_tr_fail) <- sum_n_Yc_tr_fail
  output(sum_n_Za_tr_fail) <- sum_n_Za_tr_fail
  output(sum_n_Zb_tr_fail) <- sum_n_Zb_tr_fail
  output(sum_n_Zc_tr_fail) <- sum_n_Zc_tr_fail
  output(sum_n_Ya_tr_success) <- sum_n_Ya_tr_success
  output(sum_n_Yb_tr_success) <- sum_n_Yb_tr_success
  output(sum_n_Yc_tr_success) <- sum_n_Yc_tr_success
  output(sum_n_Za_tr_success) <- sum_n_Za_tr_success
  output(sum_n_Zb_tr_success) <- sum_n_Zb_tr_success
  output(sum_n_Zc_tr_success) <- sum_n_Zc_tr_success
  output(sum_n_YaS) <- sum_n_YaS
  output(sum_n_YbS) <- sum_n_YbS
  output(sum_n_YcS) <- sum_n_YcS
  output(sum_n_ZaS) <- sum_n_ZaS
  output(sum_n_ZbS) <- sum_n_ZbS
  output(sum_n_ZcS) <- sum_n_ZcS
  output(sum_n_YaYab) <- sum_n_YaYab
  output(sum_n_YaYac) <- sum_n_YaYac
  output(sum_n_YbYab) <- sum_n_YbYab
  output(sum_n_YbYbc) <- sum_n_YbYbc
  output(sum_n_YcYac) <- sum_n_YcYac 
  output(sum_n_YcYbc) <- sum_n_YcYbc 
  output(sum_n_ZaZab) <- sum_n_ZaZab
  output(sum_n_ZaZac) <- sum_n_ZaZac
  output(sum_n_ZbZab) <- sum_n_ZbZab
  output(sum_n_ZbZbc) <- sum_n_ZbZbc
  output(sum_n_ZcZac) <- sum_n_ZcZac 
  output(sum_n_ZcZbc) <- sum_n_ZcZbc 
  output(sum_n_YabOut) <- sum_n_YabOut 
  output(sum_n_YacOut) <- sum_n_YacOut 
  output(sum_n_YbcOut) <- sum_n_YbcOut
  output(sum_n_ZabOut) <- sum_n_ZabOut 
  output(sum_n_ZacOut) <- sum_n_ZacOut 
  output(sum_n_ZbcOut) <- sum_n_ZbcOut
  output(sum_n_YabYabc) <- sum_n_YabYabc
  output(sum_n_YabS) <- sum_n_YabS
  output(sum_n_YacYabc) <- sum_n_YacYabc
  output(sum_n_YacS) <- sum_n_YacS
  output(sum_n_YbcYabc) <- sum_n_YbcYabc
  output(sum_n_YbcS) <- sum_n_YbcS
  output(sum_n_ZabZabc) <- sum_n_ZabZabc
  output(sum_n_ZabS) <- sum_n_ZabS
  output(sum_n_ZacZabc) <- sum_n_ZacZabc
  output(sum_n_ZacS) <- sum_n_ZacS
  output(sum_n_ZbcZabc) <- sum_n_ZbcZabc
  output(sum_n_ZbcS) <- sum_n_ZbcS
  output(sum_n_ZabcOut) <- sum_n_ZabcOut 
  output(p_YabcOut) <- p_YabcOut
  output(p_ZabcOut) <- p_ZabcOut
  output(sum_n_YabcOut) <- sum_n_YabcOut
})


#only necessary outputs
sequential_stochastic_6_7_brief <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] - n_SOut[i] + n_YsS[i] + n_YaS[i] + n_YbS[i] + n_YcS[i] + n_YabS[i] + n_YacS[i] + n_YbcS[i] +
    n_YabcOut[i] + n_ZsS[i] + n_ZaS[i] + n_ZbS[i] + n_ZcS[i] + n_ZabS[i] + n_ZacS[i] + n_ZbcS[i] + n_ZabcOut[i] 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + n_SYs[i] - n_YsOut[i]
  update(Z_s[]) <- Z_s[i] + n_SZs[i] - n_ZsOut[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + n_SYa[i] + n_YsYa[i] - n_YaOut[i]
  update(Z_a[]) <- Z_a[i] + n_SZa[i] + n_ZsZa[i] - n_ZaOut[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + n_SYb[i] + n_YsYb[i] - n_YbOut[i]
  update(Z_b[]) <- Z_b[i] + n_SZb[i] + n_ZsZb[i] - n_ZbOut[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + n_SYc[i] + n_YsYc[i] - n_YcOut[i] 
  update(Z_c[]) <- Z_c[i] + n_SZc[i] + n_ZsZc[i] - n_ZcOut[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + n_SYab[i] + n_YaYab[i] + n_YbYab[i] - n_YabOut[i]
  update(Z_ab[]) <- Z_ab[i] + n_SZab[i] + n_ZaZab[i] + n_ZbZab[i] - n_ZabOut[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + n_SYac[i] + n_YaYac[i] + n_YcYac[i] - n_YacOut[i]
  update(Z_ac[]) <- Z_ac[i] + n_SZac[i] + n_ZaZac[i] + n_ZcZac[i] - n_ZacOut[i] 
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + n_SYbc[i] + n_YbYbc[i] + n_YcYbc[i] - n_YbcOut[i]
  update(Z_bc[]) <- Z_bc[i] + n_SZbc[i] + n_ZbZbc[i] + n_ZcZbc[i] - n_ZbcOut[i]
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + n_SYabc[i] + n_YabYabc[i] + n_YacYabc[i] + n_YbcYabc[i] - n_YabcOut[i]
  update(Z_abc[]) <- Z_abc[i] + n_SZabc[i] + n_ZabZabc[i] + n_ZacZabc[i] + n_ZbcZabc[i] - n_ZabcOut[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  ##STOCHASTIC calculations
  #calculate the probability of leaving S from the rate
  p_SOut[] <- 1 - exp(-sum_all_BI_terms[i])
  
  #calculate the number of people leaving S
  n_SOut[] <- rbinom(S[i], p_SOut[i])
  
  #divide probabilities leaving S 
  p_SYZa[] <- sum_BI_YZa[i]/sum_all_BI_terms[i]
  p_SYZb[] <- sum_BI_YZb[i]/sum_all_BI_terms[i]
  p_SYZc[] <- sum_BI_YZc[i]/sum_all_BI_terms[i]
  p_SYZab[] <- sum_BI_YZab[i]/sum_all_BI_terms[i]
  p_SYZac[] <- sum_BI_YZac[i]/sum_all_BI_terms[i]
  p_SYZbc[] <- sum_BI_YZbc[i]/sum_all_BI_terms[i]
  p_SYZabc[] <- sum_BI_YZabc[i]/sum_all_BI_terms[i] #fixed error here 6/3
  
  #divide the people leaving S into different flows
  n_SYZa[] <- rbinom(n_SOut[i], p_SYZa[i])
  n_SYZb[] <- rbinom(n_SOut[i] - n_SYZa[i], (p_SYZb[i])/(1 - p_SYZa[i] + 1e-8)) #i) subtract those already subtracted from n_SOut ii) renormalize the probability of exit
  n_SYZc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i], (p_SYZc[i])/(1 - p_SYZa[i] - p_SYZb[i] + 1e-8))
  n_SYZab[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i], (p_SYZab[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] + 1e-8))
  n_SYZac[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i], (p_SYZac[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] + 1e-8)) #fixed n to p issue here and remaining 2 eqs 6/5
  n_SYZbc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i], (p_SYZbc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] + 1e-8)) 
  n_SYZabc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i], (p_SYZabc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] - p_SYZbc[i] + 1e-8))
  n_SYZs[] <- n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] - n_SYZabc[i]
  
  #split into symptomatic vs asymptomatics
  n_SYs[] <- rbinom(n_SYZs[i], sigma)
  n_SZs[] <- n_SYZs[i] - n_SYs[i] #fixed issue here
  n_SYa[] <- rbinom(n_SYZa[i], sigma)
  n_SZa[] <- n_SYZa[i] - n_SYa[i]
  n_SYb[] <- rbinom(n_SYZb[i], sigma)
  n_SZb[] <- n_SYZb[i] - n_SYb[i]
  n_SYc[] <- rbinom(n_SYZc[i], sigma)
  n_SZc[] <- n_SYZc[i] - n_SYc[i]
  n_SYab[] <- rbinom(n_SYZab[i], sigma)
  n_SZab[] <- n_SYZab[i] - n_SYab[i]
  n_SYac[] <- rbinom(n_SYZac[i], sigma)
  n_SZac[] <- n_SYZac[i] - n_SYac[i]
  n_SYbc[] <- rbinom(n_SYZbc[i], sigma)
  n_SZbc[] <- n_SYZbc[i] - n_SYbc[i]
  n_SYabc[] <- rbinom(n_SYZabc[i], sigma)
  n_SZabc[] <- n_SYZabc[i] - n_SYabc[i]
  
  #calculate the probability of leaving Ys and Zs by converting exit rates
  p_YsOut <- 1 - exp(-(T_s + d))
  p_ZsOut <- 1 - exp(-(T_m + d))
  
  #calculate the number of people leaving Ys and Zs based on the probabilities
  n_YsOut[] <- rbinom(Y_s[i], p_YsOut)
  n_ZsOut[] <- rbinom(Z_s[i], p_ZsOut)
  
  #split the people recovering naturally (nr) vs the people leaving due to treatment (tr)
  p_Ys_tr <- T_s / (T_s + d)
  n_Ys_tr[] <- rbinom(n_YsOut[i], p_Ys_tr)
  n_Ys_nr[] <- n_YsOut[i] - n_Ys_tr[i]
  p_Zs_tr <- T_m / (T_m + d)
  n_Zs_tr[] <- rbinom(n_ZsOut[i], p_Zs_tr)
  n_Zs_nr[] <- n_ZsOut[i] - n_Zs_tr[i]
  
  #of the people being treated, split the flows into recovery vs. resistance
  p_YZs_tr_fail <- E_a*w_a + E_b*w_b + E_c*w_c
  n_Ys_tr_fail[] <- rbinom(n_Ys_tr[i], p_YZs_tr_fail)
  n_Ys_tr_success[] <- n_Ys_tr[i] - n_Ys_tr_fail[i]
  n_Zs_tr_fail[] <- rbinom(n_Zs_tr[i], p_YZs_tr_fail)
  n_Zs_tr_success[] <- n_Zs_tr[i] - n_Zs_tr_fail[i]
  
  #combine those recovering from successful treatment and natural recovery as leaving Ys and Zs and going back to S
  n_YsS[] <- n_Ys_tr_success[i] + n_Ys_nr[i]
  n_ZsS[] <- n_Zs_tr_success[i] + n_Zs_nr[i]
  
  #split those who are treated and develop resistance into the three possible flows
  #first take the individually normalized probabilities of each flow
  p_YZs_YZa <- E_a*w_a/p_YZs_tr_fail
  p_YZs_YZb <- E_b*w_b/p_YZs_tr_fail
  p_YZs_YZc <- E_c*w_c/p_YZs_tr_fail
  
  #do stick-breaking with Y
  n_YsYa[] <- rbinom(n_Ys_tr_fail[i], p_YZs_YZa)
  n_YsYb[] <- rbinom(n_Ys_tr_fail[i] - n_YsYa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_YsYc[] <- n_Ys_tr_fail[i] - n_YsYa[i] - n_YsYb[i]
  
  #do stick-breaking with Z
  n_ZsZa[] <- rbinom(n_Zs_tr_fail[i], p_YZs_YZa)
  n_ZsZb[] <- rbinom(n_Zs_tr_fail[i] - n_ZsZa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_ZsZc[] <- n_Zs_tr_fail[i] - n_ZsZa[i] - n_ZsZb[i]
  
  ##now onto single resistance compartments. Probabilities of exit will be different for Y and Z compartments due to retreatment
  #probabilities for exit flows from Ya, Yb, Yc, Za, Zb, and Zc
  p_YaOut <- 1 - exp(-((E_b + E_c)*T_s + E_a*kappa*T_sr + d))
  p_YbOut <- 1 - exp(-((E_a + E_c)*T_s + E_b*kappa*T_sr + d))
  p_YcOut <- 1 - exp(-((E_a + E_b)*T_s + E_c*kappa*T_sr + d))
  p_ZaOut <- 1 - exp(-((E_b + E_c)*T_m + d))
  p_ZbOut <- 1 - exp(-((E_a + E_c)*T_m + d))
  p_ZcOut <- 1 - exp(-((E_a + E_b)*T_m + d))
  
  #calculate those leaving the single resistance compartments
  n_YaOut[] <- rbinom(Y_a[i], p_YaOut)
  n_YbOut[] <- rbinom(Y_b[i], p_YbOut)
  n_YcOut[] <- rbinom(Y_c[i], p_YcOut)
  n_ZaOut[] <- rbinom(Z_a[i], p_ZaOut)
  n_ZbOut[] <- rbinom(Z_b[i], p_ZbOut)
  n_ZcOut[] <- rbinom(Z_c[i], p_ZcOut)
  
  #calculate the probabilities of leaving single resistance compartments due to natural recovery
  p_Ya_nr <- d/((E_b + E_c)*T_s + E_a*kappa*T_sr + d)
  p_Yb_nr <- d/((E_a + E_c)*T_s + E_b*kappa*T_sr + d)
  p_Yc_nr <- d/((E_a + E_b)*T_s + E_c*kappa*T_sr + d)
  p_Za_nr <- d/((E_b + E_c)*T_m + d)
  p_Zb_nr <- d/((E_a + E_c)*T_m + d)
  p_Zc_nr <- d/((E_a + E_b)*T_m + d)
  
  #split the people leaving single resistance into natural recovery vs treatment
  n_Ya_nr[] <- rbinom(n_YaOut[i], p_Ya_nr)
  n_Ya_tr[] <- n_YaOut[i] - n_Ya_nr[i]
  n_Yb_nr[] <- rbinom(n_YbOut[i], p_Yb_nr)
  n_Yb_tr[] <- n_YbOut[i] - n_Yb_nr[i]
  n_Yc_nr[] <- rbinom(n_YcOut[i], p_Yc_nr)
  n_Yc_tr[] <- n_YcOut[i] - n_Yc_nr[i]
  n_Za_nr[] <- rbinom(n_ZaOut[i], p_Za_nr)
  n_Za_tr[] <- n_ZaOut[i] - n_Za_nr[i]
  n_Zb_nr[] <- rbinom(n_ZbOut[i], p_Zb_nr)
  n_Zb_tr[] <- n_ZbOut[i] - n_Zb_nr[i]
  n_Zc_nr[] <- rbinom(n_ZcOut[i], p_Zc_nr)
  n_Zc_tr[] <- n_ZcOut[i] - n_Zc_nr[i]
  
  #split those leaving single resistance after treatment into those who fail treatment vs. those who recover
  #note that here, treatment failure means people who develop resistance on treatment, not those treated with the drug they are resistance bc they stay in the resistant compartment
  p_Ya_tr_fail <- (E_b*w_b*T_s + E_c*w_c*T_s)/((E_b + E_c)*T_s + E_a*kappa*T_sr)
  p_Yb_tr_fail <- (E_a*w_a*T_s + E_c*w_c*T_s)/((E_a + E_c)*T_s + E_b*kappa*T_sr)
  p_Yc_tr_fail <- (E_a*w_a*T_s + E_b*w_b*T_s)/((E_a + E_b)*T_s + E_c*kappa*T_sr)
  p_Za_tr_fail <- (E_b*w_b + E_c*w_c)/ ((E_b + E_c) + 1e-8) #add very small term to denominator so this probability is always defined
  p_Zb_tr_fail <- (E_a*w_a + E_c*w_c)/ ((E_a + E_c) + 1e-8)
  p_Zc_tr_fail <- (E_a*w_a + E_b*w_b)/ ((E_a + E_b) + 1e-8)
  
  n_Ya_tr_fail[] <- rbinom(n_Ya_tr[i], p_Ya_tr_fail)
  n_Yb_tr_fail[] <- rbinom(n_Yb_tr[i], p_Yb_tr_fail)
  n_Yc_tr_fail[] <- rbinom(n_Yc_tr[i], p_Yc_tr_fail)
  n_Za_tr_fail[] <- rbinom(n_Za_tr[i], p_Za_tr_fail)
  n_Zb_tr_fail[] <- rbinom(n_Zb_tr[i], p_Zb_tr_fail)
  n_Zc_tr_fail[] <- rbinom(n_Zc_tr[i], p_Zc_tr_fail)
  
  n_Ya_tr_success[] <- n_Ya_tr[i] - n_Ya_tr_fail[i]
  n_Yb_tr_success[] <- n_Yb_tr[i] - n_Yb_tr_fail[i]
  n_Yc_tr_success[] <- n_Yc_tr[i] - n_Yc_tr_fail[i]
  n_Za_tr_success[] <- n_Za_tr[i] - n_Za_tr_fail[i]
  n_Zb_tr_success[] <- n_Zb_tr[i] - n_Zb_tr_fail[i]
  n_Zc_tr_success[] <- n_Zc_tr[i] - n_Zc_tr_fail[i]
  
  #combine treatment success and natural recover terms to get single resistance recovery flows
  n_YaS[] <- n_Ya_tr_success[i] + n_Ya_nr[i]
  n_YbS[] <- n_Yb_tr_success[i] + n_Yb_nr[i]
  n_YcS[] <- n_Yc_tr_success[i] + n_Yc_nr[i]
  n_ZaS[] <- n_Za_tr_success[i] + n_Za_nr[i]
  n_ZbS[] <- n_Zb_tr_success[i] + n_Zb_nr[i]
  n_ZcS[] <- n_Zc_tr_success[i] + n_Zc_nr[i]
  
  #now split the treatment failures into the two different resistance compartment flows
  p_YZaYZab <- E_b*w_b / (E_b*w_b + E_c*w_c + 1e-8) #add very small term to denominator so this probability is always defined
  p_YZbYZab <- E_a*w_a / (E_a*w_a + E_c*w_c + 1e-8)
  p_YZcYZac <- E_a*w_a / (E_a*w_a + E_b*w_b + 1e-8)
  
  n_YaYab[] <- rbinom(n_Ya_tr_fail[i], p_YZaYZab)
  n_YaYac[] <- n_Ya_tr_fail[i] - n_YaYab[i]
  n_YbYab[] <- rbinom(n_Yb_tr_fail[i], p_YZbYZab)
  n_YbYbc[] <- n_Yb_tr_fail[i] - n_YbYab[i]
  n_YcYac[] <- rbinom(n_Yc_tr_fail[i], p_YZcYZac)
  n_YcYbc[] <- n_Yc_tr_fail[i] - n_YcYac[i]
  
  n_ZaZab[] <- rbinom(n_Za_tr_fail[i], p_YZaYZab)
  n_ZaZac[] <- n_Za_tr_fail[i] - n_ZaZab[i]
  n_ZbZab[] <- rbinom(n_Zb_tr_fail[i], p_YZbYZab)
  n_ZbZbc[] <- n_Zb_tr_fail[i] - n_ZbZab[i]
  n_ZcZac[] <- rbinom(n_Zc_tr_fail[i], p_YZcYZac)
  n_ZcZbc[] <- n_Zc_tr_fail[i] - n_ZcZac[i]
  
  ##now move to the double resistance compartments
  #first calculate total outflow -- calculate exit probabilities and then draw
  p_YabOut <- 1 - exp(-(E_c*T_s + (1 - E_c)*kappa*T_sr + d))
  p_YacOut <- 1 - exp(-(E_b*T_s + (1 - E_b)*kappa*T_sr + d))
  p_YbcOut <- 1 - exp(-(E_a*T_s + (1 - E_a)*kappa*T_sr + d))
  p_ZabOut <- 1 - exp(-(E_c*T_m + d))
  p_ZacOut <- 1 - exp(-(E_b*T_m + d))
  p_ZbcOut <- 1 - exp(-(E_a*T_m + d))
  
  n_YabOut[] <- rbinom(Y_ab[i], p_YabOut)
  n_YacOut[] <- rbinom(Y_ac[i], p_YacOut)
  n_YbcOut[] <- rbinom(Y_bc[i], p_YbcOut)
  n_ZabOut[] <- rbinom(Z_ab[i], p_ZabOut)
  n_ZacOut[] <- rbinom(Z_ac[i], p_ZacOut)
  n_ZbcOut[] <- rbinom(Z_bc[i], p_ZbcOut)
  
  #because these are already dual resistance, the only ones leaving and not recovering are developing resistance to the same drug, so separate those and the rest recover
  p_Yab_tr_fail <- (E_c*w_c*T_s)/(E_c*T_s + (1 - E_c)*kappa*T_sr + d)
  p_Yac_tr_fail <- (E_b*w_b*T_s)/(E_b*T_s + (1 - E_b)*kappa*T_sr + d)
  p_Ybc_tr_fail <- (E_a*w_a*T_s)/(E_a*T_s + (1 - E_a)*kappa*T_sr + d)
  p_Zab_tr_fail <- (E_c*w_c*T_m)/(E_c*T_m + d)
  p_Zac_tr_fail <- (E_b*w_b*T_m)/(E_b*T_m + d)
  p_Zbc_tr_fail <- (E_a*w_a*T_m)/(E_a*T_m + d)
  
  #split out resistant flows to Yabc vs recoveries to S
  n_YabYabc[] <- rbinom(n_YabOut[i], p_Yab_tr_fail)
  n_YabS[] <- n_YabOut[i] - n_YabYabc[i]
  n_YacYabc[] <- rbinom(n_YacOut[i], p_Yac_tr_fail)
  n_YacS[] <- n_YacOut[i] - n_YacYabc[i]
  n_YbcYabc[] <- rbinom(n_YbcOut[i], p_Ybc_tr_fail)
  n_YbcS[] <- n_YbcOut[i] - n_YbcYabc[i]
  
  n_ZabZabc[] <- rbinom(n_ZabOut[i], p_Zab_tr_fail)
  n_ZabS[] <- n_ZabOut[i] - n_ZabZabc[i]
  n_ZacZabc[] <- rbinom(n_ZacOut[i], p_Zac_tr_fail)
  n_ZacS[] <- n_ZacOut[i] - n_ZacZabc[i]
  n_ZbcZabc[] <- rbinom(n_ZbcOut[i], p_Zbc_tr_fail)
  n_ZbcS[] <- n_ZbcOut[i] - n_ZbcZabc[i]
  
  #recovery from triple resistance states
  p_YabcOut <- 1 - exp(-(kappa*T_sr + d))
  p_ZabcOut <- 1 - exp(-d)
  n_YabcOut[] <- rbinom(Y_abc[i], p_YabcOut)
  n_ZabcOut[] <- rbinom(Z_abc[i], p_ZabcOut)
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #Binary switches between treatment
  #THIS PART IS UNIQUE TO THE MODEL STRATEGY
  E_a <- (lock_threshold_A < 1)
  E_c <- 1 - (lock_threshold_B < 1)
  E_b <- 1  - E_a - E_c
  
  sum_tx_probs = E_a + E_c + E_b 
  
  #sum array variables for book-keeping / accounting outputs
  sum_n_SOut <- sum(n_SOut)
  sum_n_SYZa <- sum(n_SYZa)
  sum_n_SYZb <- sum(n_SYZb)
  sum_n_SYZc <- sum(n_SYZc)
  sum_n_SYZab <- sum(n_SYZab)
  sum_n_SYZac <- sum(n_SYZac)
  sum_n_SYZbc <- sum(n_SYZbc)
  sum_n_SYZabc <- sum(n_SYZabc)
  sum_n_SYZs <- sum(n_SYZs)
  sum_n_SYs <- sum(n_SYs)
  sum_n_SZs <- sum(n_SZs)
  sum_n_SYa <- sum(n_SYa)
  sum_n_SZa <- sum(n_SZa)
  sum_n_SYb <- sum(n_SYb)
  sum_n_SZb <- sum(n_SZb)
  sum_n_SYc <- sum(n_SYc)
  sum_n_SZc <- sum(n_SZc)
  sum_n_SYab <- sum(n_SYab)
  sum_n_SZab <- sum(n_SZab)
  sum_n_SYac <- sum(n_SYac)
  sum_n_SZac <- sum(n_SZac)
  sum_n_SYbc <- sum(n_SYbc)
  sum_n_SZbc <- sum(n_SZbc)
  sum_n_SYabc <- sum(n_SYabc)
  sum_n_SZabc <- sum(n_SZabc)
  sum_n_YsOut <- sum(n_YsOut)
  sum_n_ZsOut <- sum(n_ZsOut)
  sum_n_Ys_tr <- sum(n_Ys_tr)
  sum_n_Ys_nr <- sum(n_Ys_nr)
  sum_n_Zs_tr <- sum(n_Zs_tr)
  sum_n_Zs_nr <- sum(n_Zs_nr)
  sum_n_Ys_tr_fail <- sum(n_Ys_tr_fail)
  sum_n_Ys_tr_success <- sum(n_Ys_tr_success)
  sum_n_Zs_tr_fail <- sum(n_Zs_tr_fail)
  sum_n_Zs_tr_success <- sum(n_Zs_tr_success)
  sum_n_YsS <- sum(n_YsS)
  sum_n_ZsS <- sum(n_ZsS)
  sum_n_YsYa <- sum(n_YsYa)
  sum_n_YsYb <- sum(n_YsYb)
  sum_n_YsYc <- sum(n_YsYc)
  sum_n_ZsZa <- sum(n_ZsZa)
  sum_n_ZsZb <- sum(n_ZsZb)
  sum_n_ZsZc <- sum(n_ZsZc)
  sum_n_YaOut <- sum(n_YaOut)
  sum_n_YbOut <- sum(n_YbOut)
  sum_n_YcOut <- sum(n_YcOut)
  sum_n_ZaOut <- sum(n_ZaOut)
  sum_n_ZbOut <- sum(n_ZbOut)
  sum_n_ZcOut <- sum(n_ZcOut)
  sum_n_Ya_nr <- sum(n_Ya_nr)
  sum_n_Ya_tr <- sum(n_Ya_tr)
  sum_n_Yb_nr <- sum(n_Yb_nr)
  sum_n_Yb_tr <- sum(n_Yb_tr)
  sum_n_Yc_nr <- sum(n_Yc_nr)
  sum_n_Yc_tr <- sum(n_Yc_tr)
  sum_n_Za_nr <- sum(n_Za_nr)
  sum_n_Za_tr <- sum(n_Za_tr)
  sum_n_Zb_nr <- sum(n_Zb_nr)
  sum_n_Zb_tr <- sum(n_Zb_tr)
  sum_n_Zc_nr <- sum(n_Zc_nr)
  sum_n_Zc_tr <- sum(n_Zc_tr)
  sum_n_Ya_tr_fail <- sum(n_Ya_tr_fail)
  sum_n_Yb_tr_fail <- sum(n_Yb_tr_fail)
  sum_n_Yc_tr_fail <- sum(n_Yc_tr_fail)
  sum_n_Za_tr_fail <- sum(n_Za_tr_fail)
  sum_n_Zb_tr_fail <- sum(n_Zb_tr_fail)
  sum_n_Zc_tr_fail <- sum(n_Zc_tr_fail)
  sum_n_Ya_tr_success <- sum(n_Ya_tr_success)
  sum_n_Yb_tr_success <- sum(n_Yb_tr_success)
  sum_n_Yc_tr_success <- sum(n_Yc_tr_success)
  sum_n_Za_tr_success <- sum(n_Za_tr_success)
  sum_n_Zb_tr_success <- sum(n_Zb_tr_success)
  sum_n_Zc_tr_success <- sum(n_Zc_tr_success)
  sum_n_YaS <- sum(n_YaS)
  sum_n_YbS <- sum(n_YbS)
  sum_n_YcS <- sum(n_YcS)
  sum_n_ZaS <- sum(n_ZaS)
  sum_n_ZbS <- sum(n_ZbS)
  sum_n_ZcS <- sum(n_ZcS)
  sum_n_YaYab <- sum(n_YaYab)
  sum_n_YaYac <- sum(n_YaYac)
  sum_n_YbYab <- sum(n_YbYab)
  sum_n_YbYbc <- sum(n_YbYbc)
  sum_n_YcYac <- sum(n_YcYac)
  sum_n_YcYbc <- sum(n_YcYbc)
  sum_n_ZaZab <- sum(n_ZaZab)
  sum_n_ZaZac <- sum(n_ZaZac)
  sum_n_ZbZab <- sum(n_ZbZab)
  sum_n_ZbZbc <- sum(n_ZbZbc)
  sum_n_ZcZac <- sum(n_ZcZac)
  sum_n_ZcZbc <- sum(n_ZcZbc)
  sum_n_YabOut <- sum(n_YabOut)
  sum_n_YacOut <- sum(n_YacOut)
  sum_n_YbcOut <- sum(n_YbcOut)
  sum_n_ZabOut <- sum(n_ZabOut)
  sum_n_ZacOut <- sum(n_ZacOut)
  sum_n_ZbcOut <- sum(n_ZbcOut)
  sum_n_YabYabc <- sum(n_YabYabc)
  sum_n_YabS <- sum(n_YabS)
  sum_n_YacYabc <- sum(n_YacYabc)
  sum_n_YacS <- sum(n_YacS)
  sum_n_YbcYabc <- sum(n_YbcYabc)
  sum_n_YbcS<- sum(n_YbcS)
  sum_n_ZabZabc <- sum(n_ZabZabc)
  sum_n_ZabS <- sum(n_ZabS)
  sum_n_ZacZabc <- sum(n_ZacZabc)
  sum_n_ZacS <- sum(n_ZacS)
  sum_n_ZbcZabc <- sum(n_ZbcZabc)
  sum_n_ZbcS<- sum(n_ZbcS)
  sum_n_YabcOut <- sum(n_YabcOut)
  sum_n_ZabcOut <- sum(n_ZabcOut)
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  dim(n_SOut) <- N_risk
  dim(p_SOut) <- N_risk
  dim(p_SYZa) <- N_risk
  dim(p_SYZb) <- N_risk
  dim(p_SYZc) <- N_risk
  dim(p_SYZab) <- N_risk
  dim(p_SYZac) <- N_risk  
  dim(p_SYZbc) <- N_risk
  dim(p_SYZabc) <- N_risk
  dim(n_SYZa) <- N_risk
  dim(n_SYZb) <- N_risk
  dim(n_SYZc) <- N_risk
  dim(n_SYZab) <- N_risk
  dim(n_SYZac) <- N_risk
  dim(n_SYZbc) <- N_risk
  dim(n_SYZabc) <- N_risk
  dim(n_SYZs) <- N_risk
  dim(n_SYs) <- N_risk
  dim(n_SZs) <- N_risk
  dim(n_SYa) <- N_risk
  dim(n_SZa) <- N_risk
  dim(n_SYb) <- N_risk
  dim(n_SZb) <- N_risk
  dim(n_SYc) <- N_risk
  dim(n_SZc) <- N_risk
  dim(n_SYab) <- N_risk
  dim(n_SZab) <- N_risk
  dim(n_SYac) <- N_risk
  dim(n_SZac) <- N_risk
  dim(n_SYbc) <- N_risk
  dim(n_SZbc) <- N_risk
  dim(n_SYabc) <- N_risk
  dim(n_SZabc) <- N_risk
  dim(n_YsOut) <- N_risk
  dim(n_ZsOut) <- N_risk
  dim(n_Ys_tr) <- N_risk
  dim(n_Ys_nr) <- N_risk
  dim(n_Zs_tr) <- N_risk
  dim(n_Zs_nr) <- N_risk
  dim(n_Ys_tr_fail) <- N_risk
  dim(n_Ys_tr_success) <- N_risk
  dim(n_Zs_tr_fail) <- N_risk
  dim(n_Zs_tr_success) <- N_risk
  dim(n_YsS) <- N_risk
  dim(n_ZsS) <- N_risk
  dim(n_YsYa) <- N_risk
  dim(n_YsYb) <- N_risk
  dim(n_YsYc) <- N_risk
  dim(n_ZsZa) <- N_risk
  dim(n_ZsZb) <- N_risk
  dim(n_ZsZc) <- N_risk
  dim(n_YaOut) <- N_risk
  dim(n_YbOut) <- N_risk
  dim(n_YcOut) <- N_risk
  dim(n_ZaOut) <- N_risk
  dim(n_ZbOut) <- N_risk
  dim(n_ZcOut) <- N_risk
  dim(n_Ya_nr) <- N_risk
  dim(n_Ya_tr) <- N_risk
  dim(n_Yb_nr) <- N_risk
  dim(n_Yb_tr) <- N_risk
  dim(n_Yc_nr) <- N_risk
  dim(n_Yc_tr) <- N_risk
  dim(n_Za_nr) <- N_risk
  dim(n_Za_tr) <- N_risk
  dim(n_Zb_nr) <- N_risk
  dim(n_Zb_tr) <- N_risk
  dim(n_Zc_nr) <- N_risk
  dim(n_Zc_tr) <- N_risk
  dim(n_Ya_tr_fail) <- N_risk
  dim(n_Yb_tr_fail) <- N_risk
  dim(n_Yc_tr_fail) <- N_risk
  dim(n_Za_tr_fail) <- N_risk
  dim(n_Zb_tr_fail) <- N_risk
  dim(n_Zc_tr_fail) <- N_risk
  dim(n_Ya_tr_success) <- N_risk
  dim(n_Yb_tr_success) <- N_risk
  dim(n_Yc_tr_success) <- N_risk
  dim(n_Za_tr_success) <- N_risk
  dim(n_Zb_tr_success) <- N_risk
  dim(n_Zc_tr_success) <- N_risk
  dim(n_YaS) <- N_risk
  dim(n_YbS) <- N_risk
  dim(n_YcS) <- N_risk
  dim(n_ZaS) <- N_risk
  dim(n_ZbS) <- N_risk
  dim(n_ZcS) <- N_risk
  dim(n_YaYab) <- N_risk
  dim(n_YaYac) <- N_risk
  dim(n_YbYab) <- N_risk
  dim(n_YbYbc) <- N_risk
  dim(n_YcYac) <- N_risk
  dim(n_YcYbc) <- N_risk
  dim(n_ZaZab) <- N_risk
  dim(n_ZaZac) <- N_risk
  dim(n_ZbZab) <- N_risk
  dim(n_ZbZbc) <- N_risk
  dim(n_ZcZac) <- N_risk
  dim(n_ZcZbc) <- N_risk
  dim(n_YabOut) <- N_risk
  dim(n_YacOut) <- N_risk
  dim(n_YbcOut) <- N_risk
  dim(n_ZabOut) <- N_risk
  dim(n_ZacOut) <- N_risk
  dim(n_ZbcOut) <- N_risk
  dim(n_YabYabc) <- N_risk
  dim(n_YabS) <- N_risk
  dim(n_YacYabc) <- N_risk
  dim(n_YacS) <- N_risk
  dim(n_YbcYabc) <- N_risk
  dim(n_YbcS) <- N_risk
  dim(n_ZabZabc) <- N_risk
  dim(n_ZabS) <- N_risk
  dim(n_ZacZabc) <- N_risk
  dim(n_ZacS) <- N_risk
  dim(n_ZbcZabc) <- N_risk
  dim(n_ZbcS) <- N_risk
  dim(n_YabcOut) <- N_risk
  dim(n_ZabcOut) <- N_risk
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  
  #define desired non-state outputs
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev

})

equal_allocation_stochastic_6_7_brief <- odin({
  ##write out all intermediate transmission terms needed for equations
  #write out intermediate equation for beta*I for the dS equation
  BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
                          f_a*(Y_a[i] + Z_a[i]) +
                          f_b*(Y_b[i] + Z_b[i]) +
                          f_c*(Y_c[i] + Z_c[i]) + 
                          f_ab*(Y_ab[i] + Z_ab[i]) + 
                          f_ac*(Y_ac[i] + Z_ac[i]) + 
                          f_bc*(Y_bc[i] + Z_bc[i]) + 
                          f_abc*(Y_abc[i] + Z_abc[i]))
  
  #write out intermediate equation for beta*I for the dY_s and dZ_s equations
  BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
  
  #write out intermediate equation for beta*I for the dY_a and dZ_a equations
  BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
  
  #write out intermediate equation for beta*I for the dY_b and dZ_b equations
  BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
  
  #write out intermediate equation for beta*I for the dY_c and dZ_c equations
  BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
  
  #write out intermediate equation for beta*I for dY_ab and dZ_ab
  BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
  
  #write out intermediate equation for beta*I for dY_ac and dZ_ac
  BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
  
  #write out intermediate equation for beta*I for dY_bc and dZ_bc
  BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
  
  #write out intermediate equation for dY_abc and dZ_abc
  BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
  
  ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
  sum_BI_S[] <- sum(BI_S[,i])
  sum_BI_YZs[] <- sum(BI_YZs[,i])
  sum_BI_YZa[] <- sum(BI_YZa[,i])
  sum_BI_YZb[] <- sum(BI_YZb[,i])
  sum_BI_YZc[] <- sum(BI_YZc[,i])
  sum_BI_YZab[] <- sum(BI_YZab[,i])
  sum_BI_YZac[] <- sum(BI_YZac[,i])
  sum_BI_YZbc[] <- sum(BI_YZbc[,i])
  sum_BI_YZabc[] <- sum(BI_YZabc[,i])
  #sum all BI terms to get the total efflux from S
  sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
    sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
  
  ##total across risk strata of BI terms for bookkeeping (to output)
  total_sum_BI_S <- sum(sum_BI_S)
  total_sum_BI_YZs <- sum(sum_BI_YZs)
  total_sum_BI_YZa <- sum(sum_BI_YZa)
  total_sum_BI_YZb <- sum(sum_BI_YZb)
  total_sum_BI_YZc <- sum(sum_BI_YZc)
  total_sum_BI_YZab <- sum(sum_BI_YZab)
  total_sum_BI_YZac <- sum(sum_BI_YZac)
  total_sum_BI_YZbc <- sum(sum_BI_YZbc)
  total_sum_BI_YZabc <- sum(sum_BI_YZabc)
  total_sum_all_BI_terms <- sum(sum_all_BI_terms)
  
  #write out equation for S
  update(S[]) <- S[i] - n_SOut[i] + n_YsS[i] + n_YaS[i] + n_YbS[i] + n_YcS[i] + n_YabS[i] + n_YacS[i] + n_YbcS[i] +
    n_YabcOut[i] + n_ZsS[i] + n_ZaS[i] + n_ZbS[i] + n_ZcS[i] + n_ZabS[i] + n_ZacS[i] + n_ZbcS[i] + n_ZabcOut[i] 
  
  #write out equations for dY_s and dZ_s
  update(Y_s[]) <- Y_s[i] + n_SYs[i] - n_YsOut[i]
  update(Z_s[]) <- Z_s[i] + n_SZs[i] - n_ZsOut[i]
  
  #write out equations for dY_a and dZ_a
  update(Y_a[]) <- Y_a[i] + n_SYa[i] + n_YsYa[i] - n_YaOut[i]
  update(Z_a[]) <- Z_a[i] + n_SZa[i] + n_ZsZa[i] - n_ZaOut[i]
  
  #write out equations for dY_b and dZ_b
  update(Y_b[]) <- Y_b[i] + n_SYb[i] + n_YsYb[i] - n_YbOut[i]
  update(Z_b[]) <- Z_b[i] + n_SZb[i] + n_ZsZb[i] - n_ZbOut[i]
  
  #write out equations for dY_c and dZ_c
  update(Y_c[]) <- Y_c[i] + n_SYc[i] + n_YsYc[i] - n_YcOut[i] 
  update(Z_c[]) <- Z_c[i] + n_SZc[i] + n_ZsZc[i] - n_ZcOut[i]
  
  #write out equations for dY_ab and dZ_ab
  update(Y_ab[]) <- Y_ab[i] + n_SYab[i] + n_YaYab[i] + n_YbYab[i] - n_YabOut[i]
  update(Z_ab[]) <- Z_ab[i] + n_SZab[i] + n_ZaZab[i] + n_ZbZab[i] - n_ZabOut[i]
  
  #write out equations for dY_ac and dZ_ac
  update(Y_ac[]) <- Y_ac[i] + n_SYac[i] + n_YaYac[i] + n_YcYac[i] - n_YacOut[i]
  update(Z_ac[]) <- Z_ac[i] + n_SZac[i] + n_ZaZac[i] + n_ZcZac[i] - n_ZacOut[i] 
  
  #write out equations for dY_ac and dZ_ac
  update(Y_bc[]) <- Y_bc[i] + n_SYbc[i] + n_YbYbc[i] + n_YcYbc[i] - n_YbcOut[i]
  update(Z_bc[]) <- Z_bc[i] + n_SZbc[i] + n_ZbZbc[i] + n_ZcZbc[i] - n_ZbcOut[i]
  
  #write out equations for dY_abc and dZ_abc
  update(Y_abc[]) <- Y_abc[i] + n_SYabc[i] + n_YabYabc[i] + n_YacYabc[i] + n_YbcYabc[i] - n_YabcOut[i]
  update(Z_abc[]) <- Z_abc[i] + n_SZabc[i] + n_ZabZabc[i] + n_ZacZabc[i] + n_ZbcZabc[i] - n_ZabcOut[i]
  
  #equations for non-state outputs
  N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
    Z_bc[i] + Y_abc[i] + Z_abc[i]
  
  sum_N <- sum(N)
  
  #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
  prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
  
  prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
             sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
  
  all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
  all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
  all_non_res <- sum(Y_s) + sum(Z_s)
  
  #sum states at any point in time across risk groups for book-keeping
  sum_S <- sum(S)
  sum_Y_s <- sum(Y_s)
  sum_Y_a <- sum(Y_a)
  sum_Y_b <- sum(Y_b)
  sum_Y_c <- sum(Y_c)
  sum_Y_ab <- sum(Y_ab)
  sum_Y_ac <- sum(Y_ac)
  sum_Y_bc <- sum(Y_bc)
  sum_Y_abc <- sum(Y_abc)
  sum_Z_s <- sum(Z_s)
  sum_Z_a <- sum(Z_a)
  sum_Z_b <- sum(Z_b)
  sum_Z_c <- sum(Z_c)
  sum_Z_ab <- sum(Z_ab)
  sum_Z_ac <- sum(Z_ac)
  sum_Z_bc <- sum(Z_bc)
  sum_Z_abc <- sum(Z_abc)
  
  #sum states across symptom status for book-keeping
  sum_YZ_s <- sum_Y_s + sum_Z_s
  sum_YZ_a <- sum_Y_a + sum_Z_a
  sum_YZ_b <- sum_Y_b + sum_Z_b
  sum_YZ_c <- sum_Y_c + sum_Z_c 
  sum_YZ_ab <- sum_Y_ab + sum_Z_ab
  sum_YZ_ac <- sum_Y_ac + sum_Z_ac
  sum_YZ_bc <- sum_Y_bc + sum_Z_bc
  sum_YZ_abc <- sum_Y_abc + sum_Z_abc
  
  ##STOCHASTIC calculations
  #calculate the probability of leaving S from the rate
  p_SOut[] <- 1 - exp(-sum_all_BI_terms[i])
  
  #calculate the number of people leaving S
  n_SOut[] <- rbinom(S[i], p_SOut[i])
  
  #divide probabilities leaving S 
  p_SYZa[] <- sum_BI_YZa[i]/sum_all_BI_terms[i]
  p_SYZb[] <- sum_BI_YZb[i]/sum_all_BI_terms[i]
  p_SYZc[] <- sum_BI_YZc[i]/sum_all_BI_terms[i]
  p_SYZab[] <- sum_BI_YZab[i]/sum_all_BI_terms[i]
  p_SYZac[] <- sum_BI_YZac[i]/sum_all_BI_terms[i]
  p_SYZbc[] <- sum_BI_YZbc[i]/sum_all_BI_terms[i]
  p_SYZabc[] <- sum_BI_YZabc[i]/sum_all_BI_terms[i] #fixed error here 6/3
  
  #divide the people leaving S into different flows
  n_SYZa[] <- rbinom(n_SOut[i], p_SYZa[i])
  n_SYZb[] <- rbinom(n_SOut[i] - n_SYZa[i], (p_SYZb[i])/(1 - p_SYZa[i] + 1e-8)) #i) subtract those already subtracted from n_SOut ii) renormalize the probability of exit
  n_SYZc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i], (p_SYZc[i])/(1 - p_SYZa[i] - p_SYZb[i] + 1e-8))
  n_SYZab[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i], (p_SYZab[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] + 1e-8))
  n_SYZac[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i], (p_SYZac[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] + 1e-8)) #fixed n to p issue here and remaining 2 eqs 6/5
  n_SYZbc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i], (p_SYZbc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] + 1e-8)) 
  n_SYZabc[] <- rbinom(n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i], (p_SYZabc[i])/(1 - p_SYZa[i] - p_SYZb[i] - p_SYZc[i] - p_SYZab[i] - p_SYZac[i] - p_SYZbc[i] + 1e-8))
  n_SYZs[] <- n_SOut[i] - n_SYZa[i] - n_SYZb[i] - n_SYZc[i] - n_SYZab[i] - n_SYZac[i] - n_SYZbc[i] - n_SYZabc[i]
  
  #split into symptomatic vs asymptomatics
  n_SYs[] <- rbinom(n_SYZs[i], sigma)
  n_SZs[] <- n_SYZs[i] - n_SYs[i] #fixed issue here
  n_SYa[] <- rbinom(n_SYZa[i], sigma)
  n_SZa[] <- n_SYZa[i] - n_SYa[i]
  n_SYb[] <- rbinom(n_SYZb[i], sigma)
  n_SZb[] <- n_SYZb[i] - n_SYb[i]
  n_SYc[] <- rbinom(n_SYZc[i], sigma)
  n_SZc[] <- n_SYZc[i] - n_SYc[i]
  n_SYab[] <- rbinom(n_SYZab[i], sigma)
  n_SZab[] <- n_SYZab[i] - n_SYab[i]
  n_SYac[] <- rbinom(n_SYZac[i], sigma)
  n_SZac[] <- n_SYZac[i] - n_SYac[i]
  n_SYbc[] <- rbinom(n_SYZbc[i], sigma)
  n_SZbc[] <- n_SYZbc[i] - n_SYbc[i]
  n_SYabc[] <- rbinom(n_SYZabc[i], sigma)
  n_SZabc[] <- n_SYZabc[i] - n_SYabc[i]
  
  #calculate the probability of leaving Ys and Zs by converting exit rates
  p_YsOut <- 1 - exp(-(T_s + d))
  p_ZsOut <- 1 - exp(-(T_m + d))
  
  #calculate the number of people leaving Ys and Zs based on the probabilities
  n_YsOut[] <- rbinom(Y_s[i], p_YsOut)
  n_ZsOut[] <- rbinom(Z_s[i], p_ZsOut)
  
  #split the people recovering naturally (nr) vs the people leaving due to treatment (tr)
  p_Ys_tr <- T_s / (T_s + d)
  n_Ys_tr[] <- rbinom(n_YsOut[i], p_Ys_tr)
  n_Ys_nr[] <- n_YsOut[i] - n_Ys_tr[i]
  p_Zs_tr <- T_m / (T_m + d)
  n_Zs_tr[] <- rbinom(n_ZsOut[i], p_Zs_tr)
  n_Zs_nr[] <- n_ZsOut[i] - n_Zs_tr[i]
  
  #of the people being treated, split the flows into recovery vs. resistance
  p_YZs_tr_fail <- E_a*w_a + E_b*w_b + E_c*w_c
  n_Ys_tr_fail[] <- rbinom(n_Ys_tr[i], p_YZs_tr_fail)
  n_Ys_tr_success[] <- n_Ys_tr[i] - n_Ys_tr_fail[i]
  n_Zs_tr_fail[] <- rbinom(n_Zs_tr[i], p_YZs_tr_fail)
  n_Zs_tr_success[] <- n_Zs_tr[i] - n_Zs_tr_fail[i]
  
  #combine those recovering from successful treatment and natural recovery as leaving Ys and Zs and going back to S
  n_YsS[] <- n_Ys_tr_success[i] + n_Ys_nr[i]
  n_ZsS[] <- n_Zs_tr_success[i] + n_Zs_nr[i]
  
  #split those who are treated and develop resistance into the three possible flows
  #first take the individually normalized probabilities of each flow
  p_YZs_YZa <- E_a*w_a/p_YZs_tr_fail
  p_YZs_YZb <- E_b*w_b/p_YZs_tr_fail
  p_YZs_YZc <- E_c*w_c/p_YZs_tr_fail
  
  #do stick-breaking with Y
  n_YsYa[] <- rbinom(n_Ys_tr_fail[i], p_YZs_YZa)
  n_YsYb[] <- rbinom(n_Ys_tr_fail[i] - n_YsYa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_YsYc[] <- n_Ys_tr_fail[i] - n_YsYa[i] - n_YsYb[i]
  
  #do stick-breaking with Z
  n_ZsZa[] <- rbinom(n_Zs_tr_fail[i], p_YZs_YZa)
  n_ZsZb[] <- rbinom(n_Zs_tr_fail[i] - n_ZsZa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
  n_ZsZc[] <- n_Zs_tr_fail[i] - n_ZsZa[i] - n_ZsZb[i]
  
  ##now onto single resistance compartments. Probabilities of exit will be different for Y and Z compartments due to retreatment
  #probabilities for exit flows from Ya, Yb, Yc, Za, Zb, and Zc
  p_YaOut <- 1 - exp(-((E_b + E_c)*T_s + E_a*kappa*T_sr + d))
  p_YbOut <- 1 - exp(-((E_a + E_c)*T_s + E_b*kappa*T_sr + d))
  p_YcOut <- 1 - exp(-((E_a + E_b)*T_s + E_c*kappa*T_sr + d))
  p_ZaOut <- 1 - exp(-((E_b + E_c)*T_m + d))
  p_ZbOut <- 1 - exp(-((E_a + E_c)*T_m + d))
  p_ZcOut <- 1 - exp(-((E_a + E_b)*T_m + d))
  
  #calculate those leaving the single resistance compartments
  n_YaOut[] <- rbinom(Y_a[i], p_YaOut)
  n_YbOut[] <- rbinom(Y_b[i], p_YbOut)
  n_YcOut[] <- rbinom(Y_c[i], p_YcOut)
  n_ZaOut[] <- rbinom(Z_a[i], p_ZaOut)
  n_ZbOut[] <- rbinom(Z_b[i], p_ZbOut)
  n_ZcOut[] <- rbinom(Z_c[i], p_ZcOut)
  
  #calculate the probabilities of leaving single resistance compartments due to natural recovery
  p_Ya_nr <- d/((E_b + E_c)*T_s + E_a*kappa*T_sr + d)
  p_Yb_nr <- d/((E_a + E_c)*T_s + E_b*kappa*T_sr + d)
  p_Yc_nr <- d/((E_a + E_b)*T_s + E_c*kappa*T_sr + d)
  p_Za_nr <- d/((E_b + E_c)*T_m + d)
  p_Zb_nr <- d/((E_a + E_c)*T_m + d)
  p_Zc_nr <- d/((E_a + E_b)*T_m + d)
  
  #split the people leaving single resistance into natural recovery vs treatment
  n_Ya_nr[] <- rbinom(n_YaOut[i], p_Ya_nr)
  n_Ya_tr[] <- n_YaOut[i] - n_Ya_nr[i]
  n_Yb_nr[] <- rbinom(n_YbOut[i], p_Yb_nr)
  n_Yb_tr[] <- n_YbOut[i] - n_Yb_nr[i]
  n_Yc_nr[] <- rbinom(n_YcOut[i], p_Yc_nr)
  n_Yc_tr[] <- n_YcOut[i] - n_Yc_nr[i]
  n_Za_nr[] <- rbinom(n_ZaOut[i], p_Za_nr)
  n_Za_tr[] <- n_ZaOut[i] - n_Za_nr[i]
  n_Zb_nr[] <- rbinom(n_ZbOut[i], p_Zb_nr)
  n_Zb_tr[] <- n_ZbOut[i] - n_Zb_nr[i]
  n_Zc_nr[] <- rbinom(n_ZcOut[i], p_Zc_nr)
  n_Zc_tr[] <- n_ZcOut[i] - n_Zc_nr[i]
  
  #split those leaving single resistance after treatment into those who fail treatment vs. those who recover
  #note that here, treatment failure means people who develop resistance on treatment, not those treated with the drug they are resistance bc they stay in the resistant compartment
  p_Ya_tr_fail <- (E_b*w_b*T_s + E_c*w_c*T_s)/((E_b + E_c)*T_s + E_a*kappa*T_sr)
  p_Yb_tr_fail <- (E_a*w_a*T_s + E_c*w_c*T_s)/((E_a + E_c)*T_s + E_b*kappa*T_sr)
  p_Yc_tr_fail <- (E_a*w_a*T_s + E_b*w_b*T_s)/((E_a + E_b)*T_s + E_c*kappa*T_sr)
  p_Za_tr_fail <- (E_b*w_b + E_c*w_c)/ ((E_b + E_c) + 1e-8) #add very small term to denominator so this probability is always defined
  p_Zb_tr_fail <- (E_a*w_a + E_c*w_c)/ ((E_a + E_c) + 1e-8)
  p_Zc_tr_fail <- (E_a*w_a + E_b*w_b)/ ((E_a + E_b) + 1e-8)
  
  n_Ya_tr_fail[] <- rbinom(n_Ya_tr[i], p_Ya_tr_fail)
  n_Yb_tr_fail[] <- rbinom(n_Yb_tr[i], p_Yb_tr_fail)
  n_Yc_tr_fail[] <- rbinom(n_Yc_tr[i], p_Yc_tr_fail)
  n_Za_tr_fail[] <- rbinom(n_Za_tr[i], p_Za_tr_fail)
  n_Zb_tr_fail[] <- rbinom(n_Zb_tr[i], p_Zb_tr_fail)
  n_Zc_tr_fail[] <- rbinom(n_Zc_tr[i], p_Zc_tr_fail)
  
  n_Ya_tr_success[] <- n_Ya_tr[i] - n_Ya_tr_fail[i]
  n_Yb_tr_success[] <- n_Yb_tr[i] - n_Yb_tr_fail[i]
  n_Yc_tr_success[] <- n_Yc_tr[i] - n_Yc_tr_fail[i]
  n_Za_tr_success[] <- n_Za_tr[i] - n_Za_tr_fail[i]
  n_Zb_tr_success[] <- n_Zb_tr[i] - n_Zb_tr_fail[i]
  n_Zc_tr_success[] <- n_Zc_tr[i] - n_Zc_tr_fail[i]
  
  #combine treatment success and natural recover terms to get single resistance recovery flows
  n_YaS[] <- n_Ya_tr_success[i] + n_Ya_nr[i]
  n_YbS[] <- n_Yb_tr_success[i] + n_Yb_nr[i]
  n_YcS[] <- n_Yc_tr_success[i] + n_Yc_nr[i]
  n_ZaS[] <- n_Za_tr_success[i] + n_Za_nr[i]
  n_ZbS[] <- n_Zb_tr_success[i] + n_Zb_nr[i]
  n_ZcS[] <- n_Zc_tr_success[i] + n_Zc_nr[i]
  
  #now split the treatment failures into the two different resistance compartment flows
  p_YZaYZab <- E_b*w_b / (E_b*w_b + E_c*w_c + 1e-8) #add very small term to denominator so this probability is always defined
  p_YZbYZab <- E_a*w_a / (E_a*w_a + E_c*w_c + 1e-8)
  p_YZcYZac <- E_a*w_a / (E_a*w_a + E_b*w_b + 1e-8)
  
  n_YaYab[] <- rbinom(n_Ya_tr_fail[i], p_YZaYZab)
  n_YaYac[] <- n_Ya_tr_fail[i] - n_YaYab[i]
  n_YbYab[] <- rbinom(n_Yb_tr_fail[i], p_YZbYZab)
  n_YbYbc[] <- n_Yb_tr_fail[i] - n_YbYab[i]
  n_YcYac[] <- rbinom(n_Yc_tr_fail[i], p_YZcYZac)
  n_YcYbc[] <- n_Yc_tr_fail[i] - n_YcYac[i]
  
  n_ZaZab[] <- rbinom(n_Za_tr_fail[i], p_YZaYZab)
  n_ZaZac[] <- n_Za_tr_fail[i] - n_ZaZab[i]
  n_ZbZab[] <- rbinom(n_Zb_tr_fail[i], p_YZbYZab)
  n_ZbZbc[] <- n_Zb_tr_fail[i] - n_ZbZab[i]
  n_ZcZac[] <- rbinom(n_Zc_tr_fail[i], p_YZcYZac)
  n_ZcZbc[] <- n_Zc_tr_fail[i] - n_ZcZac[i]
  
  ##now move to the double resistance compartments
  #first calculate total outflow -- calculate exit probabilities and then draw
  p_YabOut <- 1 - exp(-(E_c*T_s + (1 - E_c)*kappa*T_sr + d))
  p_YacOut <- 1 - exp(-(E_b*T_s + (1 - E_b)*kappa*T_sr + d))
  p_YbcOut <- 1 - exp(-(E_a*T_s + (1 - E_a)*kappa*T_sr + d))
  p_ZabOut <- 1 - exp(-(E_c*T_m + d))
  p_ZacOut <- 1 - exp(-(E_b*T_m + d))
  p_ZbcOut <- 1 - exp(-(E_a*T_m + d))
  
  n_YabOut[] <- rbinom(Y_ab[i], p_YabOut)
  n_YacOut[] <- rbinom(Y_ac[i], p_YacOut)
  n_YbcOut[] <- rbinom(Y_bc[i], p_YbcOut)
  n_ZabOut[] <- rbinom(Z_ab[i], p_ZabOut)
  n_ZacOut[] <- rbinom(Z_ac[i], p_ZacOut)
  n_ZbcOut[] <- rbinom(Z_bc[i], p_ZbcOut)
  
  #because these are already dual resistance, the only ones leaving and not recovering are developing resistance to the same drug, so separate those and the rest recover
  p_Yab_tr_fail <- (E_c*w_c*T_s)/(E_c*T_s + (1 - E_c)*kappa*T_sr + d)
  p_Yac_tr_fail <- (E_b*w_b*T_s)/(E_b*T_s + (1 - E_b)*kappa*T_sr + d)
  p_Ybc_tr_fail <- (E_a*w_a*T_s)/(E_a*T_s + (1 - E_a)*kappa*T_sr + d)
  p_Zab_tr_fail <- (E_c*w_c*T_m)/(E_c*T_m + d)
  p_Zac_tr_fail <- (E_b*w_b*T_m)/(E_b*T_m + d)
  p_Zbc_tr_fail <- (E_a*w_a*T_m)/(E_a*T_m + d)
  
  #split out resistant flows to Yabc vs recoveries to S
  n_YabYabc[] <- rbinom(n_YabOut[i], p_Yab_tr_fail)
  n_YabS[] <- n_YabOut[i] - n_YabYabc[i]
  n_YacYabc[] <- rbinom(n_YacOut[i], p_Yac_tr_fail)
  n_YacS[] <- n_YacOut[i] - n_YacYabc[i]
  n_YbcYabc[] <- rbinom(n_YbcOut[i], p_Ybc_tr_fail)
  n_YbcS[] <- n_YbcOut[i] - n_YbcYabc[i]
  
  n_ZabZabc[] <- rbinom(n_ZabOut[i], p_Zab_tr_fail)
  n_ZabS[] <- n_ZabOut[i] - n_ZabZabc[i]
  n_ZacZabc[] <- rbinom(n_ZacOut[i], p_Zac_tr_fail)
  n_ZacS[] <- n_ZacOut[i] - n_ZacZabc[i]
  n_ZbcZabc[] <- rbinom(n_ZbcOut[i], p_Zbc_tr_fail)
  n_ZbcS[] <- n_ZbcOut[i] - n_ZbcZabc[i]
  
  #recovery from triple resistance states
  p_YabcOut <- 1 - exp(-(kappa*T_sr + d))
  p_ZabcOut <- 1 - exp(-d)
  n_YabcOut[] <- rbinom(Y_abc[i], p_YabcOut)
  n_ZabcOut[] <- rbinom(Z_abc[i], p_ZabcOut)
  
  
  #check whether prevalence to A has reached the threshold
  threshold_met_A = (prevA >= threshold)
  update(lock_threshold_A) = lock_threshold_A + threshold_met_A
  
  #check whether prevalence to B has reached the threshold
  threshold_met_B = (prevB >= threshold)
  update(lock_threshold_B) = lock_threshold_B + threshold_met_B
  
  #check whether prevalence to C has reached the threshold
  threshold_met_C = (prevC >= threshold)
  update(lock_threshold_C) = lock_threshold_C + threshold_met_C
  
  #make indicator variables for which ordering we're in 
  all_off <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_A <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C < 1)
  only_B <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  only_C <- (lock_threshold_A < 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1)
  A_B_on <- (lock_threshold_A >= 1)*(lock_threshold_B >= 1)*(lock_threshold_C < 1)
  update(C_remain_lock) <- C_remain_lock + A_B_on
  A_C_on <- (lock_threshold_A >= 1)*(lock_threshold_B < 1)*(lock_threshold_C >= 1) 
  update(B_remain_lock) <- B_remain_lock + A_C_on
  B_C_on <- (lock_threshold_A < 1)*(lock_threshold_B >= 1)*(lock_threshold_C >= 1) 
  update(A_remain_lock) <- A_remain_lock + B_C_on
  
  E_a <- 1/3*all_off + 0.5*(only_B + only_C) + max((B_C_on), (A_remain_lock >= 1))
  E_b <- 1/3*all_off + 0.5*(only_A + only_C) + max((A_C_on),(B_remain_lock >= 1))
  E_c <- 1/3*all_off + 0.5*(only_A + only_B) + max((A_B_on),(C_remain_lock >= 1))
  
  sum_tx_probs = E_a + E_c + E_b 
  
  #sum array variables for book-keeping / accounting outputs
  sum_n_SOut <- sum(n_SOut)
  sum_n_SYZa <- sum(n_SYZa)
  sum_n_SYZb <- sum(n_SYZb)
  sum_n_SYZc <- sum(n_SYZc)
  sum_n_SYZab <- sum(n_SYZab)
  sum_n_SYZac <- sum(n_SYZac)
  sum_n_SYZbc <- sum(n_SYZbc)
  sum_n_SYZabc <- sum(n_SYZabc)
  sum_n_SYZs <- sum(n_SYZs)
  sum_n_SYs <- sum(n_SYs)
  sum_n_SZs <- sum(n_SZs)
  sum_n_SYa <- sum(n_SYa)
  sum_n_SZa <- sum(n_SZa)
  sum_n_SYb <- sum(n_SYb)
  sum_n_SZb <- sum(n_SZb)
  sum_n_SYc <- sum(n_SYc)
  sum_n_SZc <- sum(n_SZc)
  sum_n_SYab <- sum(n_SYab)
  sum_n_SZab <- sum(n_SZab)
  sum_n_SYac <- sum(n_SYac)
  sum_n_SZac <- sum(n_SZac)
  sum_n_SYbc <- sum(n_SYbc)
  sum_n_SZbc <- sum(n_SZbc)
  sum_n_SYabc <- sum(n_SYabc)
  sum_n_SZabc <- sum(n_SZabc)
  sum_n_YsOut <- sum(n_YsOut)
  sum_n_ZsOut <- sum(n_ZsOut)
  sum_n_Ys_tr <- sum(n_Ys_tr)
  sum_n_Ys_nr <- sum(n_Ys_nr)
  sum_n_Zs_tr <- sum(n_Zs_tr)
  sum_n_Zs_nr <- sum(n_Zs_nr)
  sum_n_Ys_tr_fail <- sum(n_Ys_tr_fail)
  sum_n_Ys_tr_success <- sum(n_Ys_tr_success)
  sum_n_Zs_tr_fail <- sum(n_Zs_tr_fail)
  sum_n_Zs_tr_success <- sum(n_Zs_tr_success)
  sum_n_YsS <- sum(n_YsS)
  sum_n_ZsS <- sum(n_ZsS)
  sum_n_YsYa <- sum(n_YsYa)
  sum_n_YsYb <- sum(n_YsYb)
  sum_n_YsYc <- sum(n_YsYc)
  sum_n_ZsZa <- sum(n_ZsZa)
  sum_n_ZsZb <- sum(n_ZsZb)
  sum_n_ZsZc <- sum(n_ZsZc)
  sum_n_YaOut <- sum(n_YaOut)
  sum_n_YbOut <- sum(n_YbOut)
  sum_n_YcOut <- sum(n_YcOut)
  sum_n_ZaOut <- sum(n_ZaOut)
  sum_n_ZbOut <- sum(n_ZbOut)
  sum_n_ZcOut <- sum(n_ZcOut)
  sum_n_Ya_nr <- sum(n_Ya_nr)
  sum_n_Ya_tr <- sum(n_Ya_tr)
  sum_n_Yb_nr <- sum(n_Yb_nr)
  sum_n_Yb_tr <- sum(n_Yb_tr)
  sum_n_Yc_nr <- sum(n_Yc_nr)
  sum_n_Yc_tr <- sum(n_Yc_tr)
  sum_n_Za_nr <- sum(n_Za_nr)
  sum_n_Za_tr <- sum(n_Za_tr)
  sum_n_Zb_nr <- sum(n_Zb_nr)
  sum_n_Zb_tr <- sum(n_Zb_tr)
  sum_n_Zc_nr <- sum(n_Zc_nr)
  sum_n_Zc_tr <- sum(n_Zc_tr)
  sum_n_Ya_tr_fail <- sum(n_Ya_tr_fail)
  sum_n_Yb_tr_fail <- sum(n_Yb_tr_fail)
  sum_n_Yc_tr_fail <- sum(n_Yc_tr_fail)
  sum_n_Za_tr_fail <- sum(n_Za_tr_fail)
  sum_n_Zb_tr_fail <- sum(n_Zb_tr_fail)
  sum_n_Zc_tr_fail <- sum(n_Zc_tr_fail)
  sum_n_Ya_tr_success <- sum(n_Ya_tr_success)
  sum_n_Yb_tr_success <- sum(n_Yb_tr_success)
  sum_n_Yc_tr_success <- sum(n_Yc_tr_success)
  sum_n_Za_tr_success <- sum(n_Za_tr_success)
  sum_n_Zb_tr_success <- sum(n_Zb_tr_success)
  sum_n_Zc_tr_success <- sum(n_Zc_tr_success)
  sum_n_YaS <- sum(n_YaS)
  sum_n_YbS <- sum(n_YbS)
  sum_n_YcS <- sum(n_YcS)
  sum_n_ZaS <- sum(n_ZaS)
  sum_n_ZbS <- sum(n_ZbS)
  sum_n_ZcS <- sum(n_ZcS)
  sum_n_YaYab <- sum(n_YaYab)
  sum_n_YaYac <- sum(n_YaYac)
  sum_n_YbYab <- sum(n_YbYab)
  sum_n_YbYbc <- sum(n_YbYbc)
  sum_n_YcYac <- sum(n_YcYac)
  sum_n_YcYbc <- sum(n_YcYbc)
  sum_n_ZaZab <- sum(n_ZaZab)
  sum_n_ZaZac <- sum(n_ZaZac)
  sum_n_ZbZab <- sum(n_ZbZab)
  sum_n_ZbZbc <- sum(n_ZbZbc)
  sum_n_ZcZac <- sum(n_ZcZac)
  sum_n_ZcZbc <- sum(n_ZcZbc)
  sum_n_YabOut <- sum(n_YabOut)
  sum_n_YacOut <- sum(n_YacOut)
  sum_n_YbcOut <- sum(n_YbcOut)
  sum_n_ZabOut <- sum(n_ZabOut)
  sum_n_ZacOut <- sum(n_ZacOut)
  sum_n_ZbcOut <- sum(n_ZbcOut)
  sum_n_YabYabc <- sum(n_YabYabc)
  sum_n_YabS <- sum(n_YabS)
  sum_n_YacYabc <- sum(n_YacYabc)
  sum_n_YacS <- sum(n_YacS)
  sum_n_YbcYabc <- sum(n_YbcYabc)
  sum_n_YbcS<- sum(n_YbcS)
  sum_n_ZabZabc <- sum(n_ZabZabc)
  sum_n_ZabS <- sum(n_ZabS)
  sum_n_ZacZabc <- sum(n_ZacZabc)
  sum_n_ZacS <- sum(n_ZacS)
  sum_n_ZbcZabc <- sum(n_ZbcZabc)
  sum_n_ZbcS<- sum(n_ZbcS)
  sum_n_YabcOut <- sum(n_YabcOut)
  sum_n_ZabcOut <- sum(n_ZabcOut)
  #incidence tracking
  Inc_cat[] = sum(BI_S[,i])*S[i]
  Inc <- sum(Inc_cat)
  ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
  ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
  ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
  ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
  ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
  ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
  
  #set initial conditions
  initial(S[]) <- S_ini[i]
  initial(Y_s[]) <- Y_s_ini[i]
  initial(Z_s[]) <- Z_s_ini[i]
  initial(Y_a[]) <- Y_a_ini[i]
  initial(Z_a[]) <- Z_a_ini[i]
  initial(Y_b[]) <- Y_b_ini[i]
  initial(Z_b[]) <- Z_b_ini[i]
  initial(Y_c[]) <- Y_c_ini[i]
  initial(Z_c[]) <- Z_c_ini[i]
  initial(Y_ab[]) <- Y_ab_ini[i]
  initial(Z_ab[]) <- Z_ab_ini[i]
  initial(Y_ac[]) <- Y_ac_ini[i]
  initial(Z_ac[]) <- Z_ac_ini[i]  
  initial(Y_bc[]) <- Y_bc_ini[i]
  initial(Z_bc[]) <- Z_bc_ini[i]   
  initial(Y_abc[]) <- Y_abc_ini[i]
  initial(Z_abc[]) <- Z_abc_ini[i]
  initial(lock_threshold_A) <- 0
  initial(lock_threshold_B) <- 0
  initial(lock_threshold_C) <- 0
  initial(A_remain_lock) <- 0
  initial(B_remain_lock) <- 0
  initial(C_remain_lock) <- 0
  
  #user-defined states
  S_ini[] <- user()
  Y_s_ini[] <- user()
  Z_s_ini[] <- user()
  Y_a_ini[] <- user()
  Z_a_ini[] <- user()
  Y_b_ini[] <- user()
  Z_b_ini[] <- user()
  Y_c_ini[] <- user()
  Z_c_ini[] <- user()
  Y_ab_ini[] <- user()
  Z_ab_ini[] <- user()
  Y_ac_ini[] <- user()
  Z_ac_ini[] <- user()
  Y_bc_ini[] <- user()
  Z_bc_ini[] <- user() 
  Y_abc_ini[] <- user()
  Z_abc_ini[] <- user()
  beta[,] <- user()
  
  #user-defined parameters
  threshold <- user()
  w_a <- user()
  w_b <- user()
  w_c <- user()
  T_s <- user()
  T_m <- user()
  T_sr <- user()
  kappa <- user()
  f_a <- user()
  f_b <- user()
  f_c <- user()
  f_ab <- user()
  f_ac <- user()
  f_bc <- user()
  f_abc <- user()
  sigma <- user()
  N_risk <- user()
  d <- user()
  
  #define variable dimensions
  ##state initial values
  dim(S_ini) <- N_risk
  dim(Y_s_ini) <- N_risk
  dim(Z_s_ini) <- N_risk
  dim(Y_a_ini) <- N_risk
  dim(Z_a_ini) <- N_risk
  dim(Y_b_ini) <- N_risk
  dim(Z_b_ini) <- N_risk
  dim(Y_c_ini) <- N_risk
  dim(Z_c_ini) <- N_risk
  dim(Y_ab_ini) <- N_risk
  dim(Z_ab_ini) <- N_risk
  dim(Y_ac_ini) <- N_risk
  dim(Z_ac_ini) <- N_risk
  dim(Y_bc_ini) <- N_risk
  dim(Z_bc_ini) <- N_risk
  dim(Y_abc_ini) <- N_risk
  dim(Z_abc_ini) <- N_risk
  dim(sum_BI_S) <- N_risk
  dim(sum_BI_YZs) <- N_risk
  dim(sum_BI_YZa) <- N_risk
  dim(sum_BI_YZb) <- N_risk
  dim(sum_BI_YZc) <- N_risk
  dim(sum_BI_YZab) <- N_risk
  dim(sum_BI_YZac) <- N_risk
  dim(sum_BI_YZbc) <- N_risk
  dim(sum_BI_YZabc) <- N_risk
  dim(sum_all_BI_terms) <- N_risk
  dim(ninf_Ys) <- N_risk
  dim(ninf_Zs) <- N_risk
  dim(ninf_Ya) <- N_risk
  dim(ninf_Za) <- N_risk
  dim(ninf_Yb) <- N_risk
  dim(ninf_Zb) <- N_risk
  dim(ninf_Yab) <- N_risk
  dim(ninf_Zab) <- N_risk
  dim(n_SOut) <- N_risk
  dim(p_SOut) <- N_risk
  dim(p_SYZa) <- N_risk
  dim(p_SYZb) <- N_risk
  dim(p_SYZc) <- N_risk
  dim(p_SYZab) <- N_risk
  dim(p_SYZac) <- N_risk  
  dim(p_SYZbc) <- N_risk
  dim(p_SYZabc) <- N_risk
  dim(n_SYZa) <- N_risk
  dim(n_SYZb) <- N_risk
  dim(n_SYZc) <- N_risk
  dim(n_SYZab) <- N_risk
  dim(n_SYZac) <- N_risk
  dim(n_SYZbc) <- N_risk
  dim(n_SYZabc) <- N_risk
  dim(n_SYZs) <- N_risk
  dim(n_SYs) <- N_risk
  dim(n_SZs) <- N_risk
  dim(n_SYa) <- N_risk
  dim(n_SZa) <- N_risk
  dim(n_SYb) <- N_risk
  dim(n_SZb) <- N_risk
  dim(n_SYc) <- N_risk
  dim(n_SZc) <- N_risk
  dim(n_SYab) <- N_risk
  dim(n_SZab) <- N_risk
  dim(n_SYac) <- N_risk
  dim(n_SZac) <- N_risk
  dim(n_SYbc) <- N_risk
  dim(n_SZbc) <- N_risk
  dim(n_SYabc) <- N_risk
  dim(n_SZabc) <- N_risk
  dim(n_YsOut) <- N_risk
  dim(n_ZsOut) <- N_risk
  dim(n_Ys_tr) <- N_risk
  dim(n_Ys_nr) <- N_risk
  dim(n_Zs_tr) <- N_risk
  dim(n_Zs_nr) <- N_risk
  dim(n_Ys_tr_fail) <- N_risk
  dim(n_Ys_tr_success) <- N_risk
  dim(n_Zs_tr_fail) <- N_risk
  dim(n_Zs_tr_success) <- N_risk
  dim(n_YsS) <- N_risk
  dim(n_ZsS) <- N_risk
  dim(n_YsYa) <- N_risk
  dim(n_YsYb) <- N_risk
  dim(n_YsYc) <- N_risk
  dim(n_ZsZa) <- N_risk
  dim(n_ZsZb) <- N_risk
  dim(n_ZsZc) <- N_risk
  dim(n_YaOut) <- N_risk
  dim(n_YbOut) <- N_risk
  dim(n_YcOut) <- N_risk
  dim(n_ZaOut) <- N_risk
  dim(n_ZbOut) <- N_risk
  dim(n_ZcOut) <- N_risk
  dim(n_Ya_nr) <- N_risk
  dim(n_Ya_tr) <- N_risk
  dim(n_Yb_nr) <- N_risk
  dim(n_Yb_tr) <- N_risk
  dim(n_Yc_nr) <- N_risk
  dim(n_Yc_tr) <- N_risk
  dim(n_Za_nr) <- N_risk
  dim(n_Za_tr) <- N_risk
  dim(n_Zb_nr) <- N_risk
  dim(n_Zb_tr) <- N_risk
  dim(n_Zc_nr) <- N_risk
  dim(n_Zc_tr) <- N_risk
  dim(n_Ya_tr_fail) <- N_risk
  dim(n_Yb_tr_fail) <- N_risk
  dim(n_Yc_tr_fail) <- N_risk
  dim(n_Za_tr_fail) <- N_risk
  dim(n_Zb_tr_fail) <- N_risk
  dim(n_Zc_tr_fail) <- N_risk
  dim(n_Ya_tr_success) <- N_risk
  dim(n_Yb_tr_success) <- N_risk
  dim(n_Yc_tr_success) <- N_risk
  dim(n_Za_tr_success) <- N_risk
  dim(n_Zb_tr_success) <- N_risk
  dim(n_Zc_tr_success) <- N_risk
  dim(n_YaS) <- N_risk
  dim(n_YbS) <- N_risk
  dim(n_YcS) <- N_risk
  dim(n_ZaS) <- N_risk
  dim(n_ZbS) <- N_risk
  dim(n_ZcS) <- N_risk
  dim(n_YaYab) <- N_risk
  dim(n_YaYac) <- N_risk
  dim(n_YbYab) <- N_risk
  dim(n_YbYbc) <- N_risk
  dim(n_YcYac) <- N_risk
  dim(n_YcYbc) <- N_risk
  dim(n_ZaZab) <- N_risk
  dim(n_ZaZac) <- N_risk
  dim(n_ZbZab) <- N_risk
  dim(n_ZbZbc) <- N_risk
  dim(n_ZcZac) <- N_risk
  dim(n_ZcZbc) <- N_risk
  dim(n_YabOut) <- N_risk
  dim(n_YacOut) <- N_risk
  dim(n_YbcOut) <- N_risk
  dim(n_ZabOut) <- N_risk
  dim(n_ZacOut) <- N_risk
  dim(n_ZbcOut) <- N_risk
  dim(n_YabYabc) <- N_risk
  dim(n_YabS) <- N_risk
  dim(n_YacYabc) <- N_risk
  dim(n_YacS) <- N_risk
  dim(n_YbcYabc) <- N_risk
  dim(n_YbcS) <- N_risk
  dim(n_ZabZabc) <- N_risk
  dim(n_ZabS) <- N_risk
  dim(n_ZacZabc) <- N_risk
  dim(n_ZacS) <- N_risk
  dim(n_ZbcZabc) <- N_risk
  dim(n_ZbcS) <- N_risk
  dim(n_YabcOut) <- N_risk
  dim(n_ZabcOut) <- N_risk
  
  ##states
  dim(S) <- N_risk
  dim(Y_s) <- N_risk
  dim(Z_s) <- N_risk
  dim(Y_a) <- N_risk
  dim(Z_a) <- N_risk
  dim(Y_b) <- N_risk
  dim(Z_b) <- N_risk
  dim(Y_c) <- N_risk
  dim(Z_c) <- N_risk
  dim(Y_ab) <- N_risk
  dim(Z_ab) <- N_risk
  dim(Y_ac) <- N_risk
  dim(Z_ac) <- N_risk
  dim(Y_bc) <- N_risk
  dim(Z_bc) <- N_risk
  dim(Y_abc) <- N_risk
  dim(Z_abc) <- N_risk
  ##other vector / matrix variables
  dim(beta) <- c(N_risk, N_risk)
  dim(N) <- N_risk
  dim(BI_S) <- c(N_risk, N_risk)
  dim(BI_YZs) <- c(N_risk, N_risk)
  dim(BI_YZa) <- c(N_risk, N_risk)
  dim(BI_YZb) <- c(N_risk, N_risk)
  dim(BI_YZc) <- c(N_risk, N_risk)
  dim(BI_YZab) <- c(N_risk, N_risk)
  dim(BI_YZac) <- c(N_risk, N_risk)
  dim(BI_YZbc) <- c(N_risk, N_risk)
  dim(BI_YZabc) <- c(N_risk, N_risk)
  dim(Inc_cat) <- N_risk
  
  
  #define desired non-state outputs
  output(N) <- N
  output(prevA) <- prevA
  output(prevB) <- prevB
  output(prevC) <- prevC
  output(Inc_cat) <- Inc_cat
  output(Inc) <- Inc
  output(E_a) <- E_a
  output(E_b) <- E_b
  output(E_c) <- E_c
  output(all_a_res) <- all_a_res
  output(all_b_res) <- all_b_res
  output(all_c_res) <- all_c_res
  output(all_non_res) <- all_non_res
  output(sum_S) <- sum_S
  output(sum_Y_s) <- sum_Y_s
  output(sum_Y_a) <- sum_Y_a
  output(sum_Y_b) <- sum_Y_b
  output(sum_Y_c) <- sum_Y_c
  output(sum_Y_ab) <- sum_Y_ab
  output(sum_Y_ac) <- sum_Y_ac
  output(sum_Y_bc) <- sum_Y_bc
  output(sum_Y_abc) <- sum_Y_abc
  output(sum_Z_s) <- sum_Z_s
  output(sum_Z_a) <- sum_Z_a
  output(sum_Z_b) <- sum_Z_b
  output(sum_Z_c) <- sum_Z_c
  output(sum_Z_ab) <- sum_Z_ab
  output(sum_Z_ac) <- sum_Z_ac
  output(sum_Z_bc) <- sum_Z_bc
  output(sum_Z_abc) <- sum_Z_abc
  output(sum_YZ_s) <- sum_YZ_s
  output(sum_YZ_a) <- sum_YZ_a
  output(sum_YZ_b) <- sum_YZ_b
  output(sum_YZ_c) <- sum_YZ_c
  output(sum_YZ_ab) <- sum_YZ_ab
  output(sum_YZ_ac) <- sum_YZ_ac
  output(sum_YZ_bc) <- sum_YZ_bc
  output(sum_YZ_abc) <- sum_YZ_abc
  output(sum_N) <- sum_N
  output(prev) <- prev
})














#Update model 6/5 to fix errors and make stick breaking less error prone
# sequential_stochastic_6_5_fail <- odin({
#   ##write out all intermediate transmission terms needed for equations
#   #write out intermediate equation for beta*I for the dS equation
#   BI_S[,] <- beta[i,j]*((Y_s[i]+Z_s[i]) +
#                           f_a*(Y_a[i] + Z_a[i]) +
#                           f_b*(Y_b[i] + Z_b[i]) +
#                           f_c*(Y_c[i] + Z_c[i]) + 
#                           f_ab*(Y_ab[i] + Z_ab[i]) + 
#                           f_ac*(Y_ac[i] + Z_ac[i]) + 
#                           f_bc*(Y_bc[i] + Z_bc[i]) + 
#                           f_abc*(Y_abc[i] + Z_abc[i]))
#   
#   #write out intermediate equation for beta*I for the dY_s and dZ_s equations
#   BI_YZs[,] <- beta[i,j]*(Y_s[i] + Z_s[i])
#   
#   #write out intermediate equation for beta*I for the dY_a and dZ_a equations
#   BI_YZa[,] <- f_a*beta[i,j]*(Y_a[i] + Z_a[i])
#   
#   #write out intermediate equation for beta*I for the dY_b and dZ_b equations
#   BI_YZb[,] <- f_b*beta[i,j]*(Y_b[i] + Z_b[i])
#   
#   #write out intermediate equation for beta*I for the dY_c and dZ_c equations
#   BI_YZc[,] <- f_c*beta[i,j]*(Y_c[i] + Z_c[i])
#   
#   #write out intermediate equation for beta*I for dY_ab and dZ_ab
#   BI_YZab[,] <- f_ab*beta[i,j]*(Y_ab[i] + Z_ab[i])
#   
#   #write out intermediate equation for beta*I for dY_ac and dZ_ac
#   BI_YZac[,] <- f_ac*beta[i,j]*(Y_ac[i] + Z_ac[i])
#   
#   #write out intermediate equation for beta*I for dY_bc and dZ_bc
#   BI_YZbc[,] <- f_bc*beta[i,j]*(Y_bc[i] + Z_bc[i])
#   
#   #write out intermediate equation for dY_abc and dZ_abc
#   BI_YZabc[,] <- f_abc*beta[i,j]*(Y_abc[i] + Z_abc[i])
#   
#   ##summarize the intermediate terms as a means of summarizing for bookkeeping purposes 
#   sum_BI_S[] <- sum(BI_S[,i])
#   sum_BI_YZs[] <- sum(BI_YZs[,i])
#   sum_BI_YZa[] <- sum(BI_YZa[,i])
#   sum_BI_YZb[] <- sum(BI_YZb[,i])
#   sum_BI_YZc[] <- sum(BI_YZc[,i])
#   sum_BI_YZab[] <- sum(BI_YZab[,i])
#   sum_BI_YZac[] <- sum(BI_YZac[,i])
#   sum_BI_YZbc[] <- sum(BI_YZbc[,i])
#   sum_BI_YZabc[] <- sum(BI_YZabc[,i])
#   #sum all BI terms to get the total efflux from S
#   sum_all_BI_terms[] <- sum_BI_YZs[i] + sum_BI_YZa[i] + sum_BI_YZb[i] + sum_BI_YZc[i] + 
#     sum_BI_YZab[i] + sum_BI_YZac[i] + sum_BI_YZbc[i] + sum_BI_YZabc[i]
#   
#   ##total across risk strata of BI terms for bookkeeping (to output)
#   total_sum_BI_S <- sum(sum_BI_S)
#   total_sum_BI_YZs <- sum(sum_BI_YZs)
#   total_sum_BI_YZa <- sum(sum_BI_YZa)
#   total_sum_BI_YZb <- sum(sum_BI_YZb)
#   total_sum_BI_YZc <- sum(sum_BI_YZc)
#   total_sum_BI_YZab <- sum(sum_BI_YZab)
#   total_sum_BI_YZac <- sum(sum_BI_YZac)
#   total_sum_BI_YZbc <- sum(sum_BI_YZbc)
#   total_sum_BI_YZabc <- sum(sum_BI_YZabc)
#   total_sum_all_BI_terms <- sum(sum_all_BI_terms)
#   
#   #write out equation for S
#   update(S[]) <- S[i] - n_SOut[i] + n_YsS[i] + n_YaS[i] + n_YbS[i] + n_YcS[i] + n_YabS[i] + n_YacS[i] + n_YbcS[i] +
#     n_YabcOut[i] + n_ZsS[i] + n_ZaS[i] + n_ZbS[i] + n_ZcS[i] + n_ZabS[i] + n_ZacS[i] + n_ZbcS[i] + n_ZabcOut[i] 
#   
#   #write out equations for dY_s and dZ_s
#   update(Y_s[]) <- Y_s[i] + n_SYs[i] - n_YsOut[i]
#   update(Z_s[]) <- Z_s[i] + n_SZs[i] - n_ZsOut[i]
#   
#   #write out equations for dY_a and dZ_a
#   update(Y_a[]) <- Y_a[i] + n_SYa[i] + n_YsYa[i] - n_YaOut[i]
#   update(Z_a[]) <- Z_a[i] + n_SZa[i] + n_ZsZa[i] - n_ZaOut[i]
#   
#   #write out equations for dY_b and dZ_b
#   update(Y_b[]) <- Y_b[i] + n_SYb[i] + n_YsYb[i] - n_YbOut[i]
#   update(Z_b[]) <- Z_b[i] + n_SZb[i] + n_ZsZb[i] - n_ZbOut[i]
#   
#   #write out equations for dY_c and dZ_c
#   update(Y_c[]) <- Y_c[i] + n_SYc[i] + n_YsYc[i] - n_YcOut[i] 
#   update(Z_c[]) <- Z_c[i] + n_SZc[i] + n_ZsZc[i] - n_ZcOut[i]
#   
#   #write out equations for dY_ab and dZ_ab
#   update(Y_ab[]) <- Y_ab[i] + n_SYab[i] + n_YaYab[i] + n_YbYab[i] - n_YabOut[i]
#   update(Z_ab[]) <- Z_ab[i] + n_SZab[i] + n_ZaZab[i] + n_ZbZab[i] - n_ZabOut[i]
#   
#   #write out equations for dY_ac and dZ_ac
#   update(Y_ac[]) <- Y_ac[i] + n_SYac[i] + n_YaYac[i] + n_YcYac[i] - n_YacOut[i]
#   update(Z_ac[]) <- Z_ac[i] + n_SZac[i] + n_ZaZac[i] + n_ZcZac[i] - n_ZacOut[i] 
#   
#   #write out equations for dY_ac and dZ_ac
#   update(Y_bc[]) <- Y_bc[i] + n_SYbc[i] + n_YbYbc[i] + n_YcYbc[i] - n_YbcOut[i]
#   update(Z_bc[]) <- Z_bc[i] + n_SZbc[i] + n_ZbZbc[i] + n_ZcZbc[i] - n_ZbcOut[i]
#   
#   #write out equations for dY_abc and dZ_abc
#   update(Y_abc[]) <- Y_abc[i] + n_SYabc[i] + n_YabYabc[i] + n_YacYabc[i] + n_YbcYabc[i] - n_YabcOut[i]
#   update(Z_abc[]) <- Z_abc[i] + n_SZabc[i] + n_ZabZabc[i] + n_ZacZabc[i] + n_ZbcZabc[i] - n_ZabcOut[i]
#   
#   #equations for non-state outputs
#   N[] <- S[i] + Y_s[i] + Z_s[i] + Y_a[i] + Z_a[i] + Y_b[i] + Z_b[i] + Y_c[i] + Z_c[i] + Y_ab[i] + Z_ab[i] + Y_ac[i] + Z_ac[i] + Y_bc[i] + 
#     Z_bc[i] + Y_abc[i] + Z_abc[i]
#   
#   sum_N <- sum(N)
#   
#   #prevalence calculatoins: any compartment with resistance to that drug (not just single resistance)
#   prevA <- (sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
#   prevB <- (sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
#   prevC <- (sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc))/(sum(N)-sum(S))
#   
#   prev <- (sum(Y_s) + sum(Z_s) + sum(Y_a) + sum(Z_a) + sum(Y_b) + sum(Z_b) + sum(Y_c) + sum(Z_c) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) +
#              sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)) / sum(N)
#   
#   all_a_res <- sum(Y_a) + sum(Z_a) + sum(Y_ab) + sum(Z_ab) + sum(Y_ac) + sum(Z_ac) + sum(Y_abc) + sum(Z_abc)
#   all_b_res <- sum(Y_b) + sum(Z_b) + sum(Y_ab) + sum(Z_ab) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
#   all_c_res <- sum(Y_c) + sum(Z_c) + sum(Y_ac) + sum(Z_ac) + sum(Y_bc) + sum(Z_bc) + sum(Y_abc) + sum(Z_abc)
#   all_non_res <- sum(Y_s) + sum(Z_s)
#   
#   #sum states at any point in time across risk groups for book-keeping
#   sum_S <- sum(S)
#   sum_Y_s <- sum(Y_s)
#   sum_Y_a <- sum(Y_a)
#   sum_Y_b <- sum(Y_b)
#   sum_Y_c <- sum(Y_c)
#   sum_Y_ab <- sum(Y_ab)
#   sum_Y_ac <- sum(Y_ac)
#   sum_Y_bc <- sum(Y_bc)
#   sum_Y_abc <- sum(Y_abc)
#   sum_Z_s <- sum(Z_s)
#   sum_Z_a <- sum(Z_a)
#   sum_Z_b <- sum(Z_b)
#   sum_Z_c <- sum(Z_c)
#   sum_Z_ab <- sum(Z_ab)
#   sum_Z_ac <- sum(Z_ac)
#   sum_Z_bc <- sum(Z_bc)
#   sum_Z_abc <- sum(Z_abc)
#   
#   #sum states across symptom status for book-keeping
#   sum_YZ_s <- sum_Y_s + sum_Z_s
#   sum_YZ_a <- sum_Y_a + sum_Z_a
#   sum_YZ_b <- sum_Y_b + sum_Z_b
#   sum_YZ_c <- sum_Y_c + sum_Z_c 
#   sum_YZ_ab <- sum_Y_ab + sum_Z_ab
#   sum_YZ_ac <- sum_Y_ac + sum_Z_ac
#   sum_YZ_bc <- sum_Y_bc + sum_Z_bc
#   sum_YZ_abc <- sum_Y_abc + sum_Z_abc
#   
#   ##STOCHASTIC calculations
#   #calculate the probability of leaving S from the rate
#   p_SOut[] <- 1 - exp(-sum_all_BI_terms[i])
#   
#   #calculate the number of people leaving S
#   n_SOut[] <- rbinom(S[i], p_SOut[i])
#   
#   #divide probabilities leaving S 
#   p_SYZa[] <- sum_BI_YZa[i]/sum_all_BI_terms[i]
#   p_SYZb[] <- sum_BI_YZb[i]/sum_all_BI_terms[i]
#   p_SYZc[] <- sum_BI_YZc[i]/sum_all_BI_terms[i]
#   p_SYZab[] <- sum_BI_YZab[i]/sum_all_BI_terms[i]
#   p_SYZac[] <- sum_BI_YZac[i]/sum_all_BI_terms[i]
#   p_SYZbc[] <- sum_BI_YZbc[i]/sum_all_BI_terms[i]
#   p_SYZabc[] <- sum_BI_YZabc[i]/sum_all_BI_terms[i] #fixed error here 6/3
#   
#   #establish a variable to keep track of n_SOut remaining and probability remaining
#   n_SOut_remaining[] <- n_SOut[i]
#   p_SOut_remaining[] <- 1 - p_SYZa[i]
#   
#   #divide the people leaving S into different flows
#   n_SYZa[] <- rbinom(n_SOut_remaining[i], p_SYZa[i])
#   #take out what has already been removed from what is remaining
#   n_SOut_remaining[] <- n_SOut[i] - n_SYZa[i]
#   n_SYZb[] <- rbinom(n_SOut_remaining[i], (p_SYZb[i])/(p_SOut_remaining[i] + 1e-8)) #i) subtract those already subtracted from n_SOut ii) renormalize the probability of exit
#   n_SOut_remaining2[] <- n_SOut_remaining[i] - n_SYZb[i]
#   p_SOut_remaining2[] <- p_SOut_remaining[i] - p_SYZb[i]
#   n_SYZc[] <- rbinom(n_SOut_remaining2[i], (p_SYZc[i])/(p_SOut_remaining2[i] + 1e-8))
#   n_SOut_remaining3[] <- n_SOut_remaining2[i] - n_SYZc[i]
#   p_SOut_remaining3[] <- p_SOut_remaining2[i] - p_SYZc[i]
#   n_SYZab[] <- rbinom(n_SOut_remaining3[i], (p_SYZab[i])/(p_SOut_remaining3[i] + 1e-8))
#   n_SOut_remaining4[] <- n_SOut_remaining3[i] - n_SYZab[i]
#   p_SOut_remaining4[] <- p_SOut_remaining3[i] - p_SYZab[i]
#   n_SYZac[] <- rbinom(n_SOut_remaining4[i], (p_SYZac[i])/(p_SOut_remaining4[i] + 1e-8)) #was an issue here!
#   n_SOut_remaining5[] <- n_SOut_remaining4[i] - n_SYZac[i]
#   p_SOut_remaining5[] <- p_SOut_remaining4[i] - p_SYZac[i]
#   n_SYZbc[] <- rbinom(n_SOut_remaining5[i], (p_SYZbc[i])/(p_SOut_remaining5[i] + 1e-8)) #FIX THIS 6/5
#   n_SOut_remaining6[] <- n_SOut_remaining5[i] - n_SYZbc[i]
#   p_SOut_remaining6[] <- p_SOut_remaining5[i] - p_SYZbc[i]
#   n_SYZabc[] <- rbinom(n_SOut_remaining6[i], (p_SYZabc[i])/(p_SOut_remaining6[i] + 1e-8))
#   n_SYZs[] <- n_SOut_remaining6[i] - n_SYZabc[i]
#   
#   #split into symptomatic vs asymptomatics
#   n_SYs[] <- rbinom(n_SYZs[i], sigma)
#   n_SZs[] <- n_SYZs[i] - n_SYs[i] #fixed issue here
#   n_SYa[] <- rbinom(n_SYZa[i], sigma)
#   n_SZa[] <- n_SYZa[i] - n_SYa[i]
#   n_SYb[] <- rbinom(n_SYZb[i], sigma)
#   n_SZb[] <- n_SYZb[i] - n_SYb[i]
#   n_SYc[] <- rbinom(n_SYZc[i], sigma)
#   n_SZc[] <- n_SYZc[i] - n_SYc[i]
#   n_SYab[] <- rbinom(n_SYZab[i], sigma)
#   n_SZab[] <- n_SYZab[i] - n_SYab[i]
#   n_SYac[] <- rbinom(n_SYZac[i], sigma)
#   n_SZac[] <- n_SYZac[i] - n_SYac[i]
#   n_SYbc[] <- rbinom(n_SYZbc[i], sigma)
#   n_SZbc[] <- n_SYZbc[i] - n_SYbc[i]
#   n_SYabc[] <- rbinom(n_SYZabc[i], sigma)
#   n_SZabc[] <- n_SYZabc[i] - n_SYabc[i]
#   
#   #calculate the probability of leaving Ys and Zs by converting exit rates
#   p_YsOut <- 1 - exp(-(T_s + d))
#   p_ZsOut <- 1 - exp(-(T_m + d))
#   
#   #calculate the number of people leaving Ys and Zs based on the probabilities
#   n_YsOut[] <- rbinom(Y_s[i], p_YsOut)
#   n_ZsOut[] <- rbinom(Z_s[i], p_ZsOut)
#   
#   #split the people recovering naturally (nr) vs the people leaving due to treatment (tr)
#   p_Ys_tr <- T_s / (T_s + d)
#   n_Ys_tr[] <- rbinom(n_YsOut[i], p_Ys_tr)
#   n_Ys_nr[] <- n_YsOut[i] - n_Ys_tr[i]
#   p_Zs_tr <- T_m / (T_m + d)
#   n_Zs_tr[] <- rbinom(n_ZsOut[i], p_Zs_tr)
#   n_Zs_nr[] <- n_ZsOut[i] - n_Zs_tr[i]
#   
#   #of the people being treated, split the flows into recovery vs. resistance
#   p_YZs_tr_fail <- E_a*w_a + E_b*w_b + E_c*w_c
#   n_Ys_tr_fail[] <- rbinom(n_Ys_tr[i], p_YZs_tr_fail)
#   n_Ys_tr_success[] <- n_Ys_tr[i] - n_Ys_tr_fail[i]
#   n_Zs_tr_fail[] <- rbinom(n_Zs_tr[i], p_YZs_tr_fail)
#   n_Zs_tr_success[] <- n_Zs_tr[i] - n_Zs_tr_fail[i]
#   
#   #combine those recovering from successful treatment and natural recovery as leaving Ys and Zs and going back to S
#   n_YsS[] <- n_Ys_tr_success[i] + n_Ys_nr[i]
#   n_ZsS[] <- n_Zs_tr_success[i] + n_Zs_nr[i]
#   
#   #split those who are treated and develop resistance into the three possible flows
#   #first take the individually normalized probabilities of each flow
#   p_YZs_YZa <- E_a*w_a/p_YZs_tr_fail
#   p_YZs_YZb <- E_b*w_b/p_YZs_tr_fail
#   p_YZs_YZc <- E_c*w_c/p_YZs_tr_fail
#   
#   #do stick-breaking with Y
#   n_YsYa[] <- rbinom(n_Ys_tr_fail[i], p_YZs_YZa)
#   n_YsYb[] <- rbinom(n_Ys_tr_fail[i] - n_YsYa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
#   n_YsYc[] <- n_Ys_tr_fail[i] - n_YsYa[i] - n_YsYb[i]
#   
#   #do stick-breaking with Z
#   n_ZsZa[] <- rbinom(n_Zs_tr_fail[i], p_YZs_YZa)
#   n_ZsZb[] <- rbinom(n_Zs_tr_fail[i] - n_ZsZa[i], (p_YZs_YZb /(1 - p_YZs_YZa + 1e-8)))#only want to draw from people who didn't go to Ya, and renormalize probability
#   n_ZsZc[] <- n_Zs_tr_fail[i] - n_ZsZa[i] - n_ZsZb[i]
#   
#   ##now onto single resistance compartments. Probabilities of exit will be different for Y and Z compartments due to retreatment
#   #probabilities for exit flows from Ya, Yb, Yc, Za, Zb, and Zc
#   p_YaOut <- 1 - exp(-((E_b + E_c)*T_s + E_a*kappa*T_sr + d))
#   p_YbOut <- 1 - exp(-((E_a + E_c)*T_s + E_b*kappa*T_sr + d))
#   p_YcOut <- 1 - exp(-((E_a + E_b)*T_s + E_c*kappa*T_sr + d))
#   p_ZaOut <- 1 - exp(-((E_b + E_c)*T_m + d))
#   p_ZbOut <- 1 - exp(-((E_a + E_c)*T_m + d))
#   p_ZcOut <- 1 - exp(-((E_a + E_b)*T_m + d))
#   
#   #calculate those leaving the single resistance compartments
#   n_YaOut[] <- rbinom(Y_a[i], p_YaOut)
#   n_YbOut[] <- rbinom(Y_b[i], p_YbOut)
#   n_YcOut[] <- rbinom(Y_c[i], p_YcOut)
#   n_ZaOut[] <- rbinom(Z_a[i], p_ZaOut)
#   n_ZbOut[] <- rbinom(Z_b[i], p_ZbOut)
#   n_ZcOut[] <- rbinom(Z_c[i], p_ZcOut)
#   
#   #calculate the probabilities of leaving single resistance compartments due to natural recovery
#   p_Ya_nr <- d/((E_b + E_c)*T_s + E_a*kappa*T_sr + d)
#   p_Yb_nr <- d/((E_a + E_c)*T_s + E_b*kappa*T_sr + d)
#   p_Yc_nr <- d/((E_a + E_b)*T_s + E_c*kappa*T_sr + d)
#   p_Za_nr <- d/((E_b + E_c)*T_m + d)
#   p_Zb_nr <- d/((E_a + E_c)*T_m + d)
#   p_Zc_nr <- d/((E_a + E_b)*T_m + d)
#   
#   #split the people leaving single resistance into natural recovery vs treatment
#   n_Ya_nr[] <- rbinom(n_YaOut[i], p_Ya_nr)
#   n_Ya_tr[] <- n_YaOut[i] - n_Ya_nr[i]
#   n_Yb_nr[] <- rbinom(n_YbOut[i], p_Yb_nr)
#   n_Yb_tr[] <- n_YbOut[i] - n_Yb_nr[i]
#   n_Yc_nr[] <- rbinom(n_YcOut[i], p_Yc_nr)
#   n_Yc_tr[] <- n_YcOut[i] - n_Yc_nr[i]
#   n_Za_nr[] <- rbinom(n_ZaOut[i], p_Za_nr)
#   n_Za_tr[] <- n_ZaOut[i] - n_Za_nr[i]
#   n_Zb_nr[] <- rbinom(n_ZbOut[i], p_Zb_nr)
#   n_Zb_tr[] <- n_ZbOut[i] - n_Zb_nr[i]
#   n_Zc_nr[] <- rbinom(n_ZcOut[i], p_Zc_nr)
#   n_Zc_tr[] <- n_ZcOut[i] - n_Zc_nr[i]
#   
#   #split those leaving single resistance after treatment into those who fail treatment vs. those who recover
#   #note that here, treatment failure means people who develop resistance on treatment, not those treated with the drug they are resistance bc they stay in the resistant compartment
#   p_Ya_tr_fail <- (E_b*w_b*T_s + E_c*w_c*T_s)/((E_b + E_c)*T_s + E_a*kappa*T_sr)
#   p_Yb_tr_fail <- (E_a*w_a*T_s + E_c*w_c*T_s)/((E_a + E_c)*T_s + E_b*kappa*T_sr)
#   p_Yc_tr_fail <- (E_a*w_a*T_s + E_b*w_b*T_s)/((E_a + E_b)*T_s + E_c*kappa*T_sr)
#   p_Za_tr_fail <- (E_b*w_b + E_c*w_c)/ ((E_b + E_c) + 1e-8) #add very small term to denominator so this probability is always defined
#   p_Zb_tr_fail <- (E_a*w_a + E_c*w_c)/ ((E_a + E_c) + 1e-8)
#   p_Zc_tr_fail <- (E_a*w_a + E_b*w_b)/ ((E_a + E_b) + 1e-8)
#   
#   n_Ya_tr_fail[] <- rbinom(n_Ya_tr[i], p_Ya_tr_fail)
#   n_Yb_tr_fail[] <- rbinom(n_Yb_tr[i], p_Yb_tr_fail)
#   n_Yc_tr_fail[] <- rbinom(n_Yc_tr[i], p_Yc_tr_fail)
#   n_Za_tr_fail[] <- rbinom(n_Za_tr[i], p_Za_tr_fail)
#   n_Zb_tr_fail[] <- rbinom(n_Zb_tr[i], p_Zb_tr_fail)
#   n_Zc_tr_fail[] <- rbinom(n_Zc_tr[i], p_Zc_tr_fail)
#   
#   n_Ya_tr_success[] <- n_Ya_tr[i] - n_Ya_tr_fail[i]
#   n_Yb_tr_success[] <- n_Yb_tr[i] - n_Yb_tr_fail[i]
#   n_Yc_tr_success[] <- n_Yc_tr[i] - n_Yc_tr_fail[i]
#   n_Za_tr_success[] <- n_Za_tr[i] - n_Za_tr_fail[i]
#   n_Zb_tr_success[] <- n_Zb_tr[i] - n_Zb_tr_fail[i]
#   n_Zc_tr_success[] <- n_Zc_tr[i] - n_Zc_tr_fail[i]
#   
#   #combine treatment success and natural recover terms to get single resistance recovery flows
#   n_YaS[] <- n_Ya_tr_success[i] + n_Ya_nr[i]
#   n_YbS[] <- n_Yb_tr_success[i] + n_Yb_nr[i]
#   n_YcS[] <- n_Yc_tr_success[i] + n_Yc_nr[i]
#   n_ZaS[] <- n_Za_tr_success[i] + n_Za_nr[i]
#   n_ZbS[] <- n_Zb_tr_success[i] + n_Zb_nr[i]
#   n_ZcS[] <- n_Zc_tr_success[i] + n_Zc_nr[i]
#   
#   #now split the treatment failures into the two different resistance compartment flows
#   p_YZaYZab <- E_b*w_b / (E_b*w_b + E_c*w_c + 1e-8) #add very small term to denominator so this probability is always defined
#   p_YZbYZab <- E_a*w_a / (E_a*w_a + E_c*w_c + 1e-8)
#   p_YZcYZac <- E_a*w_a / (E_a*w_a + E_b*w_b + 1e-8)
#   
#   n_YaYab[] <- rbinom(n_Ya_tr_fail[i], p_YZaYZab)
#   n_YaYac[] <- n_Ya_tr_fail[i] - n_YaYab[i]
#   n_YbYab[] <- rbinom(n_Yb_tr_fail[i], p_YZbYZab)
#   n_YbYbc[] <- n_Yb_tr_fail[i] - n_YbYab[i]
#   n_YcYac[] <- rbinom(n_Yc_tr_fail[i], p_YZcYZac)
#   n_YcYbc[] <- n_Yc_tr_fail[i] - n_YcYac[i]
#   
#   n_ZaZab[] <- rbinom(n_Za_tr_fail[i], p_YZaYZab)
#   n_ZaZac[] <- n_Za_tr_fail[i] - n_ZaZab[i]
#   n_ZbZab[] <- rbinom(n_Zb_tr_fail[i], p_YZbYZab)
#   n_ZbZbc[] <- n_Zb_tr_fail[i] - n_ZbZab[i]
#   n_ZcZac[] <- rbinom(n_Zc_tr_fail[i], p_YZcYZac)
#   n_ZcZbc[] <- n_Zc_tr_fail[i] - n_ZcZac[i]
#   
#   ##now move to the double resistance compartments
#   #first calculate total outflow -- calculate exit probabilities and then draw
#   p_YabOut <- 1 - exp(-(E_c*T_s + (1 - E_c)*kappa*T_sr + d))
#   p_YacOut <- 1 - exp(-(E_b*T_s + (1 - E_b)*kappa*T_sr + d))
#   p_YbcOut <- 1 - exp(-(E_a*T_s + (1 - E_a)*kappa*T_sr + d))
#   p_ZabOut <- 1 - exp(-(E_c*T_m + d))
#   p_ZacOut <- 1 - exp(-(E_b*T_m + d))
#   p_ZbcOut <- 1 - exp(-(E_a*T_m + d))
#   
#   n_YabOut[] <- rbinom(Y_ab[i], p_YabOut)
#   n_YacOut[] <- rbinom(Y_ac[i], p_YacOut)
#   n_YbcOut[] <- rbinom(Y_bc[i], p_YbcOut)
#   n_ZabOut[] <- rbinom(Z_ab[i], p_ZabOut)
#   n_ZacOut[] <- rbinom(Z_ac[i], p_ZacOut)
#   n_ZbcOut[] <- rbinom(Z_bc[i], p_ZbcOut)
#   
#   #because these are already dual resistance, the only ones leaving and not recovering are developing resistance to the same drug, so separate those and the rest recover
#   p_Yab_tr_fail <- (E_c*w_c*T_s)/(E_c*T_s + (1 - E_c)*kappa*T_sr + d)
#   p_Yac_tr_fail <- (E_b*w_b*T_s)/(E_b*T_s + (1 - E_b)*kappa*T_sr + d)
#   p_Ybc_tr_fail <- (E_a*w_a*T_s)/(E_a*T_s + (1 - E_a)*kappa*T_sr + d)
#   p_Zab_tr_fail <- (E_c*w_c*T_m)/(E_c*T_m + d)
#   p_Zac_tr_fail <- (E_b*w_b*T_m)/(E_b*T_m + d)
#   p_Zbc_tr_fail <- (E_a*w_a*T_m)/(E_a*T_m + d)
#   
#   #split out resistant flows to Yabc vs recoveries to S
#   n_YabYabc[] <- rbinom(n_YabOut[i], p_Yab_tr_fail)
#   n_YabS[] <- n_YabOut[i] - n_YabYabc[i]
#   n_YacYabc[] <- rbinom(n_YacOut[i], p_Yac_tr_fail)
#   n_YacS[] <- n_YacOut[i] - n_YacYabc[i]
#   n_YbcYabc[] <- rbinom(n_YbcOut[i], p_Ybc_tr_fail)
#   n_YbcS[] <- n_YbcOut[i] - n_YbcYabc[i]
#   
#   n_ZabZabc[] <- rbinom(n_ZabOut[i], p_Zab_tr_fail)
#   n_ZabS[] <- n_ZabOut[i] - n_ZabZabc[i]
#   n_ZacZabc[] <- rbinom(n_ZacOut[i], p_Zac_tr_fail)
#   n_ZacS[] <- n_ZacOut[i] - n_ZacZabc[i]
#   n_ZbcZabc[] <- rbinom(n_ZbcOut[i], p_Zbc_tr_fail)
#   n_ZbcS[] <- n_ZbcOut[i] - n_ZbcZabc[i]
#   
#   #recovery from triple resistance states
#   p_YabcOut <- 1 - exp(-(kappa*T_sr + d))
#   p_ZabcOut <- 1 - exp(-d)
#   n_YabcOut[] <- rbinom(Y_abc[i], p_YabcOut)
#   n_ZabcOut[] <- rbinom(Z_abc[i], p_ZabcOut)
#   
#   
#   #check whether prevalence to A has reached the threshold
#   threshold_met_A = (prevA >= threshold)
#   update(lock_threshold_A) = lock_threshold_A + threshold_met_A
#   
#   #check whether prevalence to B has reached the threshold
#   threshold_met_B = (prevB >= threshold)
#   update(lock_threshold_B) = lock_threshold_B + threshold_met_B
#   
#   #check whether prevalence to C has reached the threshold
#   threshold_met_C = (prevC >= threshold)
#   update(lock_threshold_C) = lock_threshold_C + threshold_met_C
#   
#   #Binary switches between treatment
#   #THIS PART IS UNIQUE TO THE MODEL STRATEGY
#   E_a <- (lock_threshold_A < 1)
#   E_c <- 1 - (lock_threshold_B < 1)
#   E_b <- 1  - E_a - E_c
#   
#   sum_tx_probs = E_a + E_c + E_b 
#   
#   #sum array variables for book-keeping / accounting outputs
#   sum_n_SOut <- sum(n_SOut)
#   sum_n_SYZa <- sum(n_SYZa)
#   sum_n_SYZb <- sum(n_SYZb)
#   sum_n_SYZc <- sum(n_SYZc)
#   sum_n_SYZab <- sum(n_SYZab)
#   sum_n_SYZac <- sum(n_SYZac)
#   sum_n_SYZbc <- sum(n_SYZbc)
#   sum_n_SYZabc <- sum(n_SYZabc)
#   sum_n_SYZs <- sum(n_SYZs)
#   sum_n_SYs <- sum(n_SYs)
#   sum_n_SZs <- sum(n_SZs)
#   sum_n_SYa <- sum(n_SYa)
#   sum_n_SZa <- sum(n_SZa)
#   sum_n_SYb <- sum(n_SYb)
#   sum_n_SZb <- sum(n_SZb)
#   sum_n_SYc <- sum(n_SYc)
#   sum_n_SZc <- sum(n_SZc)
#   sum_n_SYab <- sum(n_SYab)
#   sum_n_SZab <- sum(n_SZab)
#   sum_n_SYac <- sum(n_SYac)
#   sum_n_SZac <- sum(n_SZac)
#   sum_n_SYbc <- sum(n_SYbc)
#   sum_n_SZbc <- sum(n_SZbc)
#   sum_n_SYabc <- sum(n_SYabc)
#   sum_n_SZabc <- sum(n_SZabc)
#   sum_n_YsOut <- sum(n_YsOut)
#   sum_n_ZsOut <- sum(n_ZsOut)
#   sum_n_Ys_tr <- sum(n_Ys_tr)
#   sum_n_Ys_nr <- sum(n_Ys_nr)
#   sum_n_Zs_tr <- sum(n_Zs_tr)
#   sum_n_Zs_nr <- sum(n_Zs_nr)
#   sum_n_Ys_tr_fail <- sum(n_Ys_tr_fail)
#   sum_n_Ys_tr_success <- sum(n_Ys_tr_success)
#   sum_n_Zs_tr_fail <- sum(n_Zs_tr_fail)
#   sum_n_Zs_tr_success <- sum(n_Zs_tr_success)
#   sum_n_YsS <- sum(n_YsS)
#   sum_n_ZsS <- sum(n_ZsS)
#   sum_n_YsYa <- sum(n_YsYa)
#   sum_n_YsYb <- sum(n_YsYb)
#   sum_n_YsYc <- sum(n_YsYc)
#   sum_n_ZsZa <- sum(n_ZsZa)
#   sum_n_ZsZb <- sum(n_ZsZb)
#   sum_n_ZsZc <- sum(n_ZsZc)
#   sum_n_YaOut <- sum(n_YaOut)
#   sum_n_YbOut <- sum(n_YbOut)
#   sum_n_YcOut <- sum(n_YcOut)
#   sum_n_ZaOut <- sum(n_ZaOut)
#   sum_n_ZbOut <- sum(n_ZbOut)
#   sum_n_ZcOut <- sum(n_ZcOut)
#   sum_n_Ya_nr <- sum(n_Ya_nr)
#   sum_n_Ya_tr <- sum(n_Ya_tr)
#   sum_n_Yb_nr <- sum(n_Yb_nr)
#   sum_n_Yb_tr <- sum(n_Yb_tr)
#   sum_n_Yc_nr <- sum(n_Yc_nr)
#   sum_n_Yc_tr <- sum(n_Yc_tr)
#   sum_n_Za_nr <- sum(n_Za_nr)
#   sum_n_Za_tr <- sum(n_Za_tr)
#   sum_n_Zb_nr <- sum(n_Zb_nr)
#   sum_n_Zb_tr <- sum(n_Zb_tr)
#   sum_n_Zc_nr <- sum(n_Zc_nr)
#   sum_n_Zc_tr <- sum(n_Zc_tr)
#   sum_n_Ya_tr_fail <- sum(n_Ya_tr_fail)
#   sum_n_Yb_tr_fail <- sum(n_Yb_tr_fail)
#   sum_n_Yc_tr_fail <- sum(n_Yc_tr_fail)
#   sum_n_Za_tr_fail <- sum(n_Za_tr_fail)
#   sum_n_Zb_tr_fail <- sum(n_Zb_tr_fail)
#   sum_n_Zc_tr_fail <- sum(n_Zc_tr_fail)
#   sum_n_Ya_tr_success <- sum(n_Ya_tr_success)
#   sum_n_Yb_tr_success <- sum(n_Yb_tr_success)
#   sum_n_Yc_tr_success <- sum(n_Yc_tr_success)
#   sum_n_Za_tr_success <- sum(n_Za_tr_success)
#   sum_n_Zb_tr_success <- sum(n_Zb_tr_success)
#   sum_n_Zc_tr_success <- sum(n_Zc_tr_success)
#   sum_n_YaS <- sum(n_YaS)
#   sum_n_YbS <- sum(n_YbS)
#   sum_n_YcS <- sum(n_YcS)
#   sum_n_ZaS <- sum(n_ZaS)
#   sum_n_ZbS <- sum(n_ZbS)
#   sum_n_ZcS <- sum(n_ZcS)
#   sum_n_YaYab <- sum(n_YaYab)
#   sum_n_YaYac <- sum(n_YaYac)
#   sum_n_YbYab <- sum(n_YbYab)
#   sum_n_YbYbc <- sum(n_YbYbc)
#   sum_n_YcYac <- sum(n_YcYac)
#   sum_n_YcYbc <- sum(n_YcYbc)
#   sum_n_ZaZab <- sum(n_ZaZab)
#   sum_n_ZaZac <- sum(n_ZaZac)
#   sum_n_ZbZab <- sum(n_ZbZab)
#   sum_n_ZbZbc <- sum(n_ZbZbc)
#   sum_n_ZcZac <- sum(n_ZcZac)
#   sum_n_ZcZbc <- sum(n_ZcZbc)
#   sum_n_YabOut <- sum(n_YabOut)
#   sum_n_YacOut <- sum(n_YacOut)
#   sum_n_YbcOut <- sum(n_YbcOut)
#   sum_n_ZabOut <- sum(n_ZabOut)
#   sum_n_ZacOut <- sum(n_ZacOut)
#   sum_n_ZbcOut <- sum(n_ZbcOut)
#   sum_n_YabYabc <- sum(n_YabYabc)
#   sum_n_YabS <- sum(n_YabS)
#   sum_n_YacYabc <- sum(n_YacYabc)
#   sum_n_YacS <- sum(n_YacS)
#   sum_n_YbcYabc <- sum(n_YbcYabc)
#   sum_n_YbcS<- sum(n_YbcS)
#   sum_n_ZabZabc <- sum(n_ZabZabc)
#   sum_n_ZabS <- sum(n_ZabS)
#   sum_n_ZacZabc <- sum(n_ZacZabc)
#   sum_n_ZacS <- sum(n_ZacS)
#   sum_n_ZbcZabc <- sum(n_ZbcZabc)
#   sum_n_ZbcS<- sum(n_ZbcS)
#   sum_n_YabcOut <- sum(n_YabcOut)
#   sum_n_ZabcOut <- sum(n_ZabcOut)
#   sum_n_SOut_remaining <- sum(n_SOut_remaining)
#   sum_n_SOut_remaining2 <- sum(n_SOut_remaining2)
#   sum_p_SOut_remaining2 <- sum(p_SOut_remaining2)
#   sum_n_SOut_remaining3 <- sum(n_SOut_remaining3)
#   sum_p_SOut_remaining3 <- sum(p_SOut_remaining3)
#   sum_n_SOut_remaining4 <- sum(n_SOut_remaining4)
#   sum_p_SOut_remaining4 <- sum(p_SOut_remaining4)
#   sum_n_SOut_remaining5 <- sum(n_SOut_remaining5)
#   sum_p_SOut_remaining5 <- sum(p_SOut_remaining5)
#   sum_n_SOut_remaining6 <- sum(n_SOut_remaining6)
#   sum_p_SOut_remaining6 <- sum(p_SOut_remaining6)
#   
#   
#   
#   #incidence tracking
#   Inc_cat[] = sum(BI_S[,i])*S[i]
#   Inc <- sum(Inc_cat)
#   ninf_Ys[] <- sigma*sum(BI_YZs[,i])*S[i]
#   ninf_Zs[] <- (1-sigma)*sum(BI_YZs[,i])*S[i]
#   ninf_Ya[] <- sigma*f_a*sum(BI_YZa[,i])*S[i]
#   ninf_Za[] <- (1 - sigma)*f_a*sum(BI_YZa[,i])*S[i]
#   ninf_Yb[] <- sigma*f_b*sum(BI_YZb[,i])*S[i]
#   ninf_Zb[] <- (1 - sigma)*f_b*sum(BI_YZb[,i])*S[i]
#   ninf_Yab[] <- sigma*f_ab*sum(BI_YZab[,i])*S[i]
#   ninf_Zab[] <- (1-sigma)*f_ab*sum(BI_YZab[,i])*S[i]
#   
#   #set initial conditions
#   initial(S[]) <- S_ini[i]
#   initial(Y_s[]) <- Y_s_ini[i]
#   initial(Z_s[]) <- Z_s_ini[i]
#   initial(Y_a[]) <- Y_a_ini[i]
#   initial(Z_a[]) <- Z_a_ini[i]
#   initial(Y_b[]) <- Y_b_ini[i]
#   initial(Z_b[]) <- Z_b_ini[i]
#   initial(Y_c[]) <- Y_c_ini[i]
#   initial(Z_c[]) <- Z_c_ini[i]
#   initial(Y_ab[]) <- Y_ab_ini[i]
#   initial(Z_ab[]) <- Z_ab_ini[i]
#   initial(Y_ac[]) <- Y_ac_ini[i]
#   initial(Z_ac[]) <- Z_ac_ini[i]  
#   initial(Y_bc[]) <- Y_bc_ini[i]
#   initial(Z_bc[]) <- Z_bc_ini[i]   
#   initial(Y_abc[]) <- Y_abc_ini[i]
#   initial(Z_abc[]) <- Z_abc_ini[i]
#   initial(lock_threshold_A) <- 0
#   initial(lock_threshold_B) <- 0
#   initial(lock_threshold_C) <- 0
#   
#   #user-defined states
#   S_ini[] <- user()
#   Y_s_ini[] <- user()
#   Z_s_ini[] <- user()
#   Y_a_ini[] <- user()
#   Z_a_ini[] <- user()
#   Y_b_ini[] <- user()
#   Z_b_ini[] <- user()
#   Y_c_ini[] <- user()
#   Z_c_ini[] <- user()
#   Y_ab_ini[] <- user()
#   Z_ab_ini[] <- user()
#   Y_ac_ini[] <- user()
#   Z_ac_ini[] <- user()
#   Y_bc_ini[] <- user()
#   Z_bc_ini[] <- user() 
#   Y_abc_ini[] <- user()
#   Z_abc_ini[] <- user()
#   beta[,] <- user()
#   
#   #user-defined parameters
#   threshold <- user()
#   w_a <- user()
#   w_b <- user()
#   w_c <- user()
#   T_s <- user()
#   T_m <- user()
#   T_sr <- user()
#   kappa <- user()
#   f_a <- user()
#   f_b <- user()
#   f_c <- user()
#   f_ab <- user()
#   f_ac <- user()
#   f_bc <- user()
#   f_abc <- user()
#   sigma <- user()
#   N_risk <- user()
#   d <- user()
#   
#   #define variable dimensions
#   ##state initial values
#   dim(S_ini) <- N_risk
#   dim(Y_s_ini) <- N_risk
#   dim(Z_s_ini) <- N_risk
#   dim(Y_a_ini) <- N_risk
#   dim(Z_a_ini) <- N_risk
#   dim(Y_b_ini) <- N_risk
#   dim(Z_b_ini) <- N_risk
#   dim(Y_c_ini) <- N_risk
#   dim(Z_c_ini) <- N_risk
#   dim(Y_ab_ini) <- N_risk
#   dim(Z_ab_ini) <- N_risk
#   dim(Y_ac_ini) <- N_risk
#   dim(Z_ac_ini) <- N_risk
#   dim(Y_bc_ini) <- N_risk
#   dim(Z_bc_ini) <- N_risk
#   dim(Y_abc_ini) <- N_risk
#   dim(Z_abc_ini) <- N_risk
#   dim(sum_BI_S) <- N_risk
#   dim(sum_BI_YZs) <- N_risk
#   dim(sum_BI_YZa) <- N_risk
#   dim(sum_BI_YZb) <- N_risk
#   dim(sum_BI_YZc) <- N_risk
#   dim(sum_BI_YZab) <- N_risk
#   dim(sum_BI_YZac) <- N_risk
#   dim(sum_BI_YZbc) <- N_risk
#   dim(sum_BI_YZabc) <- N_risk
#   dim(sum_all_BI_terms) <- N_risk
#   dim(ninf_Ys) <- N_risk
#   dim(ninf_Zs) <- N_risk
#   dim(ninf_Ya) <- N_risk
#   dim(ninf_Za) <- N_risk
#   dim(ninf_Yb) <- N_risk
#   dim(ninf_Zb) <- N_risk
#   dim(ninf_Yab) <- N_risk
#   dim(ninf_Zab) <- N_risk
#   dim(n_SOut) <- N_risk
#   dim(p_SOut) <- N_risk
#   dim(p_SYZa) <- N_risk
#   dim(p_SYZb) <- N_risk
#   dim(p_SYZc) <- N_risk
#   dim(p_SYZab) <- N_risk
#   dim(p_SYZac) <- N_risk  
#   dim(p_SYZbc) <- N_risk
#   dim(p_SYZabc) <- N_risk
#   dim(n_SYZa) <- N_risk
#   dim(n_SYZb) <- N_risk
#   dim(n_SYZc) <- N_risk
#   dim(n_SYZab) <- N_risk
#   dim(n_SYZac) <- N_risk
#   dim(n_SYZbc) <- N_risk
#   dim(n_SYZabc) <- N_risk
#   dim(n_SYZs) <- N_risk
#   dim(n_SYs) <- N_risk
#   dim(n_SZs) <- N_risk
#   dim(n_SYa) <- N_risk
#   dim(n_SZa) <- N_risk
#   dim(n_SYb) <- N_risk
#   dim(n_SZb) <- N_risk
#   dim(n_SYc) <- N_risk
#   dim(n_SZc) <- N_risk
#   dim(n_SYab) <- N_risk
#   dim(n_SZab) <- N_risk
#   dim(n_SYac) <- N_risk
#   dim(n_SZac) <- N_risk
#   dim(n_SYbc) <- N_risk
#   dim(n_SZbc) <- N_risk
#   dim(n_SYabc) <- N_risk
#   dim(n_SZabc) <- N_risk
#   dim(n_YsOut) <- N_risk
#   dim(n_ZsOut) <- N_risk
#   dim(n_Ys_tr) <- N_risk
#   dim(n_Ys_nr) <- N_risk
#   dim(n_Zs_tr) <- N_risk
#   dim(n_Zs_nr) <- N_risk
#   dim(n_Ys_tr_fail) <- N_risk
#   dim(n_Ys_tr_success) <- N_risk
#   dim(n_Zs_tr_fail) <- N_risk
#   dim(n_Zs_tr_success) <- N_risk
#   dim(n_YsS) <- N_risk
#   dim(n_ZsS) <- N_risk
#   dim(n_YsYa) <- N_risk
#   dim(n_YsYb) <- N_risk
#   dim(n_YsYc) <- N_risk
#   dim(n_ZsZa) <- N_risk
#   dim(n_ZsZb) <- N_risk
#   dim(n_ZsZc) <- N_risk
#   dim(n_YaOut) <- N_risk
#   dim(n_YbOut) <- N_risk
#   dim(n_YcOut) <- N_risk
#   dim(n_ZaOut) <- N_risk
#   dim(n_ZbOut) <- N_risk
#   dim(n_ZcOut) <- N_risk
#   dim(n_Ya_nr) <- N_risk
#   dim(n_Ya_tr) <- N_risk
#   dim(n_Yb_nr) <- N_risk
#   dim(n_Yb_tr) <- N_risk
#   dim(n_Yc_nr) <- N_risk
#   dim(n_Yc_tr) <- N_risk
#   dim(n_Za_nr) <- N_risk
#   dim(n_Za_tr) <- N_risk
#   dim(n_Zb_nr) <- N_risk
#   dim(n_Zb_tr) <- N_risk
#   dim(n_Zc_nr) <- N_risk
#   dim(n_Zc_tr) <- N_risk
#   dim(n_Ya_tr_fail) <- N_risk
#   dim(n_Yb_tr_fail) <- N_risk
#   dim(n_Yc_tr_fail) <- N_risk
#   dim(n_Za_tr_fail) <- N_risk
#   dim(n_Zb_tr_fail) <- N_risk
#   dim(n_Zc_tr_fail) <- N_risk
#   dim(n_Ya_tr_success) <- N_risk
#   dim(n_Yb_tr_success) <- N_risk
#   dim(n_Yc_tr_success) <- N_risk
#   dim(n_Za_tr_success) <- N_risk
#   dim(n_Zb_tr_success) <- N_risk
#   dim(n_Zc_tr_success) <- N_risk
#   dim(n_YaS) <- N_risk
#   dim(n_YbS) <- N_risk
#   dim(n_YcS) <- N_risk
#   dim(n_ZaS) <- N_risk
#   dim(n_ZbS) <- N_risk
#   dim(n_ZcS) <- N_risk
#   dim(n_YaYab) <- N_risk
#   dim(n_YaYac) <- N_risk
#   dim(n_YbYab) <- N_risk
#   dim(n_YbYbc) <- N_risk
#   dim(n_YcYac) <- N_risk
#   dim(n_YcYbc) <- N_risk
#   dim(n_ZaZab) <- N_risk
#   dim(n_ZaZac) <- N_risk
#   dim(n_ZbZab) <- N_risk
#   dim(n_ZbZbc) <- N_risk
#   dim(n_ZcZac) <- N_risk
#   dim(n_ZcZbc) <- N_risk
#   dim(n_YabOut) <- N_risk
#   dim(n_YacOut) <- N_risk
#   dim(n_YbcOut) <- N_risk
#   dim(n_ZabOut) <- N_risk
#   dim(n_ZacOut) <- N_risk
#   dim(n_ZbcOut) <- N_risk
#   dim(n_YabYabc) <- N_risk
#   dim(n_YabS) <- N_risk
#   dim(n_YacYabc) <- N_risk
#   dim(n_YacS) <- N_risk
#   dim(n_YbcYabc) <- N_risk
#   dim(n_YbcS) <- N_risk
#   dim(n_ZabZabc) <- N_risk
#   dim(n_ZabS) <- N_risk
#   dim(n_ZacZabc) <- N_risk
#   dim(n_ZacS) <- N_risk
#   dim(n_ZbcZabc) <- N_risk
#   dim(n_ZbcS) <- N_risk
#   dim(n_YabcOut) <- N_risk
#   dim(n_ZabcOut) <- N_risk
#   dim(n_SOut_remaining) <- N_risk
#   dim(p_SOut_remaining) <- N_risk
#   dim(n_SOut_remaining2) <- N_risk
#   dim(p_SOut_remaining2) <- N_risk
#   dim(n_SOut_remaining3) <- N_risk
#   dim(p_SOut_remaining3) <- N_risk
#   dim(n_SOut_remaining4) <- N_risk
#   dim(p_SOut_remaining4) <- N_risk
#   dim(n_SOut_remaining5) <- N_risk
#   dim(p_SOut_remaining5) <- N_risk
#   dim(n_SOut_remaining6) <- N_risk
#   dim(p_SOut_remaining6) <- N_risk
#   
#   ##states
#   dim(S) <- N_risk
#   dim(Y_s) <- N_risk
#   dim(Z_s) <- N_risk
#   dim(Y_a) <- N_risk
#   dim(Z_a) <- N_risk
#   dim(Y_b) <- N_risk
#   dim(Z_b) <- N_risk
#   dim(Y_c) <- N_risk
#   dim(Z_c) <- N_risk
#   dim(Y_ab) <- N_risk
#   dim(Z_ab) <- N_risk
#   dim(Y_ac) <- N_risk
#   dim(Z_ac) <- N_risk
#   dim(Y_bc) <- N_risk
#   dim(Z_bc) <- N_risk
#   dim(Y_abc) <- N_risk
#   dim(Z_abc) <- N_risk
#   ##other vector / matrix variables
#   dim(beta) <- c(N_risk, N_risk)
#   dim(N) <- N_risk
#   dim(BI_S) <- c(N_risk, N_risk)
#   dim(BI_YZs) <- c(N_risk, N_risk)
#   dim(BI_YZa) <- c(N_risk, N_risk)
#   dim(BI_YZb) <- c(N_risk, N_risk)
#   dim(BI_YZc) <- c(N_risk, N_risk)
#   dim(BI_YZab) <- c(N_risk, N_risk)
#   dim(BI_YZac) <- c(N_risk, N_risk)
#   dim(BI_YZbc) <- c(N_risk, N_risk)
#   dim(BI_YZabc) <- c(N_risk, N_risk)
#   dim(Inc_cat) <- N_risk
#   
#   
#   #define desired non-state outputs
#   output(ninf_Ys) <- ninf_Ys
#   output(ninf_Zs) <- ninf_Zs
#   output(ninf_Ya) <- ninf_Ya
#   output(ninf_Za) <- ninf_Za
#   output(ninf_Yb) <- ninf_Yb
#   output(ninf_Zb) <- ninf_Zb
#   output(ninf_Yab) <- ninf_Yab
#   output(ninf_Zab) <- ninf_Zab
#   output(N) <- N
#   output(prevA) <- prevA
#   output(prevB) <- prevB
#   output(prevC) <- prevC
#   output(Inc_cat) <- Inc_cat
#   output(Inc) <- Inc
#   output(E_a) <- E_a
#   output(E_b) <- E_b
#   output(E_c) <- E_c
#   output(all_a_res) <- all_a_res
#   output(all_b_res) <- all_b_res
#   output(all_c_res) <- all_c_res
#   output(all_non_res) <- all_non_res
#   output(total_sum_BI_S) <- total_sum_BI_S
#   output(total_sum_BI_YZs) <- total_sum_BI_YZs
#   output(total_sum_BI_YZa) <- total_sum_BI_YZa
#   output(total_sum_BI_YZb) <- total_sum_BI_YZb
#   output(total_sum_BI_YZc) <- total_sum_BI_YZc
#   output(total_sum_BI_YZab) <- total_sum_BI_YZab
#   output(total_sum_BI_YZac) <- total_sum_BI_YZac
#   output(total_sum_BI_YZbc) <- total_sum_BI_YZbc
#   output(total_sum_BI_YZabc) <- total_sum_BI_YZabc
#   output(total_sum_all_BI_terms) <- total_sum_all_BI_terms
#   output(sum_S) <- sum_S
#   output(sum_Y_s) <- sum_Y_s
#   output(sum_Y_a) <- sum_Y_a
#   output(sum_Y_b) <- sum_Y_b
#   output(sum_Y_c) <- sum_Y_c
#   output(sum_Y_ab) <- sum_Y_ab
#   output(sum_Y_ac) <- sum_Y_ac
#   output(sum_Y_bc) <- sum_Y_bc
#   output(sum_Y_abc) <- sum_Y_abc
#   output(sum_Z_s) <- sum_Z_s
#   output(sum_Z_a) <- sum_Z_a
#   output(sum_Z_b) <- sum_Z_b
#   output(sum_Z_c) <- sum_Z_c
#   output(sum_Z_ab) <- sum_Z_ab
#   output(sum_Z_ac) <- sum_Z_ac
#   output(sum_Z_bc) <- sum_Z_bc
#   output(sum_Z_abc) <- sum_Z_abc
#   output(sum_YZ_s) <- sum_YZ_s
#   output(sum_YZ_a) <- sum_YZ_a
#   output(sum_YZ_b) <- sum_YZ_b
#   output(sum_YZ_c) <- sum_YZ_c
#   output(sum_YZ_ab) <- sum_YZ_ab
#   output(sum_YZ_ac) <- sum_YZ_ac
#   output(sum_YZ_bc) <- sum_YZ_bc
#   output(sum_YZ_abc) <- sum_YZ_abc
#   output(sum_N) <- sum_N
#   output(prev) <- prev
#   output(sum_tx_probs) <- sum_tx_probs
#   output(p_SOut) <- p_SOut
#   output(sum_n_SOut) <- sum_n_SOut
#   output(p_SYZa) <- p_SYZa
#   output(p_SYZb) <- p_SYZb
#   output(p_SYZc) <- p_SYZc
#   output(p_SYZab) <- p_SYZab
#   output(p_SYZac) <- p_SYZac
#   output(p_SYZbc) <- p_SYZbc
#   output(p_SYZabc) <- p_SYZabc
#   output(sum_n_SYZa) <- sum_n_SYZa
#   output(sum_n_SYZb) <- sum_n_SYZb
#   output(sum_n_SYZc) <- sum_n_SYZc
#   output(sum_n_SYZab) <- sum_n_SYZab
#   output(sum_n_SYZac) <- sum_n_SYZac
#   output(sum_n_SYZbc) <- sum_n_SYZbc
#   output(sum_n_SYZabc) <- sum_n_SYZabc
#   output(sum_n_SYZs) <- sum_n_SYZs
#   output(sum_n_SYs) <- sum_n_SYs
#   output(sum_n_SZs) <- sum_n_SZs
#   output(sum_n_SYa) <- sum_n_SYa
#   output(sum_n_SZa) <- sum_n_SZa
#   output(sum_n_SYb) <- sum_n_SYb
#   output(sum_n_SZb) <- sum_n_SZb
#   output(sum_n_SYc) <- sum_n_SYc
#   output(sum_n_SZc) <- sum_n_SZc
#   output(sum_n_SYab) <- sum_n_SYab
#   output(sum_n_SZab) <- sum_n_SZab
#   output(sum_n_SYac) <- sum_n_SYac
#   output(sum_n_SZac) <- sum_n_SZac
#   output(sum_n_SYbc) <- sum_n_SYbc
#   output(sum_n_SZbc) <- sum_n_SZbc
#   output(sum_n_SYabc) <- sum_n_SYabc
#   output(sum_n_SZabc) <- sum_n_SZabc
#   output(p_YsOut) <- p_YsOut
#   output(p_ZsOut) <- p_ZsOut
#   output(p_Ys_tr) <- p_Ys_tr
#   output(p_Zs_tr) <- p_Zs_tr
#   output(p_YZs_tr_fail) <- p_YZs_tr_fail
#   output(p_YZs_YZa) <- p_YZs_YZa
#   output(p_YZs_YZb) <- p_YZs_YZb
#   output(p_YZs_YZc) <- p_YZs_YZc
#   output(p_YaOut) <- p_YaOut
#   output(p_YbOut) <- p_YbOut
#   output(p_YcOut) <- p_YcOut
#   output(p_ZaOut) <- p_ZaOut
#   output(p_ZbOut) <- p_ZbOut
#   output(p_ZcOut) <- p_ZcOut
#   output(p_Ya_nr) <- p_Ya_nr
#   output(p_Yb_nr) <- p_Yb_nr
#   output(p_Yc_nr) <- p_Yc_nr
#   output(p_Za_nr) <- p_Za_nr
#   output(p_Zb_nr) <- p_Zb_nr
#   output(p_Zc_nr) <- p_Zc_nr
#   output(p_Ya_tr_fail) <- p_Ya_tr_fail
#   output(p_Yb_tr_fail) <- p_Yb_tr_fail
#   output(p_Yc_tr_fail) <- p_Yc_tr_fail
#   output(p_Za_tr_fail) <- p_Za_tr_fail
#   output(p_Zb_tr_fail) <- p_Zb_tr_fail
#   output(p_Zc_tr_fail) <- p_Zc_tr_fail
#   output(p_YZaYZab) <- p_YZaYZab
#   output(p_YZbYZab) <- p_YZbYZab
#   output(p_YZcYZac) <- p_YZcYZac
#   output(p_YabOut) <- p_YabOut 
#   output(p_YacOut) <- p_YacOut 
#   output(p_YbcOut) <- p_YbcOut 
#   output(p_ZabOut) <- p_ZabOut 
#   output(p_ZacOut) <- p_ZacOut 
#   output(p_ZbcOut) <- p_ZbcOut
#   output(p_Yab_tr_fail) <- p_Yab_tr_fail
#   output(p_Yac_tr_fail) <- p_Yac_tr_fail
#   output(p_Ybc_tr_fail) <- p_Ybc_tr_fail
#   output(p_Zab_tr_fail) <- p_Zab_tr_fail
#   output(p_Zac_tr_fail) <- p_Zac_tr_fail
#   output(p_Zbc_tr_fail) <- p_Zbc_tr_fail
#   output(sum_n_YsOut) <- sum_n_YsOut
#   output(sum_n_ZsOut) <- sum_n_ZsOut
#   output(sum_n_Ys_tr) <- sum_n_Ys_tr
#   output(sum_n_Ys_nr) <- sum_n_Ys_nr
#   output(sum_n_Zs_tr) <- sum_n_Zs_tr
#   output(sum_n_Zs_nr) <- sum_n_Zs_nr
#   output(sum_n_Ys_tr_fail) <- sum_n_Ys_tr_fail
#   output(sum_n_Ys_tr_success) <- sum_n_Ys_tr_success
#   output(sum_n_Zs_tr_fail) <- sum_n_Zs_tr_fail
#   output(sum_n_Zs_tr_success) <- sum_n_Zs_tr_success
#   output(sum_n_YsS) <- sum_n_YsS
#   output(sum_n_ZsS) <- sum_n_ZsS
#   output(sum_n_YsYa) <- sum_n_YsYa
#   output(sum_n_YsYb) <- sum_n_YsYb
#   output(sum_n_YsYc) <- sum_n_YsYc
#   output(sum_n_ZsZa) <- sum_n_ZsZa
#   output(sum_n_ZsZb) <- sum_n_ZsZb
#   output(sum_n_ZsZc) <- sum_n_ZsZc
#   output(sum_n_YaOut) <- sum_n_YaOut
#   output(sum_n_YbOut) <- sum_n_YbOut
#   output(sum_n_YcOut) <- sum_n_YcOut
#   output(sum_n_ZaOut) <- sum_n_ZaOut
#   output(sum_n_ZbOut) <- sum_n_ZbOut
#   output(sum_n_ZcOut) <- sum_n_ZcOut
#   output(sum_n_Ya_nr) <- sum_n_Ya_nr
#   output(sum_n_Ya_tr) <- sum_n_Ya_tr
#   output(sum_n_Yb_nr) <- sum_n_Yb_nr
#   output(sum_n_Yb_tr) <- sum_n_Yb_tr
#   output(sum_n_Yc_nr) <- sum_n_Yc_nr
#   output(sum_n_Yc_tr) <- sum_n_Yc_tr
#   output(sum_n_Za_nr) <- sum_n_Za_nr
#   output(sum_n_Za_tr) <- sum_n_Za_tr
#   output(sum_n_Zb_nr) <- sum_n_Zb_nr
#   output(sum_n_Zb_tr) <- sum_n_Zb_tr
#   output(sum_n_Zc_nr) <- sum_n_Zc_nr
#   output(sum_n_Zc_tr) <- sum_n_Zc_tr
#   output(sum_n_Ya_tr_fail) <- sum_n_Ya_tr_fail
#   output(sum_n_Yb_tr_fail) <- sum_n_Yb_tr_fail
#   output(sum_n_Yc_tr_fail) <- sum_n_Yc_tr_fail
#   output(sum_n_Za_tr_fail) <- sum_n_Za_tr_fail
#   output(sum_n_Zb_tr_fail) <- sum_n_Zb_tr_fail
#   output(sum_n_Zc_tr_fail) <- sum_n_Zc_tr_fail
#   output(sum_n_Ya_tr_success) <- sum_n_Ya_tr_success
#   output(sum_n_Yb_tr_success) <- sum_n_Yb_tr_success
#   output(sum_n_Yc_tr_success) <- sum_n_Yc_tr_success
#   output(sum_n_Za_tr_success) <- sum_n_Za_tr_success
#   output(sum_n_Zb_tr_success) <- sum_n_Zb_tr_success
#   output(sum_n_Zc_tr_success) <- sum_n_Zc_tr_success
#   output(sum_n_YaS) <- sum_n_YaS
#   output(sum_n_YbS) <- sum_n_YbS
#   output(sum_n_YcS) <- sum_n_YcS
#   output(sum_n_ZaS) <- sum_n_ZaS
#   output(sum_n_ZbS) <- sum_n_ZbS
#   output(sum_n_ZcS) <- sum_n_ZcS
#   output(sum_n_YaYab) <- sum_n_YaYab
#   output(sum_n_YaYac) <- sum_n_YaYac
#   output(sum_n_YbYab) <- sum_n_YbYab
#   output(sum_n_YbYbc) <- sum_n_YbYbc
#   output(sum_n_YcYac) <- sum_n_YcYac 
#   output(sum_n_YcYbc) <- sum_n_YcYbc 
#   output(sum_n_ZaZab) <- sum_n_ZaZab
#   output(sum_n_ZaZac) <- sum_n_ZaZac
#   output(sum_n_ZbZab) <- sum_n_ZbZab
#   output(sum_n_ZbZbc) <- sum_n_ZbZbc
#   output(sum_n_ZcZac) <- sum_n_ZcZac 
#   output(sum_n_ZcZbc) <- sum_n_ZcZbc 
#   output(sum_n_YabOut) <- sum_n_YabOut 
#   output(sum_n_YacOut) <- sum_n_YacOut 
#   output(sum_n_YbcOut) <- sum_n_YbcOut
#   output(sum_n_ZabOut) <- sum_n_ZabOut 
#   output(sum_n_ZacOut) <- sum_n_ZacOut 
#   output(sum_n_ZbcOut) <- sum_n_ZbcOut
#   output(sum_n_YabYabc) <- sum_n_YabYabc
#   output(sum_n_YabS) <- sum_n_YabS
#   output(sum_n_YacYabc) <- sum_n_YacYabc
#   output(sum_n_YacS) <- sum_n_YacS
#   output(sum_n_YbcYabc) <- sum_n_YbcYabc
#   output(sum_n_YbcS) <- sum_n_YbcS
#   output(sum_n_ZabZabc) <- sum_n_ZabZabc
#   output(sum_n_ZabS) <- sum_n_ZabS
#   output(sum_n_ZacZabc) <- sum_n_ZacZabc
#   output(sum_n_ZacS) <- sum_n_ZacS
#   output(sum_n_ZbcZabc) <- sum_n_ZbcZabc
#   output(sum_n_ZbcS) <- sum_n_ZbcS
#   output(sum_n_ZabcOut) <- sum_n_ZabcOut 
#   output(p_YabcOut) <- p_YabcOut
#   output(p_ZabcOut) <- p_ZabcOut
#   output(sum_n_YabcOut) <- sum_n_YabcOut
#   output(sum_n_SOut_remaining) <- sum_n_SOut_remaining
#   output(sum_n_SOut_remaining2) <- sum_n_SOut_remaining2
#   output(sum_p_SOut_remaining2) <- sum_p_SOut_remaining2
#   output(sum_n_SOut_remaining3) <- sum_n_SOut_remaining3
#   output(sum_p_SOut_remaining3) <- sum_p_SOut_remaining3
#   output(sum_n_SOut_remaining4) <- sum_n_SOut_remaining4
#   output(sum_p_SOut_remaining4) <- sum_p_SOut_remaining4
#   output(sum_n_SOut_remaining5) <- sum_n_SOut_remaining5
#   output(sum_p_SOut_remaining5) <- sum_p_SOut_remaining5
#   output(sum_n_SOut_remaining6) <- sum_n_SOut_remaining6
#   output(sum_p_SOut_remaining6) <- sum_p_SOut_remaining6
# })







