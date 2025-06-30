#Define parameters -- from 6/12 calibration
b <- 0.5167924 #transmission probability per partnership
w_a <- 10^-4 #rate of resistance to A emerging on treatment with A #this is a probability, over course of treatment
w_b <- 10^-5 #rate of resistance to B emerging on treatment with B
w_c <- 10^-6 #rate of resistance to C emerging on treatment with C
T_s <- 1/(14.9435648) #time to treatment for symptomatic infection
T_m <- 0.4150241/365  #screening rate (time to treatment for asymptomatic infection)
T_sr <- T_s/3 #time to retreatment for symptomatic infection, if failure ##rethink this?
kappa <- 0.90 #probability of retreatment after treatment failure, symptomatic
f_a <- 0.98 #relative fitness when resistant to A vs pan-susceptible 
f_b <- 0.95 #relative fitness when resistant to B vs pan-susceptible
f_c <- 0.95 #relative fitness when resistant to C vs pan-susceptible
f_ab <- f_a*f_b #relative fitness when resistant to A&B vs pan-susceptible
f_ac <- f_a*f_c #relative fitness when resistant to A&C vs pan-susceptible
f_bc <- f_b*f_c #relative fitness when resistant to B&C vs pan-susceptible
f_abc <- f_a*f_b*f_c #relative fitness when resistant to A,B&C vs pan-susceptible
sigma <- 0.4936291  #probability of symptomatic infection
N_risk <- 3 #number of risk strata
d <- 1/(91.8422808) #natural recovery rate from infection #this is probably shorter, Barbee study
resA <- 0.0001 #initial prev of resistance to A 
resB <- 0 
resC <- 0.0001 #initial prev of resistance to C ##GOING TO CHANGE THIS LATER FOR BASELINE SCENARIOS
pop <- 10^6 #pop size
pop.p <- c(0.3, 0.6, 0.1) #relative size of each risk group; low, M, high
c_min <- 1.3785266 #rate of partner change in lowest risk group
activities <- c(1*c_min/365, 
                5*c_min/365, 
                20*c_min/365)  #sexual contacts per day
epsilon <- 0.2357691 #mixing parameter (for non-assortativity)

#risk stratified mixing 
beta <- ((1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + epsilon*activities/(pop*pop.p)*diag(3))*b
years <- 150
tstep <- 1


#set NEW initial conditions: VERY IMPORTANT this is different from above (but probably won't change all that much)
#initial conditions and initial condition calculations
resA <- 0.0001 #initial prev of resistance to A 
resB <- 0 
resC <- 0

N_ini <- c(pop*pop.p[1],
           pop*pop.p[2],
           pop*pop.p[3])
#distribute GC cases to have overall 3% prevalence
x <- 0.03/(.3*0.029+.6*0.154+.1*0.817)
gc_lo <- round(N_ini[1]*x*0.029,0)
gc_md <- round(N_ini[2]*x*0.154,0)
gc_hi <- round(N_ini[3]*x*0.817,0)

#estimated proportion of symp. infections at cross-sectional point in time, @ equilibrium
# estimated from calibration_ER script
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

Y_b_ini <- round(c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resB)
Z_b_ini <- round(c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resB)

Y_c_ini <- round(c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resC)
Z_c_ini <- round(c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resC)

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

dt <- seq(0,365*150, tstep) 

sum(Y_s_ini + Z_s_ini + Z_a_ini) + sum(S_ini)


baseline_parms_6_12 <- list(S_ini = S_ini,
                           Y_s_ini = Y_s_ini,
                           Z_s_ini = Z_s_ini,
                           Y_a_ini = Y_a_ini,
                           Z_a_ini = Z_a_ini,
                           Y_b_ini = Y_b_ini,
                           Z_b_ini = Z_b_ini,
                           Y_c_ini = Y_c_ini,
                           Z_c_ini = Z_c_ini, 
                           Y_ab_ini = Y_ab_ini,
                           Z_ab_ini = Z_ab_ini,
                           Y_ac_ini = Y_ac_ini,
                           Z_ac_ini = Z_ac_ini,
                           Y_bc_ini = Y_bc_ini,
                           Z_bc_ini = Z_bc_ini, 
                           Y_abc_ini = Y_abc_ini,
                           Z_abc_ini = Z_abc_ini,
                           beta = beta,
                           w_a = w_a, 
                           w_b = w_b,
                           w_c = w_c,
                           T_s = T_s,
                           T_m = T_m,
                           T_sr = T_sr,
                           kappa = kappa,
                           f_a = f_a,
                           f_b = f_b,
                           f_c = f_c,
                           f_ab = f_ab,
                           f_ac = f_ac,
                           f_bc = f_bc,
                           f_abc = f_abc,
                           sigma = sigma,
                           N_risk = N_risk,
                           d = d,
                           threshold = 0.05
)


















