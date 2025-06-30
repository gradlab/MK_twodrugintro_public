#SCRIPT FROM REICHERT ET AL 2023 WITH MODIFICATIONS#
library(deSolve)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readr)
library(stringr)
library(ggplot2)
install.packages("bbmle")
library(bbmle)

#logit and inverse logit functions for processing output
logit<-function(x) {log(x/(1-x))}
ilogit <-function(x) {1/(1+exp(-x))}

# calculate GC prevalence based on current parameters, for calibration
pred_fun_er <- function(params){
  params2 <- list(pop=10^6, pop.p=c(0.3, 0.6, 0.1), rho = 1/(20*365), 
                  prA = 1, omega_a = 10^-4, fA = 0.98, #updated mutation probability for A from ER script to match baseline parameters
                  pi_s = 0.90)
  parms <- c(params, params2)
  years = 2 #can run out for more years but values don't change, gets to equilibrium quickly 
  dt <- seq(0, 365*years, 1) 
  calibration_sim <- as.data.frame(ode(inits, dt, calibration.SI, parms = parms))   
  calibration_sim <- round(calibration_sim, 0)
  end <- calibration_sim[nrow(calibration_sim),]
  #Overall prevalence of gonorrhea at t=end years
  prev_GC <- sum(end[,5:16])/sum(end[,2:16])
  prev_GC
  return(prev_GC)
}

# estimate alpha and beta for our beta distribution of prevalence of GC
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# Calibration Model (Drug A Only)
calibration.SI <- function(t, x, parms){
  with(as.list(c(t, x, parms)),{
    N = c(N1, N2, N3) #total population
    S = c(S1, S2, S3) #not infected
    Y0 = c(Y01, Y02, Y03) #symptomatic infected, no resistance
    Ya = c(Ya1, Ya2, Ya3) #symptomatic infected, resistant to A
    Z0 = c(Z01, Z02, Z03) #asymptomatic infected, no resistance
    Za = c(Za1, Za2, Za3) #asymptomatic infected, resistant to A
    #makes 3x3 contact matrix
    activities <- c(1*c_min/365, 5*c_min/365, 20*c_min/365)
    beta <- (1-epsilon)*outer(activities, activities)/sum(pop*pop.p*activities) + 
      epsilon*activities/(pop*pop.p)*diag(3)
    beta <- beta * b #contacts * transmission pr per partnership
    #susceptibles
    dS <- -beta %*% ((Y0 + Z0) + fA*(Ya+Za)) *S + 
      (1-omega_a*prA)*(Ts*Y0 + Tm*Z0) +
      pi_s*Tsr*Ya +
      g*(Y0 + Ya + Z0 + Za) #removed entry and exit from ER script to match 3 drug models
    #infections w/ no resistance
    dY0 <- sigma*(beta %*% (Y0+Z0)*S) - Ts*Y0 - g*Y0 
    dZ0 <- (1-sigma)*(beta %*% (Y0+Z0)*S) - Tm*Z0 - g*Z0 
    #infections w/ resistance to A
    dYa <- sigma*fA*(beta %*% (Ya+Za)*S) +
      omega_a*prA*Ts*Y0 - 
      pi_s*Tsr*Ya -
      g*Ya 
    dZa <- (1-sigma)*fA*(beta %*% (Ya+Za)*S) +
      omega_a*prA*Tm*Z0 -
      g*Za 
    der <- c(dS, dY0, dZ0, dYa, dZa)
    list(der)
  })
}

#SET PARAM VALUES
#parameters that we will estimate from prior lit (not MLE)
pop = 10^6                                           #pop size
pop.p = c(0.3, 0.6, 0.1)                             #relative size of each risk group; low, M, high
omega_a = 10^-4      #Pr of emergence of resistance on treatment with A (ceftriaxone) #Updated from ER script
prA = 1              #Pr of treatment with A
fA = 0.98            #relative fitness, resistant to A
pi_s = 0.9          #Pr of retreatment if initial treatment failure, symptomatic
resA = 0.0001        #initial prev of resistance to A

#Set model duration + initial conditions

#how long to run the model?
years = 2 #this was 2 -- changing to higher doesn't make a difference
tstep = 1 #in days

#set N for each sexual risk group (1 = low risk, 2 = med risk, 3 = hi risk)
N <- pop
N1 <- N*pop.p[1]
N2 <- N*pop.p[2]
N3 <- N*pop.p[3]

prev_target <- 0.03 #Can Update to desired GC prevalence target 

#distribute GC cases to have overall 3% prevalence
x <- prev_target/(.3*0.029+.6*0.154+.1*0.817)

gc_lo <- round(N1*x*0.029,0)
gc_md <- round(N2*x*0.154,0)
gc_hi <- round(N3*x*0.817,0)

#estimated proportion of symp. infections at cross-sectional point in time
prev_symp <- 0.089 #from Tuite et al. 2017 code
Y_gc_lo <- gc_lo*prev_symp
Z_gc_lo <- gc_lo*(1-prev_symp)
Y_gc_md <- gc_md*prev_symp
Z_gc_md <- gc_md*(1-prev_symp)
Y_gc_hi <- gc_hi*prev_symp
Z_gc_hi <- gc_hi*(1-prev_symp)

Y0 <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * (1-resA)
Z0 <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * (1-resA)

YA <- c(Y_gc_lo, Y_gc_md, Y_gc_hi) * resA
ZA <- c(Z_gc_lo, Z_gc_md, Z_gc_hi) * resA

inits <- round(c(S1 = N1-gc_lo, S2 = N2-gc_md, S3 = N3-gc_hi, 
           Y01 = Y0[1], Y02 = Y0[2], Y03 = Y0[3],
           Z01 = Z0[1], Z02 = Z0[2], Z03 = Z0[3], 
           Ya1 = YA[1], Ya2 = YA[2], Ya3 = YA[3], 
           Za1 = ZA[1], Za2 = ZA[2], Za3 = ZA[3]) )

#create vector w/ all timepoints
dt <- seq(0, 365*years, tstep) 


#start with parameter values from Tuite et al. 2017
theta <- c(logit.b=logit(0.5),  #transmission probability
           log.c_min=log(1.2165433),  #min rate of partner change
           logit.epsilon=logit(0.2427385), #epsilon
           logit.sigma = logit(0.50), #pr of incident symptomatic infection
           log.Ts=log(0.041*365), #duration of infectiousness if symptomatic (days) = average time to treatment
           log.g=log(105), #duration of infectiousness if asymptomatic and untreated (days)
           logit.Tm =logit(0.39)) #screening rate per year


### set parameters ###
b <- ilogit(theta["logit.b"])
c_min <- exp(theta["log.c_min"]) 
epsilon <- ilogit(theta["logit.epsilon"])
sigma <- ilogit(theta["logit.sigma"])
Ts <- exp(theta["log.Ts"])/365
g <- exp(theta["log.g"])/365
Tm <- ilogit(theta["logit.Tm"])


#######################################
#### Maximum likelihood estimation ####
#######################################
#used for estimating parameters

model.epi.loglik <- function(theta) {
  b <- ilogit(theta["logit.b"])
  c_min <- exp(theta["log.c_min"]) 
  epsilon <- ilogit(theta["logit.epsilon"])
  sigma <- ilogit(theta["logit.sigma"])
  Ts <- 1/(exp(theta["log.Ts"]))
  g <- 1/(exp(theta["log.g"]))
  Tm <- ilogit(theta["logit.Tm"])/365
  params <-list(c_min = c_min, epsilon = epsilon, sigma = sigma, b=b,Ts = Ts,Tm = Tm, g = g, Tsr = Ts/3)
  pred <- pred_fun_er(params)
  beta.params.prev <- estBetaParams(mu = pred, var = 1.47e-5)
  ll <- sum(dbeta(x= prev_target, beta.params.prev$alpha, beta.params.prev$beta, log=TRUE)) #calculate likelihood
  ll[is.na(ll)]<-(-1e20)
  print(ll)
  c(prev=pred,ll=ll)
}

f.optim <-function(theta) {
  Res <-  model.epi.loglik(theta)
  LogLL <- Res["ll"] #Model is returning the Log likelihood
  return(-LogLL)
}


values.start=theta
parnames(f.optim)<-names(values.start)
fit0 <- bbmle::mle2(f.optim, start=values.start,  vecpar=TRUE,  optimizer="optim"); fit0
fit <- mle2(f.optim, start=coef(fit0),  vecpar=TRUE, optimizer="optim"); fit
theta.fit<-coef(fit)

exp(theta.fit)
ilogit(theta.fit)


write_csv(data.frame("parms" = names(exp(theta.fit)), "values" = exp(theta.fit)), "../output/calibration/6_23_cluster_calibration_expparms.csv")
write_csv(data.frame("parms" = names(ilogit(theta.fit)), "values" = ilogit(theta.fit)), "../output/calibration/6_23_cluster_calibration_ilogitparms.csv")



