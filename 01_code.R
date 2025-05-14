################################################################################
# Code Description
################################################################################
### Yajima and Nakashima (2025)
### "A hierarchical Bayesian approach for estimating the number of groups and 
###  group sizes in group-living animals using passive detectors"
###
###
### Supporting information
###
### This script describes the NIMBLE code of the ANGST model used to analyse the field data.
### Parameters and variables in the code are given essentially the same names as in the main text
###############################################################################

# library
library(nimble)
library(MCMCvis)
library(bayesplot)

################################################################################
# NIMBLE code
################################################################################
Model <- nimbleCode(
  {
    # Priors of scalar parameters -----------------------------------------------
    # GN model
    mu  ~ dgamma(0.01, 0.01)
    tau ~ dunif(0.0, 1.0)

    # AS model
    gamma ~ dgamma(1, 1)
    gamma2 ~ dgamma(1, 1)

    # GS model
    lambda ~ dgamma(0.01, 0.01)
    phi    ~ dbeta(1.0, 1.0)
    rho    ~ dbeta(1.0, 1.0)

    a <- phi       * ((1/rho) - 1)
    b <- (1 - phi) * ((1/rho) - 1)

    # GN model (Royle-Nichols model; Royle and Nichols, 2003) ------------------
    for(i in 1:R){
      K[i] ~ dpois(mu)
      P[i] <- 1 - ((1 - tau)^K[i])
      Y[i] ~ dbin(P[i], S)
    }

    # LTSB prior ---------------------------------------------------------------
    ### 'Robs': Number of detectors with one or more group observations (Y[i]>0)
    ### The detectors (detector index is i = {1,...,i,...Robs,...,R}) are sorted so that Y[i]>0 up to 1~Robs and Y[i]=0 up to Robs+1~R.
    for(i in 1:Robs){
      # Generation of the breaking point
      for(m in 1:(M-1)){
        v[i,m] ~ dbeta(1.0, gamma * z[i,m] + (1-z[i,m]) * gamma2)
      }
      for (m in 1:M) {
        z[i,m] <- step(K[i]-m) # Indicator vector, see Data analysis section in the main text
      }
      # TSB prior
      w[i,1:M]      <- stick_breaking(v[i, 1:(M-1)])
      theta0[i,1:M] <- w[i, 1:M] * z[i, 1:M]
      theta[i,1:M]  <- theta0[i,1:M] / sum(theta0[i,1:M])
    }
    # Poisson beta-binomial model & Group assignment ---------------------------
    for(i in 1:Robs) {
      for(m in 1:M) {
        N0[i,m] ~ T(dpois(lambda), 1,) # True group size
      }
    }
    ### A modifications from Yajima and Nakashima (2024) indexing
    ### - In the main text, we used indexes detector i and passage l (e.g, the variables were represented by a matrix g[i,l] (i=1,...,Robs; l = 1,...,L[i]))
    ### - However variables corresponding to the passes were made into a single vector in this script
    ### - The length of this vector is 'nsmp', nsmp <- sum(L[1:Robs])
    ### - We now introduce a new vector St[n]
    ### - This vector contains the index of the detector i corresponding to the n_th observation (e.g. a detector (i=1) observed two groups, a detector (i=2) observed three groups and ...)
    ### - 'St' is then shown as follows: St <- c(1,1,2,2,2,...,Robs)
    for(n in 1:nsmp){
      # Assignment of true group size to observation groups
      g[n] ~ dcat(theta[St[n], 1:M])
      N[n] <- N0[St[n], g[n]]

      # Modelling of detection probability
      p[n] ~ dbeta(a, b)

      # Observed group size
      C[n] ~ T(dbin(p[n], N[n]), 1,)

      # Bayesian-P value, following Kery and Royle (2015)
      eval[n]  <- (p[n]*N[n])/(1-(1-p[n])^N[n])
      C.new[n]  ~ T(dbin(p[n], N[n]), 1,)
      E[n]     <- pow((C[n]-eval[n]),2) / (eval[n] + 0.5)
      E.new[n] <- pow((C.new[n]-eval[n]),2) / (eval[n] + 0.5)
    }

    # Statistics ---------------------------------------------------------------
    Lambda      <- lambda/(1-exp(-lambda))  # Mean group size

    # variables 'fit' and 'fit.new' are for calculating Bayesian-P value
    fit     <- sum(E[1:nsmp])
    fit.new <- sum(E.new[1:nsmp])
  }#end
)

################################################################################
# Analysis of ANGST model by NIMBLE
################################################################################
# Load field data
# This data is camera trapping data of juvenile wild boars (Sus scrofa) on the Boso Peninsula in Japan. 

# - C   : Observed group size (vectorized)
# - L   : The number of group passes at detector i
# - M   : The maximum number of observations by all the detectors plus one
# - R   : The number of detectors (camera traps)
# - Y   : The number of sampling occasions in which at least one group was detected at detector i
# - S   : The total sampling occasions
# - Robs: The number of detectors (Y[i]>0)
# - St  : The index vector of the detector i corresponding to the n_th observation (see NIMBLE code)
setwd("../../project/ANGST/")
load("case_study.Rdata")

nsmp <- sum(L) # The number of passes in all detectors
Cmax <- c()    # The maximum observation group size for each detector (this vector is used to determine the initial value of MCMC. See below)
for (i in 1:max(St)) {
  Cmax[i] <- max(C[St==i])
}
# NIMBLE data
data <- list(C = C, Y = Y)
constants <- list(M = M, R = R, S = S, Robs = Robs, nsmp = nsmp, St = St)

# Parameters to be monitored
params <- c("fit", "fit.new", "tau", "gamma", "Lambda", "phi", "rho", "g")


# Initial values
inits <- function(){
  v_init <- array(NA, dim = c(Robs, M-1))
  N_init <- array(NA, dim = c(Robs, M))
  z_init <- array(NA, dim = c(Robs, M))

  for (i in 1:Robs) {
    v_init[i,] <- rbeta(M-1, 1.0, 1.0)
    ### For C ~ dbin(p,N), C <= N must always be satisfied,
    ### so the initial value of N0[i,1:M] was initially chosen to be 0 to 3 larger than 
    ### the maximum observation group size (Cmax[i]) for detector i.
    N_init[i,] <- Cmax[i]+sample(c(0, 1, 2), M, replace = TRUE)
    z_init[i,] <- rep(1, M)
  }
  list(
    # GN model
    mu  = runif(1, 2.0, 6.0),
    K   = rep(1, R),
    P   = runif(R, 0.2, 0.9),
    tau = runif(1, 0.0, 1.0),

    # AS model
    g      = rep(1, nsmp),
    v      = v_init,
    gamma  = runif(1, 0.2, 0.7),
    gamma2 = runif(1, 0.2, 0.7),
    z      = z_init,

    # GS model
    lambda = runif(1, 2.0, 6.0),
    N0     = N_init,
    a      = runif(1, 1, 5),
    b      = runif(1, 1, 5),
    phi    = runif(1, 0.4, 0.9),
    rho    = runif(1, 0.4, 0.9),
    p      = runif(nsmp, 0.4, 0.9),

    # Bayesian-P value
    C.new = C
  )
}

# Setting for the MCMC
ni <- 200000 # Number of iterations
nb <- 100000 # Number of initial, pre-thinning, MCMC iterations to discard
nt <- 50     # Thinning interval for collecting MCMC samples
nc <- 3     # Number of MCMC chains

MCMCsteps <- ((ni-nb)*nc)/nt
sprintf("Resulting in %d samples", MCMCsteps)

# run
# Processing time is approximately 14 hours (5 hours if parallelized)
# out <- nimbleMCMC(
#   Model,
#   data = data, constants = constants,
#   inits = inits, monitors = params,
#   niter = ni, nburnin = nb, thin = nt,nchains = nc,
#   samplesAsCodaMCMC = TRUE
# )

myModel <- nimbleModel(
  code = Model,
  data = data,
  constants = constants,
  inits = inits()
)

CmyModel    <- compileNimble(myModel)
configModel <- configureMCMC(myModel, monitors = params)

configModel$removeSamplers(c("gamma", "mu", "tau", "phi", "rho", "lambda"))
# change samplers
configModel$addSampler(c("gamma"), type = "slice")
configModel$addSampler(c("mu", "tau"), type = "AF_slice")
configModel$addSampler(c("phi", "rho"), type = "AF_slice")
configModel$addSampler(c("lambda"), type = "slice")

myMCMC <- buildMCMC(configModel, monitors = params)
CmyMCMC <- compileNimble(myMCMC)

seed <- rep(712, nc)
out  <- runMCMC(CmyMCMC, niter = ni, nburnin = nb, thin = nt, nchains = nc, setSeed = seed, samplesAsCodaMCMC = TRUE)

# save(out, file = "case_study_results.Rdata")

# results of posterior distributions
MCMCsummary(object = out, round = 2, params = params[-length(params)])

# mu^(detected) samples
K_detected_samples <- array(NA, dim = c(MCMCsteps, Robs))
for(i in 1:Robs) {
  det_index  <- which(St == i)
  if(length(det_index) > 1) {
    K_detected_samples[, i] <- apply(out[[1]][, paste0("g[", det_index, "]")], 1, function(x) length(unique(x)))
    K_detected_samples[, i] <- apply(out[[1]][, paste0("g[", det_index, "]")], 1, function(x) length(unique(x)))
  }else{
    K_detected_samples[, i] <- rep(1, MCMCsteps)
    K_detected_samples[, i] <- rep(1, MCMCsteps)
  }
}
mu_detected_samples <- apply(K_detected_samples, 1, mean)

# mu^(detected) results
# mean
mean(mu_detected_samples)
quantile(mu_detected_samples, 0.500)
quantile(mu_detected_samples, 0.025)
quantile(mu_detected_samples, 0.975)


# Trace plot & Density plot
MCMCtrace(
  object = out,
  pdf = FALSE, # no export to PDF
  ind = TRUE, # separate density lines per chain
  params = c("mu", "tau", "gamma", "Lambda", "phi", "rho")
)

# Bayesian-P value
fits <- c(out$chain1[,"fit"],out$chain2[,"fit"], out$chain3[,"fit"])
fits_new <- c(out$chain1[,"fit.new"],out$chain2[,"fit.new"], out$chain3[,"fit.new"])
sprintf("Bayesian-p value is = %.2f", sum(fits>fits_new)/length(fits))

# end of the script