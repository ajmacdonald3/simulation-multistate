################################################################################
# SEASONAL SURVIVAL SIMULATION SCENARIOS
#
# trying to improve parameter identifiability
################################################################################

library(rjags)
library(jagsUI)
library(tidyverse)
library(viridis)
library(cowplot)

set.seed(42)

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
  # Unobservable: number of state that is unobservable
  n.occasions <- dim(PSI.STATE)[4] + 1
  CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
  g <- colSums(marked)
  for (s in 1:dim(PSI.STATE)[1]){
    if (g[s]==0) next  # To avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
  } #s
  for (i in 1:sum(marked)){
    for (s in 1:dim(PSI.STATE)[1]){
      if (mark.occ[i,s]==0) next
      first <- mark.occ[i,s]
      CH[i,first] <- s
      CH.TRUE[i,first] <- s
    } #s
    for (t in (first+1):n.occasions){
      # Multinomial trials for state transitions
      if (first==n.occasions) next
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t] <- state
      # Multinomial trials for observation process
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t] <- event
    } #t
  } #i
  # Replace the NA and the highest state number (dead) in the file by 0
  CH[is.na(CH)] <- 0
  CH[CH==dim(PSI.STATE)[1]] <- 0
  CH[CH==unobservable] <- 0
  id <- numeric(0)
  for (i in 1:dim(CH)[1]){
    z <- min(which(CH[i,]!=0))
    ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
  }
  return(list(CH=CH[-id,], CH.TRUE=CH.TRUE[-id,]))
  # CH: capture histories to be used
  # CH.TRUE: capture histories with perfect observation
}

known.state.ms <- function(ms, notseen){
  # notseen: label for “not seen”
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}

################################################################################

#### SCENARIO 9 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: temporal random effect, following fieldwork schedule
# psi: following annual cycle and time-dependent with impossible transitions set to 0

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15 # 5 years
n.states <- 6
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93
phiUN <- 0.92

pDB <- c(0, 0, 0.5, 0, 0, 0.6, 0, 0, 0.7, 0, 0, 0.6, 0, 0)
pJB <- c(0, 0, 0, 0.6, 0, 0, 0.4, 0, 0, 0.7, 0, 0, 0.6, 0)
pMI <- c(0.5, 0, 0, 0.4, 0, 0, 0.8, 0, 0, 0.6, 0, 0, 0.3, 0)
pAR <- c(0, 0.4, 0, 0, 0.8, 0, 0, 0.5, 0, 0, 0.3, 0, 0, 0.5)

psiDB.DB <- rep(0, n.occasions)
psiDB.JB <- rep(c(0.4, 0, 0), n.occasions/3)
psiDB.MI <- rep(c(0.4, 0, 0), n.occasions/3)
psiDB.AR <- rep(0, n.occasions)

psiJB.DB <- rep(0, n.occasions)
psiJB.JB <- rep(0, n.occasions)
psiJB.MI <- rep(0, n.occasions)
psiJB.AR <- rep(c(0, 0.7, 0), n.occasions/3)

psiMI.DB <- rep(0, n.occasions)
psiMI.JB <- rep(0, n.occasions)
psiMI.MI <- rep(0, n.occasions)
psiMI.AR <- rep(c(0, 0.6, 0), n.occasions/3)

psiAR.DB <- rep(c(0, 0, 0.8), n.occasions/3)
psiAR.JB <- rep(0, n.occasions)
psiAR.MI <- rep(0, n.occasions)
psiAR.AR <- rep(0, n.occasions)

psiUN.DB <- rep(c(0, 0, 0.8), n.occasions/3)
psiUN.JB <- rep(c(0.4, 0, 0), n.occasions/3)
psiUN.MI <- rep(c(0.4, 0, 0), n.occasions/3)
psiUN.AR <- rep(c(0, 0.7, 0), n.occasions/3)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(700, 0, 0), n.occasions/3)
marked[,2] <- rep(c(0, 50, 0), n.occasions/3)
marked[,3] <- rep(c(0, 100, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 100), n.occasions/3)
marked[,5] <- rep(0, n.occasions)
marked[,6] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiDB*psiDB.DB[t], phiDB*psiDB.JB[t], phiDB*psiDB.MI[t], phiDB*psiDB.AR[t], phiDB*(1-psiDB.JB[t]-psiDB.MI[t]),                         1-phiDB,
      phiJB*psiJB.DB[t], phiJB*psiJB.JB[t], phiJB*psiJB.MI[t], phiJB*psiJB.AR[t], phiJB*(1-psiJB.AR[t]),                                     1-phiJB,
      phiMI*psiMI.DB[t], phiMI*psiMI.JB[t], phiMI*psiMI.MI[t], phiMI*psiMI.AR[t], phiMI*(1-psiMI.AR[t]),                                     1-phiMI,
      phiAR*psiAR.DB[t], phiAR*psiAR.JB[t], phiAR*psiAR.MI[t], phiAR*psiAR.AR[t], phiAR*(1-psiAR.DB[t]),                                     1-phiAR,
      phiUN*psiUN.DB[t], phiUN*psiUN.JB[t], phiUN*psiUN.MI[t], phiUN*psiUN.AR[t], phiUN*(1-psiUN.DB[t]-psiUN.JB[t]-psiUN.MI[t]-psiUN.AR[t]), 1-phiUN,
      0,                 0,                 0,                 0,                 0,                                                         1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pDB[t], 0,      0,      0,      1-pDB[t],
      0,      pJB[t], 0,      0,      1-pJB[t],
      0,      0,      pMI[t], 0,      1-pMI[t],
      0,      0,      0,      pAR[t], 1-pAR[t],
      0,      0,      0,      0,      1,
      0,      0,      0,      0,      1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = 5)
CH <- sim$CH
CH.TRUE <- sim$CH.TRUE

# Compute vector with occasions of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in DB, 2 = seen alive in JB, 3 = seen alive in MI, 4 = seen alive in AR, 5 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 5

sink("simulation_scenarios.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiDB: survival probability at Delaware Bay
# phiJB: survival probability at James Bay
# phiMI: survival probability at Mingan
# phiAR: survival probability in Argentina
# phiUN: survival probability when unobservable

# psiDB.JB: movement probability from Delaware Bay to James Bay
# psiDB.MI: movement probability from Delaware Bay to James Bay
# psiDB.UN: movement probability from Delaware Bay to unobservable state (1-psiDB.JB-psiDB.MI)
# psiJB.AR: movement probability from James Bay to Argentina
# psiJB.UN: movement probability from James Bay to unobservable state (1-psiJB.AR)
# psiMI.AR: movement probability from Mingan to Argentina
# psiMI.UN: movement probability from Mingan to unobservable state (1-psiMI.AR)
# psiAR.DB: movement probability from Argentina to Delaware Bay
# psiAR.UN: movement probability from Argentina to unobservable state (1-psiAR.DB)
# psiUN.DB: movement probability from unobservable to Delaware Bay
# psiUN.JB: movement probability from unobservable to James Bay
# psiUN.MI: movement probability from unobservable to Mingan
# psiUN.AR: movement probability from unobservable to Argentina
# all other movement probabilities set to zero because not biologically plausible

# pDB: resight probability at Delaware Bay
# pJB: resight probability at James Bay
# pMI: resight probability at Mingan
# pAR: resight probability in Argentina
# -------------------------------------------------
# States (S):
# 1 alive at Delaware Bay
# 2 alive at James Bay
# 3 alive at Mingan
# 4 alive in Argentina
# 5 alive but unobservable
# 6 dead

# Observations (O):
# 1 seen at Delaware Bay 
# 2 seen at James Bay
# 3 seen at Mingan
# 4 seen in Argentina
# 5 not seen
# -------------------------------------------------

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phi[1]
   phiJB[t] <- mean.phi[2]
   phiMI[t] <- mean.phi[3]
   phiAR[t] <- mean.phi[4]
   phiUN[t] <- mean.phi[5]
   }
   
   for (u in 1:5){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   pDB[1] <- 0
   pDB[2] <- 0
   pDB[3] <- ptDB[1]
   pDB[4] <- 0
   pDB[5] <- 0
   pDB[6] <- ptDB[2]
   pDB[7] <- 0
   pDB[8] <- 0
   pDB[9] <- ptDB[3]
   pDB[10] <- 0
   pDB[11] <- 0
   pDB[12] <- ptDB[4]
   pDB[13] <- 0
   pDB[14] <- 0
   
   for (u in 1:4){
   logit(ptDB[u]) <- muDB + epsilonDB[u]
   epsilonDB[u] ~ dnorm(0, tauDB)T(-15, 15)
   }
  
   muDB <- log(mean.pDB / (1-mean.pDB))     # Logit transformation
   mean.pDB ~ dunif(0, 1)                   # Prior for mean survival
   tauDB <- pow(sigmaDB, -2)
   sigmaDB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2DB <- pow(sigmaDB, 2)              # Temporal variance
   
   pJB[1] <- 0
   pJB[2] <- 0
   pJB[3] <- 0
   pJB[4] <- ptJB[1]
   pJB[5] <- 0
   pJB[6] <- 0
   pJB[7] <- ptJB[2]
   pJB[8] <- 0
   pJB[9] <- 0
   pJB[10] <- ptJB[3]
   pJB[11] <- 0
   pJB[12] <- 0
   pJB[13] <- ptJB[4]
   pJB[14] <- 0
   
   for (u in 1:4){
   logit(ptJB[u]) <- muJB + epsilonJB[u]
   epsilonJB[u] ~ dnorm(0, tauJB)T(-15, 15)
   }
   
   muJB <- log(mean.pJB / (1-mean.pJB))     # Logit transformation
   mean.pJB ~ dunif(0, 1)                   # Prior for mean survival
   tauJB <- pow(sigmaJB, -2)
   sigmaJB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2JB <- pow(sigmaJB, 2)              # Temporal variance
   
   pMI[1] <- ptMI[1]
   pMI[2] <- 0
   pMI[3] <- 0
   pMI[4] <- ptMI[2]
   pMI[5] <- 0
   pMI[6] <- 0
   pMI[7] <- ptMI[3]
   pMI[8] <- 0
   pMI[9] <- 0
   pMI[10] <- ptMI[4]
   pMI[11] <- 0
   pMI[12] <- 0
   pMI[13] <- ptMI[5]
   pMI[14] <- 0
   
   for (u in 1:5){
   logit(ptMI[u]) <- muMI + epsilonMI[u]
   epsilonMI[u] ~ dnorm(0, tauMI)T(-15, 15)
   }
   
   muMI <- log(mean.pMI / (1-mean.pMI))     # Logit transformation
   mean.pMI ~ dunif(0, 1)                   # Prior for mean survival
   tauMI <- pow(sigmaMI, -2)
   sigmaMI ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2MI <- pow(sigmaMI, 2)              # Temporal variance
   
   pAR[1] <- 0
   pAR[2] <- ptAR[1]
   pAR[3] <- 0
   pAR[4] <- 0
   pAR[5] <- ptAR[2]
   pAR[6] <- 0
   pAR[7] <- 0
   pAR[8] <- ptAR[3]
   pAR[9] <- 0
   pAR[10] <- 0
   pAR[11] <- ptAR[4]
   pAR[12] <- 0
   pAR[13] <- 0
   pAR[14] <- ptAR[5]
   
   for (u in 1:5){
   logit(ptAR[u]) <- muAR + epsilonAR[u]
   epsilonAR[u] ~ dnorm(0, tauAR)T(-15, 15)
   }
   
   muAR <- log(mean.pAR / (1-mean.pAR))     # Logit transformation
   mean.pAR ~ dunif(0, 1)                   # Prior for mean survival
   tauAR <- pow(sigmaAR, -2)
   sigmaAR ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2AR <- pow(sigmaAR, 2)              # Temporal variance

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probs
      for (i in 1:4){
         lpsiDB[i] ~ dnorm(0, 0.001)
         lpsiJB[i] ~ dnorm(0, 0.001)
         lpsiMI[i] ~ dnorm(0, 0.001)
         lpsiAR[i] ~ dnorm(0, 0.001)
         lpsiUN[i] ~ dnorm(0, 0.001)
      }

      # Constrain the transitions such that their sum is < 1
      # Delaware Bay
      for (t in 1:(n.occasions-1)){
      psiDB[1,t] <- 0
      psiDB[4,t] <- 0
      }
      
      psiDB[2,1] <- exp(lpsiDB[2]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[2,2] <- 0
      psiDB[2,3] <- 0
      psiDB[2,4] <- exp(lpsiDB[2]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[2,5] <- 0
      psiDB[2,6] <- 0
      psiDB[2,7] <- exp(lpsiDB[2]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[2,8] <- 0
      psiDB[2,9] <- 0
      psiDB[2,10] <- exp(lpsiDB[2]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[2,11] <- 0
      psiDB[2,12] <- 0
      psiDB[2,13] <- exp(lpsiDB[2]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[2,14] <- 0
      #psiDB[2,15] <- 0

      psiDB[3,1] <- exp(lpsiDB[3]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[3,2] <- 0
      psiDB[3,3] <- 0
      psiDB[3,4] <- exp(lpsiDB[3]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[3,5] <- 0
      psiDB[3,6] <- 0
      psiDB[3,7] <- exp(lpsiDB[3]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[3,8] <- 0
      psiDB[3,9] <- 0
      psiDB[3,10] <- exp(lpsiDB[3]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[3,11] <- 0
      psiDB[3,12] <- 0
      psiDB[3,13] <- exp(lpsiDB[3]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
      psiDB[3,14] <- 0
      #psiDB[3,15] <- 0
      
      # James Bay
      for (t in 1:(n.occasions-1)){
      psiJB[1,t] <- 0
      psiJB[2,t] <- 0
      psiJB[3,t] <- 0
      }
      
      psiJB[4,1] <- 0
      psiJB[4,2] <- exp(lpsiJB[4]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
      psiJB[4,3] <- 0
      psiJB[4,4] <- 0
      psiJB[4,5] <- exp(lpsiJB[4]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
      psiJB[4,6] <- 0
      psiJB[4,7] <- 0
      psiJB[4,8] <- exp(lpsiJB[4]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
      psiJB[4,9] <- 0
      psiJB[4,10] <- 0
      psiJB[4,11] <- exp(lpsiJB[4]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
      psiJB[4,12] <- 0
      psiJB[4,13] <- 0
      psiJB[4,14] <- exp(lpsiJB[4]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
      #psiJB[4,15] <- 0
      
      # Mingan
      for (t in 1:(n.occasions-1)){
      psiMI[1,t] <- 0
      psiMI[2,t] <- 0
      psiMI[3,t] <- 0
      }
      
      psiMI[4,1] <- 0
      psiMI[4,2] <- exp(lpsiMI[4]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
      psiMI[4,3] <- 0
      psiMI[4,4] <- 0
      psiMI[4,5] <- exp(lpsiMI[4]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
      psiMI[4,6] <- 0
      psiMI[4,7] <- 0
      psiMI[4,8] <- exp(lpsiMI[4]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
      psiMI[4,9] <- 0
      psiMI[4,10] <- 0
      psiMI[4,11] <- exp(lpsiMI[4]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
      psiMI[4,12] <- 0
      psiMI[4,13] <- 0
      psiMI[4,14] <- exp(lpsiMI[4]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
      #psiMI[4,15] <- 0

      # Argentina
      psiAR[1,1] <- 0
      psiAR[1,2] <- 0
      psiAR[1,3] <- exp(lpsiAR[1]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
      psiAR[1,4] <- 0
      psiAR[1,5] <- 0
      psiAR[1,6] <- exp(lpsiAR[1]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
      psiAR[1,7] <- 0
      psiAR[1,8] <- 0
      psiAR[1,9] <- exp(lpsiAR[1]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
      psiAR[1,10] <- 0
      psiAR[1,11] <- 0
      psiAR[1,12] <- exp(lpsiAR[1]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
      psiAR[1,13] <- 0
      psiAR[1,14] <- 0
      #psiAR[1,15] <- exp(lpsiAR[1]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
      
      for (t in 1:(n.occasions-1)){
      psiAR[2,t] <- 0
      psiAR[3,t] <- 0
      psiAR[4,t] <- 0
      }
      
      # Unobservable state
      psiUN[1,1] <- 0
      psiUN[1,2] <- 0
      psiUN[1,3] <- exp(lpsiUN[1]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[1,4] <- 0
      psiUN[1,5] <- 0
      psiUN[1,6] <- exp(lpsiUN[1]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[1,7] <- 0
      psiUN[1,8] <- 0
      psiUN[1,9] <- exp(lpsiUN[1]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[1,10] <- 0
      psiUN[1,11] <- 0
      psiUN[1,12] <- exp(lpsiUN[1]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[1,13] <- 0
      psiUN[1,14] <- 0
      #psiUN[1,15] <- exp(lpsiUN[1]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      
      psiUN[2,1] <- exp(lpsiUN[2]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[2,2] <- 0
      psiUN[2,3] <- 0
      psiUN[2,4] <- exp(lpsiUN[2]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[2,5] <- 0
      psiUN[2,6] <- 0
      psiUN[2,7] <- exp(lpsiUN[2]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[2,8] <- 0
      psiUN[2,9] <- 0
      psiUN[2,10] <- exp(lpsiUN[2]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[2,11] <- 0
      psiUN[2,12] <- 0
      psiUN[2,13] <- exp(lpsiUN[2]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[2,14] <- 0
      #psiUN[2,15] <- 0

      psiUN[3,1] <- exp(lpsiUN[3]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[3,2] <- 0
      psiUN[3,3] <- 0
      psiUN[3,4] <- exp(lpsiUN[3]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[3,5] <- 0
      psiUN[3,6] <- 0
      psiUN[3,7] <- exp(lpsiUN[3]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[3,8] <- 0
      psiUN[3,9] <- 0
      psiUN[3,10] <- exp(lpsiUN[3]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[3,11] <- 0
      psiUN[3,12] <- 0
      psiUN[3,13] <- exp(lpsiUN[3]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[3,14] <- 0
      #psiUN[3,15] <- 0

      psiUN[4,1] <- 0
      psiUN[4,2] <- exp(lpsiUN[4]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[4,3] <- 0
      psiUN[4,4] <- 0
      psiUN[4,5] <- exp(lpsiUN[4]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[4,6] <- 0
      psiUN[4,7] <- 0
      psiUN[4,8] <- exp(lpsiUN[4]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[4,9] <- 0
      psiUN[4,10] <- 0
      psiUN[4,11] <- exp(lpsiUN[4]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      psiUN[4,12] <- 0
      psiUN[4,13] <- 0
      psiUN[4,14] <- exp(lpsiUN[4]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      #psiUN[4,15] <- 0

    # Calculate the last transition probability
      psiDB[5,1] <- 1-psiDB[1,1]-psiDB[2,1]-psiDB[3,1]-psiDB[4,1]
      psiDB[5,2] <- 0
      psiDB[5,3] <- 0
      psiDB[5,4] <- 1-psiDB[1,4]-psiDB[2,4]-psiDB[3,4]-psiDB[4,4]
      psiDB[5,5] <- 0
      psiDB[5,6] <- 0
      psiDB[5,7] <- 1-psiDB[1,7]-psiDB[2,7]-psiDB[3,7]-psiDB[4,7]
      psiDB[5,8] <- 0
      psiDB[5,9] <- 0
      psiDB[5,10] <- 1-psiDB[1,10]-psiDB[2,10]-psiDB[3,10]-psiDB[4,10]
      psiDB[5,11] <- 0
      psiDB[5,12] <- 0
      psiDB[5,13] <- 1-psiDB[1,13]-psiDB[2,13]-psiDB[3,13]-psiDB[4,13]
      psiDB[5,14] <- 0
      #psiDB[5,15] <- 0
      
      psiJB[5,1] <- 0
      psiJB[5,2] <- 1-psiJB[1,2]-psiJB[2,2]-psiJB[3,2]-psiJB[4,2]
      psiJB[5,3] <- 0
      psiJB[5,4] <- 0
      psiJB[5,5] <- 1-psiJB[1,5]-psiJB[2,5]-psiJB[3,5]-psiJB[4,5]
      psiJB[5,6] <- 0
      psiJB[5,7] <- 0
      psiJB[5,8] <- 1-psiJB[1,8]-psiJB[2,8]-psiJB[3,8]-psiJB[4,8]
      psiJB[5,9] <- 0
      psiJB[5,10] <- 0
      psiJB[5,11] <- 1-psiJB[1,11]-psiJB[2,11]-psiJB[3,11]-psiJB[4,11]
      psiJB[5,12] <- 0
      psiJB[5,13] <- 0
      psiJB[5,14] <- 1-psiJB[1,14]-psiJB[2,14]-psiJB[3,14]-psiJB[4,14]
      #psiJB[5,15] <- 0

      psiMI[5,1] <- 0
      psiMI[5,2] <- 1-psiMI[1,2]-psiMI[2,2]-psiMI[3,2]-psiMI[4,2]
      psiMI[5,3] <- 0
      psiMI[5,4] <- 0
      psiMI[5,5] <- 1-psiMI[1,5]-psiMI[2,5]-psiMI[3,5]-psiMI[4,5]
      psiMI[5,6] <- 0
      psiMI[5,7] <- 0
      psiMI[5,8] <- 1-psiMI[1,8]-psiMI[2,8]-psiMI[3,8]-psiMI[4,8]
      psiMI[5,9] <- 0
      psiMI[5,10] <- 0
      psiMI[5,11] <- 1-psiMI[1,11]-psiMI[2,11]-psiMI[3,11]-psiMI[4,11]
      psiMI[5,12] <- 0
      psiMI[5,13] <- 0
      psiMI[5,14] <- 1-psiMI[1,14]-psiMI[2,14]-psiMI[3,14]-psiMI[4,14]
      #psiMI[5,15] <- 0

      psiAR[5,1] <- 0
      psiAR[5,2] <- 0
      psiAR[5,3] <- 1-psiAR[1,3]-psiAR[2,3]-psiAR[3,3]-psiAR[4,3]
      psiAR[5,4] <- 0
      psiAR[5,5] <- 0
      psiAR[5,6] <- 1-psiAR[1,6]-psiAR[2,6]-psiAR[3,6]-psiAR[4,6]
      psiAR[5,7] <- 0
      psiAR[5,8] <- 0
      psiAR[5,9] <- 1-psiAR[1,9]-psiAR[2,9]-psiAR[3,9]-psiAR[4,9]
      psiAR[5,10] <- 0
      psiAR[5,11] <- 0
      psiAR[5,12] <- 1-psiAR[1,12]-psiAR[2,12]-psiAR[3,12]-psiAR[4,12]
      psiAR[5,13] <- 0
      psiAR[5,14] <- 0
      #psiAR[5,15] <- 1-psiAR[1,15]-psiAR[2,15]-psiAR[3,15]-psiAR[4,15]

      for (t in 1:(n.occasions-1)){
      psiUN[5,t] <- 1-psiUN[1,t]-psiUN[2,t]-psiUN[3,t]-psiUN[4,t]
      }

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phiDB[t] * psiDB[2,t]
      ps[1,i,t,3] <- phiDB[t] * psiDB[3,t]
      ps[1,i,t,4] <- 0
      ps[1,i,t,5] <- phiDB[t] * psiDB[5,t]
      ps[1,i,t,6] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 0
      ps[2,i,t,4] <- phiJB[t] * psiJB[4,t]
      ps[2,i,t,5] <- phiJB[t] * psiJB[5,t]
      ps[2,i,t,6] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phiMI[t] * psiMI[4,t]
      ps[3,i,t,5] <- phiMI[t] * psiMI[5,t]
      ps[3,i,t,6] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[1,t]
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- phiAR[t] * psiAR[5,t]
      ps[4,i,t,6] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- phiUN[t] * psiUN[1,t]
      ps[5,i,t,2] <- phiUN[t] * psiUN[2,t]
      ps[5,i,t,3] <- phiUN[t] * psiUN[3,t]
      ps[5,i,t,4] <- phiUN[t] * psiUN[4,t]
      ps[5,i,t,5] <- phiUN[t] * psiUN[5,t]
      ps[5,i,t,6] <- 1-phiUN[t]     
      
      ps[6,i,t,1] <- 0
      ps[6,i,t,2] <- 0
      ps[6,i,t,3] <- 0
      ps[6,i,t,4] <- 0
      ps[6,i,t,5] <- 0
      ps[6,i,t,6] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pDB[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[1,i,t,4] <- 0
      po[1,i,t,5] <- 1-pDB[t]
      
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pJB[t]
      po[2,i,t,3] <- 0
      po[2,i,t,4] <- 0
      po[2,i,t,5] <- 1-pJB[t]
      
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- pMI[t]
      po[3,i,t,4] <- 0
      po[3,i,t,5] <- 1-pMI[t]
      
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 0
      po[4,i,t,4] <- pAR[t]
      po[4,i,t,5] <- 1-pAR[t]
      
      po[5,i,t,1] <- 0
      po[5,i,t,2] <- 0
      po[5,i,t,3] <- 0
      po[5,i,t,4] <- 0
      po[5,i,t,5] <- 1
      
      po[6,i,t,1] <- 0
      po[6,i,t,2] <- 0
      po[6,i,t,3] <- 0
      po[6,i,t,4] <- 0
      po[6,i,t,5] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture   
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 5))

# modify function for initial z values to follow annual cycle
mod.init.z2 <- function(ch, notseen){
  z.known <- ch # get known states
  z.known[z.known==notseen] <- NA
  z.init <- matrix(NA, ncol = dim(rCH)[1], nrow = n.occasions)
  z.init[1:dim(z.init)[1],] <- rep(c(1, 2, 4), n.occasions/3)
  z.init <- t(z.init)
  z.init[!is.na(z.known)] <- NA # replace known values with NA
  for (i in 1:dim(z.init)[1]){z.init[i,1:f[i]] <- NA}
  return(z.init)
}

# # Function to create initial values for psi if time-varying
# init.psi <- function(cols){
#   psi.array <- array(NA,dim=c(cols,n.occasions-1))
#   for (i in 1:cols){
#     for (t in 1:(n.occasions-1)){
#       psi.array[i,t] <- rnorm(1)
#     } #t
#   } #i
#   return(psi.array)
# }

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = rnorm(4), lpsiJB = rnorm(4), lpsiMI = rnorm(4), lpsiAR = rnorm(4), lpsiUN = rnorm(4),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 1), sigmaJB = runif(1, 0, 1), sigmaMI = runif(1, 0, 1), sigmaAR = runif(1, 0, 1),
                         z = mod.init.z2(rCH, 5))} 


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "pDB","mean.pDB", "sigma2DB",
                "pJB", "mean.pJB", "sigma2JB",
                "pMI", "mean.pMI", "sigma2MI",
                "pAR", "mean.pAR", "sigma2AR")

# MCMC settings
ni <- 100000
nt <- 5
nb <- 50000
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(rekn.ms, digits = 3)

rekn.ms.update <- autojags(jags.data, inits, parameters, "simulation_scenarios.jags",
                           n.chains = nc, iter.increment = 50000, n.burnin = nb,
                           n.thin = nt, save.all.iter = FALSE, Rhat.limit = 1.1,
                           max.iter = 500000, parallel = TRUE)

print(rekn.ms.update, digits = 3)

rekn.ms.summary <- rekn.ms$summary

saveRDS(rekn.ms$summary, file = paste0("./analysis-output/ms-simulation-summary", Sys.Date(), ".rds"))

saveRDS(rekn.ms$sims.list, file = paste0("./analysis-output/ms-simulation-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

theme_set(theme_bw())

# parameter identifiability checks
sims.list <- readRDS("analysis-output/ms-simulation-simslist2021-02-08.rds")
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
   
   png(filename = paste0("analysis-output/parameter-identifiability/2021-02-08/",
                         i,"-","check.png"),
       width=4, height=3, units="in", res=600)
   
   print(ggplot(sims.list, aes(sims.list[,i])) +
            geom_density() +
            geom_hline(yintercept = 1, linetype = "dashed") +
            xlab(i))
   
   dev.off()
   
}

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "AR", "UN"),
                  value = c(phiDB, phiJB, phiMI, phiAR, phiUN))

# resighting
p.sim <- tibble(site = c(rep("DB", n.occasions-1),
                         rep("JB", n.occasions-1),
                         rep("MI", n.occasions-1),
                         rep("AR", n.occasions-1)),
                occasion = rep(1:14, 4),
                value = c(pDB, pJB, pMI, pAR)) %>% 
   na_if(0)

# transition
psi.sim <- tibble(transition = c("DB-JB", "DB-MI", "JB-AR", "MI-AR", "AR-DB"),
                  value = c(0.4, 0.4, 0.7, 0.6, 0.8))

psiUN.sim <- tibble(transition = c("UN-DB", "UN-JB", "UN-MI", "UN-AR"),
                    value = c(0.8, 0.4, 0.4, 0.7))

# format survival
phi.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4, mean.phi.5) %>% 
   pivot_longer(cols = 1:5, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
   mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
   mutate(site = str_replace(site, "mean.phi.3", "MI")) %>%
   mutate(site = str_replace(site, "mean.phi.4", "AR")) %>%
   mutate(site = str_replace(site, "mean.phi.5", "UN"))

# plot survival
phi.plot <- ggplot() +
   geom_violin(phi.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phi.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

# format resighting
p.mod.DB <- sims.list %>% 
   as.data.frame() %>% 
   select(pDB.3, pDB.6, pDB.9, pDB.12) %>% 
   pivot_longer(cols = 1:4, names_to = "occasion", values_to = "estimate") %>% 
   mutate(occasion = str_sub(occasion, 5)) %>% 
   mutate(site = "DB")


p.mod.JB <- sims.list %>% 
   as.data.frame() %>% 
   select(pJB.4, pJB.7, pJB.10, pJB.13) %>% 
   pivot_longer(cols = 1:4, names_to = "occasion", values_to = "estimate") %>% 
   mutate(occasion = str_sub(occasion, 5)) %>% 
   mutate(site = "JB")

p.mod.MI <- sims.list %>% 
   as.data.frame() %>% 
   select(pMI.1, pMI.4, pMI.7, pMI.10, pMI.13) %>% 
   pivot_longer(cols = 1:5, names_to = "occasion", values_to = "estimate") %>% 
   mutate(occasion = str_sub(occasion, 5)) %>% 
   mutate(site = "MI")

p.mod.AR <- sims.list %>% 
   as.data.frame() %>% 
   select(pAR.2, pAR.5, pAR.8, pAR.11, pAR.14) %>% 
   pivot_longer(cols = 1:5, names_to = "occasion", values_to = "estimate") %>% 
   mutate(occasion = str_sub(occasion, 5)) %>% 
   mutate(site = "AR")

p.mod <- bind_rows(p.mod.DB, p.mod.JB, p.mod.MI, p.mod.AR)

# plot resighting
p.plot <- ggplot() +
   geom_violin(p.mod, mapping = aes(x = as.numeric(occasion), y = estimate,
                                    group = occasion, fill = occasion), alpha = 0.6) +
   geom_point(p.sim, mapping = aes(x = as.numeric(occasion), y = value),
              size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   scale_x_continuous(breaks = seq(1, 14, 1)) +
   #coord_flip() +
   facet_wrap(site ~ ., nrow = 4) +
   xlab("Occasion") +
   ylab("Resighting probability") +
   theme(legend.position = "none")

# format transitions
psi.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiDB.2, psiDB.3, psiJB.9, psiMI.9, psiAR.11) %>% 
   pivot_longer(cols = 1:5, names_to = "transition", values_to = "estimate") %>% 
   mutate(transition = str_replace(transition, "psiDB.2", "DB-JB")) %>%
   mutate(transition = str_replace(transition, "psiDB.3", "DB-MI")) %>%
   mutate(transition = str_replace(transition, "psiJB.9", "JB-AR")) %>%
   mutate(transition = str_replace(transition, "psiMI.9", "MI-AR")) %>%
   mutate(transition = str_replace(transition, "psiAR.11", "AR-DB"))

# plot transitions
psiDB.plot <- ggplot() +
   geom_violin(psi.mod %>% filter(transition %in% c("DB-JB", "DB-MI")),
               mapping = aes(x = transition, y = estimate, group = transition,
                                      fill = transition), alpha = 0.6) +
   geom_point(psi.sim %>% filter(transition %in% c("DB-JB", "DB-MI")),
                                 mapping = aes(x = transition, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Transition") +
   ylab("Transition probability") +
   theme(legend.position = "none")

psi.plot <- ggplot() +
   geom_violin(psi.mod %>% filter(!transition %in% c("DB-JB", "DB-MI")),
               mapping = aes(x = transition, y = estimate, group = transition,
                             fill = transition), alpha = 0.6) +
   geom_point(psi.sim %>% filter(!transition %in% c("DB-JB", "DB-MI")),
              mapping = aes(x = transition, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Transition") +
   ylab("Transition probability") +
   theme(legend.position = "none")

# plot all in grid
bottom.row <- plot_grid(phi.plot,psiDB.plot, psi.plot, labels = c("A", "B", "C"), ncol = 3)

png(filename = "figures/plot.sim.results.png", width = 8, height = 10,
    units = "in", res = 600)

plot_grid(bottom.row, p.plot, nrow = 2, labels = c("", "D"), rel_heights = c(2/5, 3/5))

dev.off()
