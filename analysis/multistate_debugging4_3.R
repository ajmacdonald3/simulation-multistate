################################################################################
# MULTISTATE SEASONAL SURVIVAL MODEL
# Gradual debugging
#
################################################################################

library(rjags)
library(jagsUI)
library(tidyverse)
library(cowplot)
library(viridis)

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

# function to create multistate m-array
marray <- function(CH, unobs = 0){ # unobs = number of unobservable states
  n.states <- max(CH) + unobs
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  out <- matrix(0, ncol = n.states*(n.occasions-1)+1, nrow = n.states*(n.occasions-1))
  
  # remove capture histories of individuals marked in last occasion
  get.first <- function(x) min(which(x!=0))
  first <- apply(CH, 1, get.first)
  last.only <- which(first==n.occasions)
  if (length(last.only) > 0) CH <- CH[-last.only,]
  
  # create m-array
  for (i in 1:nind){
    cap.occ <- which(CH[i,]!=0)
    state <- CH[i,cap.occ]
    if (length(state) == 1) {
      out[state[1] + n.states*(cap.occ[1]-1), n.states*(n.occasions-1)+1] <- 
        out[state[1] + n.states*(cap.occ[1]-1), n.states*(n.occasions-1)+1] + 1
    }
    
    if (length(state) > 1) {
      for (t in 2:length(cap.occ)){
        out[(cap.occ[t-1]-1)*n.states + state[t-1], (cap.occ[t]-2)*n.states + state[t]] <-
          out[(cap.occ[t-1]-1)*n.states + state[t-1], (cap.occ[t]-2)*n.states + state[t]] + 1
      } # t
      
      if (max(cap.occ) < n.occasions){
        out[(cap.occ[t]-1)*n.states + state[t], n.states*(n.occasions-1)+1] <-
          out[(cap.occ[t]-1)*n.states + state[t], n.states*(n.occasions-1)+1] + 1
      } # i
    }
  }
  
  return(out)
  
}

################################################################################

# Survival: constant
# Resighting: constant
# Transitions: constant with JB-JB set to 0

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 40
n.years <- 10
n.states <- 8
n.obs <- 8

# DB
phiDB1 <- 0.9

phiDB <- rep(c(phiDB1, 0, 0, 0), length.out = n.occasions-1)

pDB4 <- 0.65

pDB <- rep(c(0, 0, 0, pDB4), length.out = n.occasions-1)

# JB
phiJB2 <- 0.7

phiJB <- rep(c(0, phiJB2, 0, 0), length.out = n.occasions-1)

pJB1 <- 0.4

pJB <- rep(c(pJB1, 0, 0, 0), length.out = n.occasions-1)

# MI
phiMI2 <- 0.6

phiMI <- rep(c(0, phiMI2, 0, 0), length.out = n.occasions-1)

pMI1 <- 0.5

pMI <- rep(c(pMI1, 0, 0, 0), length.out = n.occasions-1)

# CC
phiCC2 <- 0.85
phiCC3 <- 0.75

phiCC <- rep(c(0, phiCC2, phiCC3, 0), length.out = n.occasions-1)

pCC1 <- 0.6
pCC2 <- 0.4

pCC <- rep(c(pCC1, pCC2, 0, 0), length.out = n.occasions-1)

# SE
phiSE1 <- 0.85
phiSE2 <- 0.9
phiSE3 <- 0.9
phiSE4 <- 0.95

phiSE <- rep(c(phiSE1, phiSE2, phiSE3, phiSE4), length.out = n.occasions-1)

pSE1 <- 0.6
pSE2 <- 0.6
pSE3 <- 0.4
pSE4 <- 0.3

pSE <- rep(c(pSE1, pSE2, pSE3, pSE4), length.out = n.occasions-1)

# BR
phiBR1 <- 0.7
phiBR2 <- 0.6
phiBR3 <- 0.85
phiBR4 <- 0.92

phiBR <- rep(c(phiBR1, phiBR2, phiBR3, phiBR4), length.out = n.occasions-1)

pBR1 <- 0.7
pBR2 <- 0.5
pBR3 <- 0.3
pBR4 <- 0.6

pBR <- rep(c(pBR1, pBR2, pBR3, pBR4), length.out = n.occasions-1)

# AR
phiAR1 <- 0.8
phiAR4 <- 0.93

phiAR <- rep(c(phiAR1, 0, 0, phiAR4), length.out = n.occasions-1)

pAR3 <- 0.6
pAR4 <- 0.4

pAR <- rep(c(0, 0, pAR3, pAR4), length.out = n.occasions-1)

# transition probabilities (sum to 1)
psiDB.DB <- rep(0, n.occasions)
psiDB.JB <- rep(c(0.35, 0, 0, 0), length.out = n.occasions)
psiDB.MI <- rep(c(0.25, 0, 0, 0), length.out = n.occasions)
psiDB.CC <- rep(c(0.20, 0, 0, 0), length.out = n.occasions)
psiDB.SE <- rep(c(0.10, 0, 0, 0), length.out = n.occasions)
psiDB.BR <- rep(c(0.10, 0, 0, 0), length.out = n.occasions)
psiDB.AR <- rep(0, n.occasions)

psiJB.DB <- rep(0, n.occasions)
psiJB.JB <- rep(0, n.occasions)
psiJB.MI <- rep(0, n.occasions)
psiJB.CC <- rep(c(0, 0.25, 0, 0), length.out = n.occasions)
psiJB.SE <- rep(c(0, 0.35, 0, 0), length.out = n.occasions)
psiJB.BR <- rep(c(0, 0.40, 0, 0), length.out = n.occasions)
psiJB.AR <- rep(0, n.occasions)

psiMI.DB <- rep(0, n.occasions)
psiMI.JB <- rep(0, n.occasions)
psiMI.MI <- rep(0, n.occasions)
psiMI.CC <- rep(c(0, 0.40, 0, 0), length.out = n.occasions)
psiMI.SE <- rep(c(0, 0.40, 0, 0), length.out = n.occasions)
psiMI.BR <- rep(c(0, 0.20, 0, 0), length.out = n.occasions)
psiMI.AR <- rep(0, n.occasions)

psiCC.DB <- rep(0, n.occasions)
psiCC.JB <- rep(0, n.occasions)
psiCC.MI <- rep(0, n.occasions)
psiCC.CC <- rep(c(0, 0.5, 0, 0), length.out = n.occasions) 
psiCC.SE <- rep(c(0, 0.2, 0.3, 0), length.out = n.occasions)
psiCC.BR <- rep(c(0, 0.1, 0.4, 0), length.out = n.occasions)
psiCC.AR <- rep(c(0, 0.2, 0.3, 0), length.out = n.occasions)

psiSE.DB <- rep(c(0, 0, 0, 0.6), length.out = n.occasions)
psiSE.JB <- rep(c(0.3, 0, 0, 0), length.out = n.occasions)
psiSE.MI <- rep(c(0.4, 0, 0, 0), length.out = n.occasions)
psiSE.CC <- rep(c(0.1, 0.2, 0, 0), length.out = n.occasions)
psiSE.SE <- rep(c(0.05, 0.3, 0.5, 0.1), length.out = n.occasions)
psiSE.BR <- rep(c(0.05, 0.5, 0.2, 0.1), length.out = n.occasions)
psiSE.AR <- rep(c(0, 0, 0.3, 0.2), length.out = n.occasions)

psiBR.DB <- rep(c(0, 0, 0, 0.3), length.out = n.occasions)
psiBR.JB <- rep(c(0.4, 0, 0, 0), length.out = n.occasions)
psiBR.MI <- rep(c(0.3, 0, 0, 0), length.out = n.occasions)
psiBR.CC <- rep(c(0.05, 0.2, 0, 0), length.out = n.occasions)
psiBR.SE <- rep(c(0.1, 0.05, 0.05, 0.2), length.out = n.occasions)
psiBR.BR <- rep(c(0.05, 0.75, 0.45, 0.4), length.out = n.occasions)
psiBR.AR <- rep(c(0, 0, 0.5, 0.1), length.out = n.occasions)

psiAR.DB <- rep(c(0, 0, 0, 0.7), length.out = n.occasions)
psiAR.JB <- rep(c(0.4, 0, 0, 0), length.out = n.occasions)
psiAR.MI <- rep(c(0.4, 0, 0, 0), length.out = n.occasions)
psiAR.CC <- rep(c(0.10, 0, 0, 0), length.out = n.occasions)
psiAR.SE <- rep(c(0.05, 0, 0, 0.1), length.out = n.occasions)
psiAR.BR <- rep(c(0.05, 0, 0, 0.1), length.out = n.occasions)
psiAR.AR <- rep(c(0, 0, 0, 0.1), length.out = n.occasions)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(100, 0, 0, 0), length.out = n.occasions) # DB
marked[,2] <- rep(c(0, 30, 0, 0), length.out = n.occasions) # JB
marked[,3] <- rep(c(0, 30, 0, 0), length.out = n.occasions) # MI
marked[,4] <- rep(c(0, 20, 20, 0), length.out = n.occasions) # CC
marked[,5] <- rep(c(10, 10, 10, 10), length.out = n.occasions) # SE
marked[,6] <- rep(c(5, 5, 5, 5), length.out = n.occasions) # BR
marked[,7] <- rep(c(20, 0, 0, 20), length.out = n.occasions) # AR
marked[,8] <- rep(0, n.occasions)

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
      phiDB[t]*psiDB.DB[t], phiDB[t]*psiDB.JB[t], phiDB[t]*psiDB.MI[t], phiDB[t]*psiDB.CC[t], phiDB[t]*psiDB.SE[t], phiDB[t]*psiDB.BR[t], phiDB[t]*psiDB.AR[t], 1-phiDB[t],
      phiJB[t]*psiJB.DB[t], phiJB[t]*psiJB.JB[t], phiJB[t]*psiJB.MI[t], phiJB[t]*psiJB.CC[t], phiJB[t]*psiJB.SE[t], phiJB[t]*psiJB.BR[t], phiJB[t]*psiJB.AR[t], 1-phiJB[t],
      phiMI[t]*psiMI.DB[t], phiMI[t]*psiMI.JB[t], phiMI[t]*psiMI.MI[t], phiMI[t]*psiMI.CC[t], phiMI[t]*psiMI.SE[t], phiMI[t]*psiMI.BR[t], phiMI[t]*psiMI.AR[t], 1-phiMI[t],
      phiCC[t]*psiCC.DB[t], phiCC[t]*psiCC.JB[t], phiCC[t]*psiCC.MI[t], phiCC[t]*psiCC.CC[t], phiCC[t]*psiCC.SE[t], phiCC[t]*psiCC.BR[t], phiCC[t]*psiCC.AR[t], 1-phiCC[t],
      phiSE[t]*psiSE.DB[t], phiSE[t]*psiSE.JB[t], phiSE[t]*psiSE.MI[t], phiSE[t]*psiSE.CC[t], phiSE[t]*psiSE.SE[t], phiSE[t]*psiSE.BR[t], phiSE[t]*psiSE.AR[t], 1-phiSE[t],
      phiBR[t]*psiBR.DB[t], phiBR[t]*psiBR.JB[t], phiBR[t]*psiBR.MI[t], phiBR[t]*psiBR.CC[t], phiBR[t]*psiBR.SE[t], phiBR[t]*psiBR.BR[t], phiBR[t]*psiBR.AR[t], 1-phiBR[t],
      phiAR[t]*psiAR.DB[t], phiAR[t]*psiAR.JB[t], phiAR[t]*psiAR.MI[t], phiAR[t]*psiAR.CC[t], phiAR[t]*psiAR.SE[t], phiAR[t]*psiAR.BR[t], phiAR[t]*psiAR.AR[t], 1-phiAR[t],
      0,              0,              0,              0,              0,              0,              0,              1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pDB[t], 0,   0,   0,   0,   0,   0,   1-pDB[t],
      0,   pJB[t], 0,   0,   0,   0,   0,   1-pJB[t],
      0,   0,   pMI[t], 0,   0,   0,   0,   1-pMI[t],
      0,   0,   0,   pCC[t], 0,   0,   0,   1-pCC[t],
      0,   0,   0,   0,   pSE[t], 0,   0,   1-pSE[t],
      0,   0,   0,   0,   0,   pBR[t], 0,   1-pBR[t],
      0,   0,   0,   0,   0,   0,   pAR[t], 1-pAR[t],
      0,   0,   0,   0,   0,   0,   0,   1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = NA)
CH <- sim$CH
CH.TRUE <- sim$CH.TRUE

marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

# JAGS model code
sink("simulation_scenarios.jags")
cat("
model {
  
  # Priors and constraints
  # Survival and recapture: uniform
  
  # DB
  for (t in 1:(n.occasions-1)){
    phiDB[t] <- mean.phiDB[season[t]]
    pDB[t] <- mean.pDB[season[t]]
  }
  
  for (s in 1:4){
    mean.phiDB[s] ~ dunif(0, 1)
    mean.pDB[s] ~ dunif(0, 1)
  }
  
  # JB
  for (t in 1:(n.occasions-1)){
    phiJB[t] <- mean.phiJB[season[t]]
    pJB[t] <- mean.pJB[season[t]]
  }
  
  for (s in 1:4){
    mean.phiJB[s] ~ dunif(0, 1)
    mean.pJB[s] ~ dunif(0, 1)
  }
  
  # MI
  for (t in 1:(n.occasions-1)){
    phiMI[t] <- mean.phiMI[season[t]]
    pMI[t] <- mean.pMI[season[t]]
  }
  
  for (s in 1:4){
    mean.phiMI[s] ~ dunif(0, 1)
    mean.pMI[s] ~ dunif(0, 1)
  }
  
  # CC
    for (t in 1:(n.occasions-1)){
    phiCC[t] <- mean.phiCC[season[t]]
    pCC[t] <- mean.pCC[season[t]]
  }
  
  for (s in 1:4){
    mean.phiCC[s] ~ dunif(0, 1)
    mean.pCC[s] ~ dunif(0, 1)
  }
  
  # SE
  for (t in 1:(n.occasions-1)){
    phiSE[t] <- mean.phiSE[season[t]]
    pSE[t] <- mean.pSE[season[t]]
  }
  
  for (s in 1:4){
    mean.phiSE[s] ~ dunif(0, 1)
    mean.pSE[s] ~ dunif(0, 1)
  }
  
  # BR
  for (t in 1:(n.occasions-1)){
    phiBR[t] <- mean.phiBR[season[t]]
    pBR[t] <- mean.pBR[season[t]]
  }
  
  for (s in 1:4){
    mean.phiBR[s] ~ dunif(0, 1)
    mean.pBR[s] ~ dunif(0, 1)
  }
  
  # AR
  for (t in 1:(n.occasions-1)){
    phiAR[t] <- mean.phiAR[season[t]]
    pAR[t] <- mean.pAR[season[t]]
  }
  
  for (s in 1:4){
    mean.phiAR[s] ~ dunif(0, 1)
    mean.pAR[s] ~ dunif(0, 1)
  }
  
  # Transitions: gamma priors
  
  for (i in 1:7){
    db[i] ~ dgamma(1, 1)
    jb[i] ~ dgamma(1, 1)
    mi[i] ~ dgamma(1, 1)
    cc[i] ~ dgamma(1, 1)
    se[i] ~ dgamma(1, 1)
    br[i] ~ dgamma(1, 1)
    ar[i] ~ dgamma(1, 1)
  }

  for (t in 1:(n.occasions-1)){
    psiDB[1,t] <- (db[1]/sum(db[]))
    psiDB[2,t] <- (db[2]/sum(db[])) * season[t]
    psiDB[3,t] <- (db[3]/sum(db[])) * season[t]
    psiDB[4,t] <- (db[4]/sum(db[])) * season[t]
    psiDB[5,t] <- (db[5]/sum(db[])) * season[t]
    psiDB[6,t] <- (db[6]/sum(db[])) * season[t]
    psiDB[7,t] <- (db[7]/sum(db[]))
  }
    
  for (t in 1:(n.occasions-1)){
    psiJB[1,t] <- (jb[1]/sum(jb[]))
    psiJB[2,t] <- (jb[2]/sum(jb[]))
    psiJB[3,t] <- (jb[3]/sum(jb[]))
    psiJB[4,t] <- (jb[4]/sum(jb[])) * season[t]
    psiJB[5,t] <- (jb[5]/sum(jb[])) * season[t]
    psiJB[6,t] <- (jb[6]/sum(jb[])) * season[t]
    psiJB[7,t] <- (jb[7]/sum(jb[]))
  }
  
  for (t in 1:(n.occasions-1)){
    psiMI[1,t] <- (mi[1]/sum(mi[]))
    psiMI[2,t] <- (mi[2]/sum(mi[]))
    psiMI[3,t] <- (mi[3]/sum(mi[]))
    psiMI[4,t] <- (mi[4]/sum(mi[])) * season[t]
    psiMI[5,t] <- (mi[5]/sum(mi[])) * season[t]
    psiMI[6,t] <- (mi[6]/sum(mi[])) * season[t]
    psiMI[7,t] <- (mi[7]/sum(mi[]))
  }
  
for (t in 1:(n.occasions-1)){
    psiCC[1,t] <- (cc[1]/sum(cc[]))
    psiCC[2,t] <- (cc[2]/sum(cc[]))
    psiCC[3,t] <- (cc[3]/sum(cc[]))
    psiCC[4,t] <- (cc[4]/sum(cc[])) * season[t]
    psiCC[5,t] <- (cc[5]/sum(cc[])) * season[t]
    psiCC[6,t] <- (cc[6]/sum(cc[])) * season[t]
    psiCC[7,t] <- (cc[7]/sum(cc[])) * season[t]
  }

  for (t in 1:(n.occasions-1)){
    psiSE[1,t] <- (se[1]/sum(se[])) * season[t]
    psiSE[2,t] <- (se[2]/sum(se[])) * season[t]
    psiSE[3,t] <- (se[3]/sum(se[])) * season[t]
    psiSE[4,t] <- (se[4]/sum(se[])) * season[t]
    psiSE[5,t] <- (se[5]/sum(se[])) * season[t]
    psiSE[6,t] <- (se[6]/sum(se[])) * season[t]
    psiSE[7,t] <- (se[7]/sum(se[])) * season[t]
  }

  for (t in 1:(n.occasions-1)){
    psiBR[1,t] <- (br[1]/sum(br[])) * season[t]
    psiBR[2,t] <- (br[2]/sum(br[])) * season[t]
    psiBR[3,t] <- (br[3]/sum(br[])) * season[t]
    psiBR[4,t] <- (br[4]/sum(br[])) * season[t]
    psiBR[5,t] <- (br[5]/sum(br[])) * season[t]
    psiBR[6,t] <- (br[6]/sum(br[])) * season[t]
    psiBR[7,t] <- (br[7]/sum(br[])) * season[t]
  }
  
  for (t in 1:(n.occasions-1)){
    psiAR[1,t] <- (ar[1]/sum(ar[])) * season[t]
    psiAR[2,t] <- (ar[2]/sum(ar[])) * season[t]
    psiAR[3,t] <- (ar[3]/sum(ar[])) * season[t]
    psiAR[4,t] <- (ar[4]/sum(ar[])) * season[t]
    psiAR[5,t] <- (ar[5]/sum(ar[])) * season[t]
    psiAR[6,t] <- (ar[6]/sum(ar[])) * season[t]
    psiAR[7,t] <- (ar[7]/sum(ar[])) * season[t]
  }
    
  # Define state-transition and reencounter probabilities - note no i index - is no longer individual 
  for (t in 1:(n.occasions-1)){
    psi[1,t,1] <- phiDB[t] * psiDB[1,t]
    psi[1,t,2] <- phiDB[t] * psiDB[2,t]
    psi[1,t,3] <- phiDB[t] * psiDB[3,t]
    psi[1,t,4] <- phiDB[t] * psiDB[4,t]
    psi[1,t,5] <- phiDB[t] * psiDB[5,t]
    psi[1,t,6] <- phiDB[t] * psiDB[6,t]
    psi[1,t,7] <- phiDB[t] * psiDB[7,t]
    
    psi[2,t,1] <- phiJB[t] * psiJB[1,t]
    psi[2,t,2] <- phiJB[t] * psiJB[2,t]
    psi[2,t,3] <- phiJB[t] * psiJB[3,t]
    psi[2,t,4] <- phiJB[t] * psiJB[4,t]
    psi[2,t,5] <- phiJB[t] * psiJB[5,t]
    psi[2,t,6] <- phiJB[t] * psiJB[6,t]
    psi[2,t,7] <- phiJB[t] * psiJB[7,t]
    
    psi[3,t,1] <- phiMI[t] * psiMI[1,t]
    psi[3,t,2] <- phiMI[t] * psiMI[2,t]
    psi[3,t,3] <- phiMI[t] * psiMI[3,t]
    psi[3,t,4] <- phiMI[t] * psiMI[4,t]
    psi[3,t,5] <- phiMI[t] * psiMI[5,t]
    psi[3,t,6] <- phiMI[t] * psiMI[6,t]
    psi[3,t,7] <- phiMI[t] * psiMI[7,t]
    
    psi[4,t,1] <- phiCC[t] * psiCC[1,t]
    psi[4,t,2] <- phiCC[t] * psiCC[2,t]
    psi[4,t,3] <- phiCC[t] * psiCC[3,t]
    psi[4,t,4] <- phiCC[t] * psiCC[4,t]
    psi[4,t,5] <- phiCC[t] * psiCC[5,t]
    psi[4,t,6] <- phiCC[t] * psiCC[6,t]
    psi[4,t,7] <- phiCC[t] * psiCC[7,t]
    
    psi[5,t,1] <- phiSE[t] * psiSE[1,t]
    psi[5,t,2] <- phiSE[t] * psiSE[2,t]
    psi[5,t,3] <- phiSE[t] * psiSE[3,t]
    psi[5,t,4] <- phiSE[t] * psiSE[4,t]
    psi[5,t,5] <- phiSE[t] * psiSE[5,t]
    psi[5,t,6] <- phiSE[t] * psiSE[6,t]
    psi[5,t,7] <- phiSE[t] * psiSE[7,t]
    
    psi[6,t,1] <- phiBR[t] * psiBR[1,t]
    psi[6,t,2] <- phiBR[t] * psiBR[2,t]
    psi[6,t,3] <- phiBR[t] * psiBR[3,t]
    psi[6,t,4] <- phiBR[t] * psiBR[4,t]
    psi[6,t,5] <- phiBR[t] * psiBR[5,t]
    psi[6,t,6] <- phiBR[t] * psiBR[6,t]
    psi[6,t,7] <- phiBR[t] * psiBR[7,t]
    
    psi[7,t,1] <- phiAR[t] * psiAR[1,t]
    psi[7,t,2] <- phiAR[t] * psiAR[2,t]
    psi[7,t,3] <- phiAR[t] * psiAR[3,t]
    psi[7,t,4] <- phiAR[t] * psiAR[4,t]
    psi[7,t,5] <- phiAR[t] * psiAR[5,t]
    psi[7,t,6] <- phiAR[t] * psiAR[6,t]
    psi[7,t,7] <- phiAR[t] * psiAR[7,t]
    
    po[1,t] <- pDB[t]
    po[2,t] <- pJB[t]
    po[3,t] <- pMI[t]
    po[4,t] <- pCC[t]
    po[5,t] <- pSE[t]
    po[6,t] <- pBR[t]
    po[7,t] <- pAR[t]
    
    
    # Calculate probability of non-encounter (dq) and reshape the array for the encounter probabilities - below here says the same     
    for (s in 1:ns){
      dp[s,t,s] <- po[s,t]
      dq[s,t,s] <- 1-po[s,t]
    } # s
    for (s in 1:(ns-1)){
      for (m in (s+1):ns){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } # s
    } # m
    for (s in 2:ns){
      for (m in 1:(s-1)){
        dp[s,t,m] <- 0
        dq[s,t,m] <- 0
      } # s
    } # m
  } # t
  
  # Define the multinomial likelihood
  for (t in 1:((n.occasions-1)*ns)){
    marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t,], rel[t])
  }
  
  # Define the cell probabilities of the multistate m-array   
  # Define matrix U: product of probabilities of state-transition and non-encounter (this is just done because there is no product function for matrix multiplication in JAGS)
  for (t in 1:(n.occasions-2)){
    U[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones
    for (j in (t+1):(n.occasions-1)){
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% psi[,t,] %*% dq[,t,]
    }
  }
  U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
  # Diagonal
  for (t in 1:(n.occasions-2)){
    pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% psi[,t,] %*% dp[,t,]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% psi[,j,] %*% dp[,j,]
    }
  }
  pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- psi[,n.occasions-1,] %*% dp[,n.occasions-1,]
  
  # Below main diagonal
  for (t in 2:(n.occasions-1)){
    for (j in 1:(t-1)){
      pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
    } #j
  } #t
  
  # Last column: probability of non-recapture
  for (t in 1:((n.occasions-1)*ns)){
    pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-1)*ns)])
  } #t
  
  # Assess model fit using Freeman-Tukey statistic
  # Compute fit statistics for observed data
  for (t in 1:((n.occasions-1)*ns)){
    for (j in 1:(n.occasions*ns-(ns-1))){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
    } #j
  } #t
  # Generate replicate data and compute fit stats from them
  for (t in 1:((n.occasions-1)*ns)){
    marr.new[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t, ], rel[t])
    for (j in 1:(n.occasions*ns-(ns-1))){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
    } #j
  } #t
  fit <- sum(E.org[,])
  fit.new <- sum(E.new[,])
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr = marr, n.occasions = ncol(CH), rel = rowSums(marr),
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns),
                  season = rep(c(1, 2, 3, 4), length.out = n.occasions-1))

# Initial values 
inits <- function(){list(mean.phiDB = runif(4, 0, 1), mean.phiJB = runif(4, 0, 1), mean.phiMI = runif(4, 0, 1),
                         mean.phiCC = runif(4, 0, 1), mean.phiSE = runif(4, 0, 1), mean.phiBR = runif(4, 0, 1), mean.phiAR = runif(4, 0, 1),
                         mean.pDB = runif(4, 0, 1), mean.pJB = runif(4, 0, 1), mean.pMI = runif(4, 0, 1),
                         mean.pCC = runif(4, 0, 1), mean.pSE = runif(4, 0, 1), mean.pBR = runif(4, 0, 1), mean.pAR = runif(4, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phiDB", "mean.phiJB", "mean.phiMI", "mean.phiCC", "mean.phiSE", "mean.phiBR", "mean.phiAR",
                "mean.pDB", "mean.pJB", "mean.pMI", "mean.pCC", "mean.pSE", "mean.pBR", "mean.pAR",
                "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
                "fit", "fit.new")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 2 days)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/debugging/multistate-debugging4-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/debugging/multistate-debugging4-simslist", Sys.Date(), ".rds"))

# parameter identifiability checks
summary <- readRDS("./analysis-output/debugging/multistate-debugging4-summary2021-11-04.rds")
sims.list <- readRDS("./analysis-output/debugging/multistate-debugging4-simslist2021-11-04.rds")

sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/debugging4/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          #geom_hline(yintercept = 1, linetype = "dashed") +
          xlab(i))
  
  dev.off()
  
}
rm(i)

# evaluation of fit
mean(sims.list$fit.new > sims.list$fit)

ppcheck <- ggplot(sims.list, aes(x = fit, y = fit.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

png(filename = "figures/debugging/debugging4/multistate-debugging4.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phiDBJBMICC.sim <- tibble(site = c("DB-PreB", "JB-PostB1", "MI-PostB1", "CC-PostB1", "CC-PostB2"),
                          value = c(phiDB1, phiJB2, phiMI2, phiCC2, phiCC3))

phiSE.sim <- tibble(site = c("SE-PreB", "SE-PostB1", "SE-PostB2", "SE-NonB"),
                    value = c(mean.phiSE1, mean.phiSE2, mean.phiSE3, mean.phiSE4))

phiBR.sim <- tibble(site = c("BR-PreB", "BR-PostB1", "BR-PostB2", "BR-NonB"),
                    value = c(mean.phiBR1, mean.phiBR2, mean.phiBR3, mean.phiBR4))

phiAR.sim <- tibble(site = c("AR-PreB", "AR-NonB"),
                    value = c(phiAR1, phiAR4))

# resighting
pDBJBMICC.sim <- tibble(site = c("DB-PreB", "JB-PostB1", "MI-PostB1", "CC-PostB1", "CC-PostB2"),
                        value = c(pDB4, pJB1, pMI1, pCC1, pCC2))

pSE.sim <- tibble(site = c("SE-PreB", "SE-PostB1", "SE-PostB2", "SE-NonB"),
                  value = c(mean.pSE1, mean.pSE2, mean.pSE3, mean.pSE4))

pBR.sim <- tibble(site = c("BR-PreB", "BR-PostB1", "BR-PostB2", "BR-NonB"),
                  value = c(mean.pBR1, mean.pBR2, mean.pBR3, mean.pBR4))

pAR.sim <- tibble(site = c("AR-PreB", "AR-NonB"),
                  value = c(pAR4, pAR3))

# transitions
psiDB.sim <- tibble(site = c("DB-JB", "DB-MI", "DB-CC", "DB-SE", "DB-BR"),
                    value = c(0.35, 0.25, 0.20, 0.10, 0.10))

psiJB.sim <- tibble(site = c("JB-CC", "JB-SE", "JB-BR"),
                    value = c(0.25, 0.35, 0.40))

psiMI.sim <- tibble(site = c("MI-CC", "MI-SE", "MI-BR"),
                    value = c(0.4, 0.4, 0.2))

psiCC.sim <- tibble(site = c("CC-CC2", "CC-SE2", "CC-SE3", "CC-BR2", "CC-BR3", "CC-AR2", "CC-AR3"),
                    value = c(0.5, 0.2, 0.3, 0.1, 0.4, 0.2, 0.3))

psiSE.sim <- tibble(site = c("SE-DB4", "SE-JB1", "SE-MI1", "SE-CC1", "SE-CC2", "SE-SE1", "SE-SE2",
                             "SE-SE3", "SE-SE4", "SE-BR1", "SE-BR2", "SE-BR3", "SE-BR4", "SE-AR3", "SE-AR4"),
                    value = c(0.6, 0.3, 0.4, 0.1, 0.2, 0.05, 0.3, 0.5, 0.1, 0.05, 0.5, 0.2, 0.1, 0.3, 0.2))

psiBR.sim <- tibble(site = c("BR-DB4", "BR-JB1", "BR-MI1", "BR-CC1", "BR-CC2", "BR-SE1", "BR-SE2",
                             "BR-SE3", "BR-SE4", "BR-BR1", "BR-BR2", "BR-BR3", "BR-BR4", "BR-AR3", "BR-AR4"),
                    value = c(0.3, 0.4, 0.3, 0.05, 0.2, 0.1, 0.05, 0.05, 0.2, 0.05, 0.75, 0.45, 0.4, 0.5, 0.1))

psiAR.sim <- tibble(site = c("AR-DB4", "AR-JB1", "AR-MI1", "AR-CC1", "AR-SE1", "AR-SE4", "AR-BR1", "AR-BR4", "AR-AR4"),
                    value = c(0.7, 0.4, 0.4, 0.1, 0.05, 0.1, 0.05, 0.1, 0.1))

# format survival
phiDBJBMICC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiDB.1, mean.phiJB.2, mean.phiMI.2, mean.phiCC.2, mean.phiCC.3) %>% 
  pivot_longer(cols = 1:5, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiDB.1", "DB-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiJB.2", "JB-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiMI.2", "MI-PostB1")) %>% 
  mutate(site = str_replace(site, "mean.phiCC.2", "CC-PostB1")) %>% 
  mutate(site = str_replace(site, "mean.phiCC.3", "CC-PostB2"))

phiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiSE.1, mean.phiSE.2, mean.phiSE.3, mean.phiSE.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiSE.1", "SE-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiSE.2", "SE-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiSE.3", "SE-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.phiSE.4", "SE-NonB"))

phiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiBR.1, mean.phiBR.2, mean.phiBR.3, mean.phiBR.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiBR.1", "BR-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiBR.2", "BR-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiBR.3", "BR-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.phiBR.4", "BR-NonB"))

phiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiAR.1, mean.phiAR.4) %>% 
  pivot_longer(cols = 1:2, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiAR.1", "AR-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiAR.4", "AR-NonB"))

# plot survival
phiDBJBMICC.plot <- ggplot() +
  geom_violin(phiDBJBMICC.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiDBJBMICC.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiDBJBMICC.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiDBJBMICC.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

phiSE.plot <- ggplot() +
  geom_violin(phiSE.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiSE.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiSE.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiSE.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

phiBR.plot <- ggplot() +
  geom_violin(phiBR.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiBR.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiBR.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiBR.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

phiAR.plot <- ggplot() +
  geom_violin(phiAR.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiAR.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiAR.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiAR.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

# format resighting
pDBJBMICC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pDB.4, mean.pJB.1, mean.pMI.1, mean.pCC.1, mean.pCC.2) %>% 
  pivot_longer(cols = 1:5, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pDB.4", "DB-PreB")) %>%
  mutate(site = str_replace(site, "mean.pJB.1", "JB-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pMI.1", "MI-PostB1")) %>% 
  mutate(site = str_replace(site, "mean.pCC.1", "CC-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pCC.2", "CC-PostB2"))

pSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pSE.1, mean.pSE.2, mean.pSE.3, mean.pSE.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pSE.1", "SE-PreB")) %>%
  mutate(site = str_replace(site, "mean.pSE.2", "SE-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pSE.3", "SE-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.pSE.4", "SE-NonB"))

pBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pBR.1, mean.pBR.2, mean.pBR.3, mean.pBR.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pBR.1", "BR-PreB")) %>%
  mutate(site = str_replace(site, "mean.pBR.2", "BR-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pBR.3", "BR-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.pBR.4", "BR-NonB"))

pAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pAR.3, mean.pAR.4) %>% 
  pivot_longer(cols = 1:2, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pAR.4", "AR-PreB")) %>%
  mutate(site = str_replace(site, "mean.pAR.3", "AR-NonB"))

# plot resighting
pDBJBMICC.plot <- ggplot() +
  geom_violin(pDBJBMICC.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pDBJBMICC.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pDBJBMICC.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pDBJBMICC.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

pSE.plot <- ggplot() +
  geom_violin(pSE.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pSE.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pSE.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pSE.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

pBR.plot <- ggplot() +
  geom_violin(pBR.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pBR.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pBR.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pBR.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

pAR.plot <- ggplot() +
  geom_violin(pAR.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pAR.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pAR.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pAR.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

# format transitions
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB.2, psiDB.3, psiDB.4, psiDB.5, psiDB.6) %>% 
  pivot_longer(cols = 1:5, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiDB.2", "DB-JB")) %>%
  mutate(site = str_replace(site, "psiDB.3", "DB-MI")) %>%
  mutate(site = str_replace(site, "psiDB.4", "DB-CC")) %>%
  mutate(site = str_replace(site, "psiDB.5", "DB-SE")) %>%
  mutate(site = str_replace(site, "psiDB.6", "DB-BR"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiJB.11, psiJB.12, psiJB.13) %>% 
  pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiJB.11", "JB-CC")) %>%
  mutate(site = str_replace(site, "psiJB.12", "JB-SE")) %>%
  mutate(site = str_replace(site, "psiJB.13", "JB-BR"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiMI.11, psiMI.12, psiMI.13) %>% 
  pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiMI.11", "MI-CC")) %>%
  mutate(site = str_replace(site, "psiMI.12", "MI-SE")) %>%
  mutate(site = str_replace(site, "psiMI.13", "MI-BR"))

psiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiCC.11, psiCC.12, psiCC.13, psiCC.19, psiCC.20, psiCC.21) %>% 
  pivot_longer(cols = 1:6, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiCC.11", "CC-CC2")) %>%
  mutate(site = str_replace(site, "psiCC.12", "CC-SE2")) %>%
  mutate(site = str_replace(site, "psiCC.13", "CC-BR2")) %>% 
  mutate(site = str_replace(site, "psiCC.19", "CC-SE3")) %>%
  mutate(site = str_replace(site, "psiCC.20", "CC-BR3")) %>%
  mutate(site = str_replace(site, "psiCC.21", "CC-AR3"))

psiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiSE.2, psiSE.3, psiSE.4, psiSE.5, psiSE.6, psiSE.11, psiSE.12, psiSE.13,
         psiSE.19, psiSE.20, psiSE.21, psiSE.22, psiSE.26, psiSE.27, psiSE.28) %>% 
  pivot_longer(cols = 1:15, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiSE.2$", "SE-JB1")) %>%
  mutate(site = str_replace(site, "psiSE.3", "SE-MI1")) %>%
  mutate(site = str_replace(site, "psiSE.4", "SE-CC1")) %>% 
  mutate(site = str_replace(site, "psiSE.5", "SE-SE1")) %>%
  mutate(site = str_replace(site, "psiSE.6", "SE-BR1")) %>% 
  mutate(site = str_replace(site, "psiSE.11", "SE-CC2")) %>%
  mutate(site = str_replace(site, "psiSE.12", "SE-SE2")) %>% 
  mutate(site = str_replace(site, "psiSE.13", "SE-BR2")) %>%
  mutate(site = str_replace(site, "psiSE.19", "SE-SE3")) %>% 
  mutate(site = str_replace(site, "psiSE.20", "SE-BR3")) %>% 
  mutate(site = str_replace(site, "psiSE.21", "SE-AR3")) %>% 
  mutate(site = str_replace(site, "psiSE.22", "SE-DB4")) %>% 
  mutate(site = str_replace(site, "psiSE.26", "SE-SE4")) %>% 
  mutate(site = str_replace(site, "psiSE.27", "SE-BR4")) %>% 
  mutate(site = str_replace(site, "psiSE.28", "SE-AR4")) 

psiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiBR.2, psiBR.3, psiBR.4, psiBR.5, psiBR.6, psiBR.11, psiBR.12, psiBR.13,
         psiBR.19, psiBR.20, psiBR.21, psiBR.22, psiBR.26, psiBR.27, psiBR.28) %>% 
  pivot_longer(cols = 1:15, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiBR.2$", "BR-JB1")) %>%
  mutate(site = str_replace(site, "psiBR.3", "BR-MI1")) %>%
  mutate(site = str_replace(site, "psiBR.4", "BR-CC1")) %>% 
  mutate(site = str_replace(site, "psiBR.5", "BR-SE1")) %>%
  mutate(site = str_replace(site, "psiBR.6", "BR-BR1")) %>% 
  mutate(site = str_replace(site, "psiBR.11", "BR-CC2")) %>%
  mutate(site = str_replace(site, "psiBR.12", "BR-SE2")) %>% 
  mutate(site = str_replace(site, "psiBR.13", "BR-BR2")) %>%
  mutate(site = str_replace(site, "psiBR.19", "BR-SE3")) %>% 
  mutate(site = str_replace(site, "psiBR.20", "BR-BR3")) %>% 
  mutate(site = str_replace(site, "psiBR.21", "BR-AR3")) %>% 
  mutate(site = str_replace(site, "psiBR.22", "BR-DB4")) %>% 
  mutate(site = str_replace(site, "psiBR.26", "BR-SE4")) %>% 
  mutate(site = str_replace(site, "psiBR.27", "BR-BR4")) %>% 
  mutate(site = str_replace(site, "psiBR.28", "BR-AR4")) 

psiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiAR.2, psiAR.3, psiAR.4, psiAR.5, psiAR.6, psiAR.22, psiAR.26, psiAR.27, psiAR.28) %>% 
  pivot_longer(cols = 1:9, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiAR.2$", "AR-JB1")) %>%
  mutate(site = str_replace(site, "psiAR.3", "AR-MI1")) %>% 
  mutate(site = str_replace(site, "psiAR.4", "AR-CC1")) %>%
  mutate(site = str_replace(site, "psiAR.5", "AR-SE1")) %>%
  mutate(site = str_replace(site, "psiAR.6", "AR-BR1")) %>%
  mutate(site = str_replace(site, "psiAR.22", "AR-DB4")) %>% 
  mutate(site = str_replace(site, "psiAR.26", "AR-SE4")) %>%
  mutate(site = str_replace(site, "psiAR.27", "AR-BR4")) %>%
  mutate(site = str_replace(site, "psiAR.28", "AR-AR4"))

# plot transitions
psiDB.plot <- ggplot() +
  geom_violin(psiDB.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
  geom_boxplot(psiDB.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiDB.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiDB.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiJB.plot <- ggplot() +
  geom_violin(psiJB.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
  geom_boxplot(psiJB.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiJB.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiJB.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiMI.plot <- ggplot() +
  geom_violin(psiMI.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
  geom_boxplot(psiMI.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiMI.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiMI.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiCC.plot <- ggplot() +
  geom_violin(psiCC.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
  geom_boxplot(psiCC.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiCC.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiCC.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiSE.plot <- ggplot() +
  geom_violin(psiSE.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
  geom_boxplot(psiSE.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiSE.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiSE.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiBR.plot <- ggplot() +
  geom_violin(psiBR.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
  geom_boxplot(psiBR.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiBR.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiBR.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiAR.plot <- ggplot() +
  geom_violin(psiAR.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
  geom_boxplot(psiAR.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiAR.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiAR.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/debugging/debugging4/phiDBJBMICC.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiDBJBMICC.plot)

dev.off()

png(filename = "figures/debugging/debugging4/phiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiSE.plot)

dev.off()

png(filename = "figures/debugging/debugging4/phiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiBR.plot)

dev.off()

png(filename = "figures/debugging/debugging4/phiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiAR.plot)

dev.off()

# resighting
png(filename = "figures/debugging/debugging4/pDBJBMICC.png", width = 8, height = 8,
    units = "in", res = 600)

print(pDBJBMICC.plot)

dev.off()

png(filename = "figures/debugging/debugging4/pSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(pSE.plot)

dev.off()

png(filename = "figures/debugging/debugging4/pBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(pBR.plot)

dev.off()

png(filename = "figures/debugging/debugging4/pAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(pAR.plot)

dev.off()

# transition
png(filename = "figures/debugging/debugging4/psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/debugging/debugging4/psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/debugging/debugging4/psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

png(filename = "figures/debugging/debugging4/psiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiCC.plot)

dev.off()

png(filename = "figures/debugging/debugging4/psiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiSE.plot)

dev.off()

png(filename = "figures/debugging/debugging4/psiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiBR.plot)

dev.off()

png(filename = "figures/debugging/debugging4/psiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiAR.plot)

dev.off()

