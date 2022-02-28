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
mean.phiSE1 <- 0.85
mean.phiSE2 <- 0.9
mean.phiSE3 <- 0.9
mean.phiSE4 <- 0.95

phiSE <- rep(c(mean.phiSE1, mean.phiSE2, mean.phiSE3, mean.phiSE4), length.out = n.occasions-1)

mean.pSE1 <- 0.6
mean.pSE2 <- 0.6
mean.pSE3 <- 0.4
mean.pSE4 <- 0.3

pSE <- rep(c(mean.pSE1, mean.pSE2, mean.pSE3, mean.pSE4), length.out = n.occasions-1)

# BR
mean.phiBR1 <- 0.7
mean.phiBR2 <- 0.6
mean.phiBR3 <- 0.85
mean.phiBR4 <- 0.92

phiBR <- rep(c(mean.phiBR1, mean.phiBR2, mean.phiBR3, mean.phiBR4), length.out = n.occasions-1)

mean.pBR1 <- 0.7
mean.pBR2 <- 0.5
mean.pBR3 <- 0.3
mean.pBR4 <- 0.6

pBR <- rep(c(mean.pBR1, mean.pBR2, mean.pBR3, mean.pBR4), length.out = n.occasions-1)

# AR
phiAR1 <- 0.8
phiAR4 <- 0.93

phiAR <- rep(c(phiAR1, 0, 0, phiAR4), length.out = n.occasions-1)

pAR3 <- 0.6
pAR4 <- 0.4

pAR <- rep(c(0, 0, pAR3, pAR4), length.out = n.occasions-1)

# transition probabilities (sum to 1)
psiDB.DB <- rep(0, n.occasions-1)
psiDB.JB <- rep(c(0.35, 0, 0, 0), length.out = n.occasions-1)
psiDB.MI <- rep(c(0.25, 0, 0, 0), length.out = n.occasions-1)
psiDB.CC <- rep(c(0.20, 0, 0, 0), length.out = n.occasions-1)
psiDB.SE <- rep(c(0.10, 0, 0, 0), length.out = n.occasions-1)
psiDB.BR <- rep(c(0.10, 0, 0, 0), length.out = n.occasions-1)
psiDB.AR <- rep(0, n.occasions-1)

psiJB.DB <- rep(0, n.occasions-1)
psiJB.JB <- rep(0, n.occasions-1)
psiJB.MI <- rep(0, n.occasions-1)
psiJB.CC <- rep(c(0, 0.25, 0, 0), length.out = n.occasions-1)
psiJB.SE <- rep(c(0, 0.35, 0, 0), length.out = n.occasions-1)
psiJB.BR <- rep(c(0, 0.40, 0, 0), length.out = n.occasions-1)
psiJB.AR <- rep(0, n.occasions-1)

psiMI.DB <- rep(0, n.occasions-1)
psiMI.JB <- rep(0, n.occasions-1)
psiMI.MI <- rep(0, n.occasions-1)
psiMI.CC <- rep(c(0, 0.40, 0, 0), length.out = n.occasions-1)
psiMI.SE <- rep(c(0, 0.40, 0, 0), length.out = n.occasions-1)
psiMI.BR <- rep(c(0, 0.20, 0, 0), length.out = n.occasions-1)
psiMI.AR <- rep(0, n.occasions-1)

psiCC.DB <- rep(0, n.occasions-1)
psiCC.JB <- rep(0, n.occasions-1)
psiCC.MI <- rep(0, n.occasions-1)
psiCC.CC <- rep(c(0, 0.5, 0, 0), length.out = n.occasions-1) 
psiCC.SE <- rep(c(0, 0.2, 0.3, 0), length.out = n.occasions-1)
psiCC.BR <- rep(c(0, 0.3, 0.4, 0), length.out = n.occasions-1)
psiCC.AR <- rep(c(0, 0, 0.3, 0), length.out = n.occasions-1)

psiSE.DB <- rep(c(0, 0, 0, 0.6), length.out = n.occasions-1)
psiSE.JB <- rep(c(0.3, 0, 0, 0), length.out = n.occasions-1)
psiSE.MI <- rep(c(0.4, 0, 0, 0), length.out = n.occasions-1)
psiSE.CC <- rep(c(0.1, 0.2, 0, 0), length.out = n.occasions-1)
psiSE.SE <- rep(c(0.05, 0.3, 0.5, 0.1), length.out = n.occasions-1)
psiSE.BR <- rep(c(0.05, 0.5, 0.2, 0.1), length.out = n.occasions-1)
psiSE.AR <- rep(c(0, 0, 0.3, 0.2), length.out = n.occasions-1)

psiBR.DB <- rep(c(0, 0, 0, 0.3), length.out = n.occasions-1)
psiBR.JB <- rep(c(0.4, 0, 0, 0), length.out = n.occasions-1)
psiBR.MI <- rep(c(0.3, 0, 0, 0), length.out = n.occasions-1)
psiBR.CC <- rep(c(0.05, 0.2, 0, 0), length.out = n.occasions-1)
psiBR.SE <- rep(c(0.1, 0.05, 0.05, 0.2), length.out = n.occasions-1)
psiBR.BR <- rep(c(0.05, 0.75, 0.45, 0.4), length.out = n.occasions-1)
psiBR.AR <- rep(c(0, 0, 0.5, 0.1), length.out = n.occasions-1)

psiAR.DB <- rep(c(0, 0, 0, 0.7), length.out = n.occasions-1)
psiAR.JB <- rep(c(0.4, 0, 0, 0), length.out = n.occasions-1)
psiAR.MI <- rep(c(0.4, 0, 0, 0), length.out = n.occasions-1)
psiAR.CC <- rep(c(0.10, 0, 0, 0), length.out = n.occasions-1)
psiAR.SE <- rep(c(0.05, 0, 0, 0.1), length.out = n.occasions-1)
psiAR.BR <- rep(c(0.05, 0, 0, 0.1), length.out = n.occasions-1)
psiAR.AR <- rep(c(0, 0, 0, 0.1), length.out = n.occasions-1)

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
  
  for(s in 1:4){
  db[1,s] <- 0
  for (i in 2:6){
    db[i,s] ~ dgamma(1, 1)
  }
  db[7,s] <- 0
  
  psiDB1[s] <- (db[1,s]/sum(db[,s])) 
  psiDB2[s] <- (db[2,s]/sum(db[,s]))
  psiDB3[s] <- (db[3,s]/sum(db[,s]))
  psiDB4[s] <- (db[4,s]/sum(db[,s])) 
  psiDB5[s] <- (db[5,s]/sum(db[,s]))
  psiDB6[s] <- (db[6,s]/sum(db[,s]))
  psiDB7[s] <- (db[7,s]/sum(db[,s])) 

  }
  
  for (t in 1:(n.occasions-1)){
    psiDB[1,t] <- psiDB1[season[t]]
    psiDB[2,t] <- psiDB2[season[t]]
    psiDB[3,t] <- psiDB3[season[t]]
    psiDB[4,t] <- psiDB4[season[t]]
    psiDB[5,t] <- psiDB5[season[t]]
    psiDB[6,t] <- psiDB6[season[t]]
    psiDB[7,t] <- psiDB7[season[t]]
  }
  
  for(s in 1:4){
  for (i in 1:3){
    jb[i,s] <- 0
  }
  
  for (i in 4:6){
    jb[i,s] ~ dgamma(1, 1)
  }
  
    jb[7,s] <- 0
  
  psiJB1[s] <- (jb[1,s]/sum(jb[,s])) 
  psiJB2[s] <- (jb[2,s]/sum(jb[,s]))
  psiJB3[s] <- (jb[3,s]/sum(jb[,s]))
  psiJB4[s] <- (jb[4,s]/sum(jb[,s])) 
  psiJB5[s] <- (jb[5,s]/sum(jb[,s]))
  psiJB6[s] <- (jb[6,s]/sum(jb[,s]))
  psiJB7[s] <- (jb[7,s]/sum(jb[,s])) 
    
  }  
    
  for (t in 1:(n.occasions-1)){
    psiJB[1,t] <- psiJB1[season[t]]
    psiJB[2,t] <- psiJB2[season[t]]
    psiJB[3,t] <- psiJB3[season[t]]
    psiJB[4,t] <- psiJB4[season[t]]
    psiJB[5,t] <- psiJB5[season[t]]
    psiJB[6,t] <- psiJB6[season[t]]
    psiJB[7,t] <- psiJB7[season[t]]
  }
  
  for (s in 1:4){
  for (i in 1:3){
    mi[i,s] <- 0
  }
  
  for (i in 4:6){
    mi[i,s] ~ dgamma(1, 1)
  }
  
    mi[7,s] <- 0
  
  psiMI1[s] <- (mi[1,s]/sum(mi[,s])) 
  psiMI2[s] <- (mi[2,s]/sum(mi[,s]))
  psiMI3[s] <- (mi[3,s]/sum(mi[,s]))
  psiMI4[s] <- (mi[4,s]/sum(mi[,s])) 
  psiMI5[s] <- (mi[5,s]/sum(mi[,s]))
  psiMI6[s] <- (mi[6,s]/sum(mi[,s]))
  psiMI7[s] <- (mi[7,s]/sum(mi[,s]))
  
  }
  
  for (t in 1:(n.occasions-1)){
    psiMI[1,t] <- psiMI1[season[t]]
    psiMI[2,t] <- psiMI2[season[t]]
    psiMI[3,t] <- psiMI3[season[t]]
    psiMI[4,t] <- psiMI4[season[t]]
    psiMI[5,t] <- psiMI5[season[t]]
    psiMI[6,t] <- psiMI6[season[t]]
    psiMI[7,t] <- psiMI7[season[t]]
  }
  
  for (s in 1:4){
  for (i in 1:3){
    cc[i,s] <- 0
  }
  
  for (i in 4:7){
    cc[i,s] ~ dgamma(1, 1)
  }
  
  psiCC1[s] <- (cc[1,s]/sum(cc[,s])) 
  psiCC2[s] <- (cc[2,s]/sum(cc[,s]))
  psiCC3[s] <- (cc[3,s]/sum(cc[,s]))
  psiCC4[s] <- (cc[4,s]/sum(cc[,s])) 
  psiCC5[s] <- (cc[5,s]/sum(cc[,s]))
  psiCC6[s] <- (cc[6,s]/sum(cc[,s]))
  psiCC7[s] <- (cc[7,s]/sum(cc[,s]))
  
  }
  
for (t in 1:(n.occasions-1)){
    psiCC[1,t] <- psiCC1[season[t]]
    psiCC[2,t] <- psiCC2[season[t]]
    psiCC[3,t] <- psiCC3[season[t]]
    psiCC[4,t] <- psiCC4[season[t]]
    psiCC[5,t] <- psiCC5[season[t]]
    psiCC[6,t] <- psiCC6[season[t]]
    psiCC[7,t] <- psiCC7[season[t]]
  }
  
  for (s in 1:4){
  for (i in 1:7){
    se[i,s] ~ dgamma(1, 1)
  }
  
  psiSE1[s] <- (se[1,s]/sum(se[,s])) 
  psiSE2[s] <- (se[2,s]/sum(se[,s]))
  psiSE3[s] <- (se[3,s]/sum(se[,s]))
  psiSE4[s] <- (se[4,s]/sum(se[,s])) 
  psiSE5[s] <- (se[5,s]/sum(se[,s]))
  psiSE6[s] <- (se[6,s]/sum(se[,s]))
  psiSE7[s] <- (se[7,s]/sum(se[,s]))
  
  }

  for (t in 1:(n.occasions-1)){
    psiSE[1,t] <- psiSE1[season[t]]
    psiSE[2,t] <- psiSE2[season[t]]
    psiSE[3,t] <- psiSE3[season[t]]
    psiSE[4,t] <- psiSE4[season[t]]
    psiSE[5,t] <- psiSE5[season[t]]
    psiSE[6,t] <- psiSE6[season[t]]
    psiSE[7,t] <- psiSE7[season[t]]
  }
  
  for (s in 1:4){
  for (i in 1:7){
    br[i,s] ~ dgamma(1, 1)
  }
  
  psiBR1[s] <- (br[1,s]/sum(br[,s])) 
  psiBR2[s] <- (br[2,s]/sum(br[,s]))
  psiBR3[s] <- (br[3,s]/sum(br[,s]))
  psiBR4[s] <- (br[4,s]/sum(br[,s])) 
  psiBR5[s] <- (br[5,s]/sum(br[,s]))
  psiBR6[s] <- (br[6,s]/sum(br[,s]))
  psiBR7[s] <- (br[7,s]/sum(br[,s]))
  
  }

  for (t in 1:(n.occasions-1)){
    psiBR[1,t] <- psiBR1[season[t]]
    psiBR[2,t] <- psiBR2[season[t]]
    psiBR[3,t] <- psiBR3[season[t]]
    psiBR[4,t] <- psiBR4[season[t]]
    psiBR[5,t] <- psiBR5[season[t]]
    psiBR[6,t] <- psiBR6[season[t]]
    psiBR[7,t] <- psiBR7[season[t]]
  }
  
  for (s in 1:4){
  for (i in 1:7){
    ar[i,s] ~ dgamma(1, 1)
  }
  
  psiAR1[s] <- (ar[1,s]/sum(ar[,s])) 
  psiAR2[s] <- (ar[2,s]/sum(ar[,s]))
  psiAR3[s] <- (ar[3,s]/sum(ar[,s]))
  psiAR4[s] <- (ar[4,s]/sum(ar[,s])) 
  psiAR5[s] <- (ar[5,s]/sum(ar[,s]))
  psiAR6[s] <- (ar[6,s]/sum(ar[,s]))
  psiAR7[s] <- (ar[7,s]/sum(ar[,s]))
  
  }
  
  for (t in 1:(n.occasions-1)){
    psiAR[1,t] <- psiAR1[season[t]]
    psiAR[2,t] <- psiAR2[season[t]]
    psiAR[3,t] <- psiAR3[season[t]]
    psiAR[4,t] <- psiAR4[season[t]]
    psiAR[5,t] <- psiAR5[season[t]]
    psiAR[6,t] <- psiAR6[season[t]]
    psiAR[7,t] <- psiAR7[season[t]]
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
                "psiDB1", "psiDB2", "psiDB3", "psiDB4", "psiDB5", "psiDB6", "psiDB7",
                "psiJB1", "psiJB2", "psiJB3", "psiJB4", "psiJB5", "psiJB6", "psiJB7",
                "psiMI1", "psiMI2", "psiMI3", "psiMI4", "psiMI5", "psiMI6", "psiMI7",
                "psiCC1", "psiCC2", "psiCC3", "psiCC4", "psiCC5", "psiCC6", "psiCC7",
                "psiSE1", "psiSE2", "psiSE3", "psiSE4", "psiSE5", "psiSE6", "psiSE7",
                "psiBR1", "psiBR2", "psiBR3", "psiBR4", "psiBR5", "psiBR6", "psiBR7",
                "psiAR1", "psiAR2", "psiAR3", "psiAR4", "psiAR5", "psiAR6", "psiAR7",
                "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 2 days)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/debugging/multistate-debugging7-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/debugging/multistate-debugging7-simslist", Sys.Date(), ".rds"))

sim.summ <- sim.marr$summary

# parameter identifiability checks

sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/debugging7/",
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

png(filename = "figures/debugging/debugging7/ppcheck.png", width = 8, height = 6,
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
psiDB.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 28),
                    transition = c(rep("DB-DB", 4), rep("DB-JB", 4), rep("DB-MI", 4), rep("DB-CC", 4),
                                   rep("DB-SE", 4), rep("DB-BR", 4), rep("DB-AR", 4)),
                    value = c(psiDB.DB[1:4], psiDB.JB[1:4], psiDB.MI[1:4], psiDB.CC[1:4],
                              psiDB.SE[1:4], psiDB.BR[1:4], psiDB.AR[1:4]))

psiJB.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 28),
                    transition = c(rep("JB-DB", 4), rep("JB-JB", 4), rep("JB-MI", 4), rep("JB-CC", 4),
                                   rep("JB-SE", 4), rep("JB-BR", 4), rep("JB-AR", 4)),
                    value = c(psiJB.DB[1:4], psiJB.JB[1:4], psiJB.MI[1:4], psiJB.CC[1:4],
                              psiJB.SE[1:4], psiJB.BR[1:4], psiJB.AR[1:4]))

psiMI.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 28),
                    transition = c(rep("MI-DB", 4), rep("MI-JB", 4), rep("MI-MI", 4), rep("MI-CC", 4),
                                   rep("MI-SE", 4), rep("MI-BR", 4), rep("MI-AR", 4)),
                    value = c(psiMI.DB[1:4], psiMI.JB[1:4], psiMI.MI[1:4], psiMI.CC[1:4],
                              psiMI.SE[1:4], psiMI.BR[1:4], psiMI.AR[1:4]))

psiCC.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 28),
                    transition = c(rep("CC-DB", 4), rep("CC-JB", 4), rep("CC-MI", 4), rep("CC-CC", 4),
                                   rep("CC-SE", 4), rep("CC-BR", 4), rep("CC-AR", 4)),
                    value = c(psiCC.DB[1:4], psiCC.JB[1:4], psiCC.MI[1:4], psiCC.CC[1:4],
                              psiCC.SE[1:4], psiCC.BR[1:4], psiCC.AR[1:4]))

psiSE.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 28),
                    transition = c(rep("SE-DB", 4), rep("SE-JB", 4), rep("SE-MI", 4), rep("SE-CC", 4),
                                   rep("SE-SE", 4), rep("SE-BR", 4), rep("SE-AR", 4)),
                    value = c(psiSE.DB[1:4], psiSE.JB[1:4], psiSE.MI[1:4], psiSE.CC[1:4],
                              psiSE.SE[1:4], psiSE.BR[1:4], psiSE.AR[1:4]))

psiBR.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 28),
                    transition = c(rep("BR-DB", 4), rep("BR-JB", 4), rep("BR-MI", 4), rep("BR-CC", 4),
                                   rep("BR-SE", 4), rep("BR-BR", 4), rep("BR-AR", 4)),
                    value = c(psiBR.DB[1:4], psiBR.JB[1:4], psiBR.MI[1:4], psiBR.CC[1:4],
                              psiBR.SE[1:4], psiBR.BR[1:4], psiBR.AR[1:4]))

psiAR.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 28),
                    transition = c(rep("AR-DB", 4), rep("AR-JB", 4), rep("AR-MI", 4), rep("AR-CC", 4),
                                   rep("AR-SE", 4), rep("AR-BR", 4), rep("AR-AR", 4)),
                    value = c(psiAR.DB[1:4], psiAR.JB[1:4], psiAR.MI[1:4], psiAR.CC[1:4],
                              psiAR.SE[1:4], psiAR.BR[1:4], psiAR.AR[1:4]))

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
  ylim(0, 1) +
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
  ylim(0, 1) +
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
  ylim(0, 1) +
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
  ylim(0, 1) +
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
  ylim(0, 1) +
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
  ylim(0, 1) +
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
  ylim(0, 1) +
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
  ylim(0, 1) +
  xlab("Site") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

# format transitions
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB2.1, psiDB3.1, psiDB4.1, psiDB5.1, psiDB6.1) %>% 
  pivot_longer(cols = 1:5, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "psiDB2.1" ~ "1",
                            transition == "psiDB3.1" ~ "1",
                            transition == "psiDB4.1" ~ "1",
                            transition == "psiDB5.1" ~ "1",
                            transition == "psiDB6.1" ~ "1")) %>% 
  mutate(transition = str_replace(transition, "psiDB2.1", "DB-JB")) %>%
  mutate(transition = str_replace(transition, "psiDB3.1", "DB-MI")) %>%
  mutate(transition = str_replace(transition, "psiDB4.1", "DB-CC")) %>%
  mutate(transition = str_replace(transition, "psiDB5.1", "DB-SE")) %>%
  mutate(transition = str_replace(transition, "psiDB6.1", "DB-BR"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiJB4.2, psiJB5.2, psiJB6.2) %>% 
  pivot_longer(cols = 1:3, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "psiJB4.2" ~ "2",
                            transition == "psiJB5.2" ~ "2",
                            transition == "psiJB6.2" ~ "2")) %>% 
  mutate(transition = str_replace(transition, "psiJB4.2", "JB-CC")) %>%
  mutate(transition = str_replace(transition, "psiJB5.2", "JB-SE")) %>%
  mutate(transition = str_replace(transition, "psiJB6.2", "JB-BR"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiMI4.2, psiMI5.2, psiMI6.2) %>% 
  pivot_longer(cols = 1:3, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "psiMI4.2" ~ "2",
                            transition == "psiMI5.2" ~ "2",
                            transition == "psiMI6.2" ~ "2")) %>% 
  mutate(transition = str_replace(transition, "psiMI4.2", "MI-CC")) %>%
  mutate(transition = str_replace(transition, "psiMI5.2", "MI-SE")) %>%
  mutate(transition = str_replace(transition, "psiMI6.2", "MI-BR"))

psiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiCC4.2, psiCC5.2, psiCC5.3, psiCC6.2, psiCC6.3, psiCC7.3) %>% 
  pivot_longer(cols = 1:6, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "psiCC4.2" ~ "2",
                            transition == "psiCC5.2" ~ "2",
                            transition == "psiCC5.3" ~ "3",
                            transition == "psiCC6.2" ~ "2",
                            transition == "psiCC6.3" ~ "3",
                            transition == "psiCC7.3" ~ "3")) %>% 
  mutate(transition = str_replace(transition, "psiCC4.2", "CC-CC")) %>%
  mutate(transition = str_replace(transition, "psiCC5.2", "CC-SE")) %>%
  mutate(transition = str_replace(transition, "psiCC5.3", "CC-SE")) %>% 
  mutate(transition = str_replace(transition, "psiCC6.2", "CC-BR")) %>%
  mutate(transition = str_replace(transition, "psiCC6.3", "CC-BR")) %>% 
  mutate(transition = str_replace(transition, "psiCC7.3", "CC-AR"))

psiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiSE2.1, psiSE3.1, psiSE4.1, psiSE4.2, psiSE5.1, psiSE5.2, psiSE5.3,
         psiSE5.4, psiSE6.1, psiSE6.2, psiSE6.3, psiSE6.4, psiSE7.3, psiSE7.4) %>% 
  pivot_longer(cols = 1:14, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "psiSE2.1" ~ "1",
                            transition == "psiSE3.1" ~ "1",
                            transition == "psiSE4.1" ~ "1",
                            transition == "psiSE4.2" ~ "2",
                            transition == "psiSE5.1" ~ "1",
                            transition == "psiSE5.2" ~ "2",
                            transition == "psiSE5.3" ~ "3",
                            transition == "psiSE5.4" ~ "4",
                            transition == "psiSE6.1" ~ "1",
                            transition == "psiSE6.2" ~ "2",
                            transition == "psiSE6.3" ~ "3",
                            transition == "psiSE6.4" ~ "4",
                            transition == "psiSE7.3" ~ "3",
                            transition == "psiSE7.4" ~ "4")) %>% 
  mutate(transition = str_replace(transition, "psiSE2.1", "SE-JB")) %>%
  mutate(transition = str_replace(transition, "psiSE3.1", "SE-MI")) %>%
  mutate(transition = str_replace(transition, "psiSE4.1", "SE-CC")) %>%
  mutate(transition = str_replace(transition, "psiSE4.2", "SE-CC")) %>%
  mutate(transition = str_replace(transition, "psiSE5.1", "SE-SE")) %>%
  mutate(transition = str_replace(transition, "psiSE5.2", "SE-SE")) %>%
  mutate(transition = str_replace(transition, "psiSE5.3", "SE-SE")) %>% 
  mutate(transition = str_replace(transition, "psiSE5.4", "SE-SE")) %>% 
  mutate(transition = str_replace(transition, "psiSE6.1", "SE-BR")) %>%
  mutate(transition = str_replace(transition, "psiSE6.2", "SE-BR")) %>%
  mutate(transition = str_replace(transition, "psiSE6.3", "SE-BR")) %>% 
  mutate(transition = str_replace(transition, "psiSE6.4", "SE-BR")) %>% 
  mutate(transition = str_replace(transition, "psiSE7.3", "SE-AR")) %>% 
  mutate(transition = str_replace(transition, "psiSE7.4", "SE-AR"))

psiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiBR2.1, psiBR3.1, psiBR4.1, psiBR4.2, psiBR5.1, psiBR5.2, psiBR5.3,
         psiBR5.4, psiBR6.1, psiBR6.2, psiBR6.3, psiBR6.4, psiBR7.3, psiBR7.4) %>% 
  pivot_longer(cols = 1:14, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "psiBR2.1" ~ "1",
                            transition == "psiBR3.1" ~ "1",
                            transition == "psiBR4.1" ~ "1",
                            transition == "psiBR4.2" ~ "2",
                            transition == "psiBR5.1" ~ "1",
                            transition == "psiBR5.2" ~ "2",
                            transition == "psiBR5.3" ~ "3",
                            transition == "psiBR5.4" ~ "4",
                            transition == "psiBR6.1" ~ "1",
                            transition == "psiBR6.2" ~ "2",
                            transition == "psiBR6.3" ~ "3",
                            transition == "psiBR6.4" ~ "4",
                            transition == "psiBR7.3" ~ "3",
                            transition == "psiBR7.4" ~ "4")) %>% 
  mutate(transition = str_replace(transition, "psiBR2.1", "BR-JB")) %>%
  mutate(transition = str_replace(transition, "psiBR3.1", "BR-MI")) %>%
  mutate(transition = str_replace(transition, "psiBR4.1", "BR-CC")) %>%
  mutate(transition = str_replace(transition, "psiBR4.2", "BR-CC")) %>%
  mutate(transition = str_replace(transition, "psiBR5.1", "BR-SE")) %>%
  mutate(transition = str_replace(transition, "psiBR5.2", "BR-SE")) %>%
  mutate(transition = str_replace(transition, "psiBR5.3", "BR-SE")) %>% 
  mutate(transition = str_replace(transition, "psiBR5.4", "BR-SE")) %>% 
  mutate(transition = str_replace(transition, "psiBR6.1", "BR-BR")) %>%
  mutate(transition = str_replace(transition, "psiBR6.2", "BR-BR")) %>%
  mutate(transition = str_replace(transition, "psiBR6.3", "BR-BR")) %>% 
  mutate(transition = str_replace(transition, "psiBR6.4", "BR-BR")) %>% 
  mutate(transition = str_replace(transition, "psiBR7.3", "BR-AR")) %>% 
  mutate(transition = str_replace(transition, "psiBR7.4", "BR-AR"))

psiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiAR1.4, psiAR2.1, psiAR3.1, psiAR4.1, psiAR5.1, psiAR5.4, psiAR6.1,
         psiAR6.4, psiAR7.4) %>% 
  pivot_longer(cols = 1:9, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "psiAR1.4" ~ "4",
                            transition == "psiAR2.1" ~ "1",
                            transition == "psiAR3.1" ~ "1",
                            transition == "psiAR4.1" ~ "1",
                            transition == "psiAR5.1" ~ "1",
                            transition == "psiAR5.4" ~ "4",
                            transition == "psiAR6.1" ~ "1",
                            transition == "psiAR6.4" ~ "4",
                            transition == "psiAR7.4" ~ "4")) %>% 
  mutate(transition = str_replace(transition, "psiAR1.4", "AR-DB")) %>% 
  mutate(transition = str_replace(transition, "psiAR2.1", "AR-JB")) %>%
  mutate(transition = str_replace(transition, "psiAR3.1", "AR-MI")) %>%
  mutate(transition = str_replace(transition, "psiAR4.1", "AR-CC")) %>%
  mutate(transition = str_replace(transition, "psiAR5.1", "AR-SE")) %>%
  mutate(transition = str_replace(transition, "psiAR5.4", "AR-SE")) %>% 
  mutate(transition = str_replace(transition, "psiAR6.1", "AR-BR")) %>%
  mutate(transition = str_replace(transition, "psiAR6.4", "AR-BR")) %>% 
  mutate(transition = str_replace(transition, "psiAR7.4", "AR-AR"))

# plot transitions
psiDB.plot <- ggplot() +
  geom_violin(psiDB.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiDB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiDB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiDB.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0, 1) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiJB.plot <- ggplot() +
  geom_violin(psiJB.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiJB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiJB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiJB.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0, 1) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiMI.plot <- ggplot() +
  geom_violin(psiMI.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiMI.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiMI.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiMI.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0, 1) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiCC.plot <- ggplot() +
  geom_violin(psiCC.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiCC.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiCC.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiCC.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0, 1) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiSE.plot <- ggplot() +
  geom_violin(psiSE.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiSE.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiSE.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiSE.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0, 1) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiBR.plot <- ggplot() +
  geom_violin(psiBR.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiBR.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiBR.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiBR.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0, 1) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiAR.plot <- ggplot() +
  geom_violin(psiAR.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiAR.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiAR.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiAR.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0, 1) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/debugging/debugging7/phiDBJBMICC.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiDBJBMICC.plot)

dev.off()

png(filename = "figures/debugging/debugging7/phiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiSE.plot)

dev.off()

png(filename = "figures/debugging/debugging7/phiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiBR.plot)

dev.off()

png(filename = "figures/debugging/debugging7/phiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiAR.plot)

dev.off()

# resighting
png(filename = "figures/debugging/debugging7/pDBJBMICC.png", width = 8, height = 8,
    units = "in", res = 600)

print(pDBJBMICC.plot)

dev.off()

png(filename = "figures/debugging/debugging7/pSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(pSE.plot)

dev.off()

png(filename = "figures/debugging/debugging7/pBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(pBR.plot)

dev.off()

png(filename = "figures/debugging/debugging7/pAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(pAR.plot)

dev.off()

# transition
png(filename = "figures/debugging/debugging7/psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/debugging/debugging7/psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/debugging/debugging7/psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

png(filename = "figures/debugging/debugging7/psiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiCC.plot)

dev.off()

png(filename = "figures/debugging/debugging7/psiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiSE.plot)

dev.off()

png(filename = "figures/debugging/debugging7/psiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiBR.plot)

dev.off()

png(filename = "figures/debugging/debugging7/psiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiAR.plot)

dev.off()

