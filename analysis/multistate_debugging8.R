################################################################################
# MULTISTATE SEASONAL SURVIVAL MODEL
# Gradual debugging
#
################################################################################

library(rjags)
library(jagsUI)
library(tidyverse)

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
ni <- 12000
nt <- 3
nb <- 6000
nc <- 3

# Call JAGS from R (BRT 2 days)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/debugging/multistate-debugging8-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/debugging/multistate-debugging8-simslist", Sys.Date(), ".rds"))

sim.summ <- sim.marr$summary

# parameter identifiability checks

sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/debugging8/",
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

png(filename = "figures/debugging/debugging8/ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()
