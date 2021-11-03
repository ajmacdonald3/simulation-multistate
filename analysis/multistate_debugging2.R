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
# Transitions: constant with all that never happen set to 0

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 40
n.years <- 10
n.states <- 8
n.obs <- 8

# DB
mean.phiDB1 <- 0.95
mean.phiDB2 <- 0.9
mean.phiDB3 <- 0.8
mean.phiDB4 <- 0.85

phiDB <- rep(c(mean.phiDB1, mean.phiDB2, mean.phiDB3, mean.phiDB4), length.out = n.occasions-1)

mean.pDB1 <- 0.65
mean.pDB2 <- 0.7
mean.pDB3 <- 0.5
mean.pDB4 <- 0.7

pDB <- rep(c(mean.pDB1, mean.pDB2, mean.pDB3, mean.pDB4), length.out = n.occasions-1)

# JB
mean.phiJB1 <- 0.7
mean.phiJB2 <- 0.9
mean.phiJB3 <- 0.8
mean.phiJB4 <- 0.7

phiJB <- rep(c(mean.phiJB1, mean.phiJB2, mean.phiJB3, mean.phiJB4), length.out = n.occasions-1)

mean.pJB1 <- 0.4
mean.pJB2 <- 0.6
mean.pJB3 <- 0.3
mean.pJB4 <- 0.7

pJB <- rep(c(mean.pJB1, mean.pJB2, mean.pJB3, mean.pJB4), length.out = n.occasions-1)

# MI
mean.phiMI1 <- 0.6
mean.phiMI2 <- 0.8
mean.phiMI3 <- 0.7
mean.phiMI4 <- 0.9

phiMI <- rep(c(mean.phiMI1, mean.phiMI2, mean.phiMI3, mean.phiMI4), length.out = n.occasions-1)

mean.pMI1 <- 0.5
mean.pMI2 <- 0.4
mean.pMI3 <- 0.7
mean.pMI4 <- 0.6

pMI <- rep(c(mean.pMI1, mean.pMI2, mean.pMI3, mean.pMI4), length.out = n.occasions-1)

# CC
mean.phiCC1 <- 0.75
mean.phiCC2 <- 0.85
mean.phiCC3 <- 0.85
mean.phiCC4 <- 0.8

phiCC <- rep(c(mean.phiCC1, mean.phiCC2, mean.phiCC3, mean.phiCC4), length.out = n.occasions-1)

mean.pCC1 <- 0.6
mean.pCC2 <- 0.6
mean.pCC3 <- 0.4
mean.pCC4 <- 0.5

pCC <- rep(c(mean.pCC1, mean.pCC2, mean.pCC3, mean.pCC4), length.out = n.occasions-1)

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
mean.phiAR1 <- 0.8
mean.phiAR2 <- 0.9
mean.phiAR3 <- 0.7
mean.phiAR4 <- 0.93

phiAR <- rep(c(mean.phiAR1, mean.phiAR2, mean.phiAR3, mean.phiAR4), length.out = n.occasions-1)

mean.pAR1 <- 0.5
mean.pAR2 <- 0.7
mean.pAR3 <- 0.6
mean.pAR4 <- 0.4

pAR <- rep(c(mean.pAR1, mean.pAR2, mean.pAR3, mean.pAR4), length.out = n.occasions-1)

# transition probabilities (sum to 1)
psiDB.DB <- 0
psiDB.JB <- 0.25
psiDB.MI <- 0.15
psiDB.CC <- 0.1
psiDB.SE <- 0.3
psiDB.BR <- 0.2
psiDB.AR <- 0

psiJB.DB <- 0
psiJB.JB <- 0
psiJB.MI <- 0
psiJB.CC <- 0.25
psiJB.SE <- 0.35
psiJB.BR <- 0.4
psiJB.AR <- 0

psiMI.DB <- 0
psiMI.JB <- 0
psiMI.MI <- 0
psiMI.CC <- 0.4
psiMI.SE <- 0.4
psiMI.BR <- 0.2
psiMI.AR <- 0

psiCC.DB <- 0
psiCC.JB <- 0
psiCC.MI <- 0
psiCC.CC <- 0.2
psiCC.SE <- 0.3
psiCC.BR <- 0.1
psiCC.AR <- 0.4

psiSE.DB <- 0.2
psiSE.JB <- 0.1
psiSE.MI <- 0.1
psiSE.CC <- 0.05
psiSE.SE <- 0.25
psiSE.BR <- 0.15
psiSE.AR <- 0.15

psiBR.DB <- 0.15
psiBR.JB <- 0.05
psiBR.MI <- 0.1
psiBR.CC <- 0.1
psiBR.SE <- 0.25
psiBR.BR <- 0.2
psiBR.AR <- 0.15

psiAR.DB <- 0.25
psiAR.JB <- 0.1
psiAR.MI <- 0.1
psiAR.CC <- 0.05
psiAR.SE <- 0.15
psiAR.BR <- 0.2
psiAR.AR <- 0.15

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(10, n.occasions) # DB
marked[,2] <- rep(10, n.occasions) # JB
marked[,3] <- rep(10, n.occasions) # MI
marked[,4] <- rep(10, n.occasions) # CC
marked[,5] <- rep(10, n.occasions) # SE
marked[,6] <- rep(10, n.occasions) # BR
marked[,7] <- rep(10, n.occasions) # AR
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
      phiDB[t]*psiDB.DB, phiDB[t]*psiDB.JB, phiDB[t]*psiDB.MI, phiDB[t]*psiDB.CC, phiDB[t]*psiDB.SE, phiDB[t]*psiDB.BR, phiDB[t]*psiDB.AR, 1-phiDB[t],
      phiJB[t]*psiJB.DB, phiJB[t]*psiJB.JB, phiJB[t]*psiJB.MI, phiJB[t]*psiJB.CC, phiJB[t]*psiJB.SE, phiJB[t]*psiJB.BR, phiJB[t]*psiJB.AR, 1-phiJB[t],
      phiMI[t]*psiMI.DB, phiMI[t]*psiMI.JB, phiMI[t]*psiMI.MI, phiMI[t]*psiMI.CC, phiMI[t]*psiMI.SE, phiMI[t]*psiMI.BR, phiMI[t]*psiMI.AR, 1-phiMI[t],
      phiCC[t]*psiCC.DB, phiCC[t]*psiCC.JB, phiCC[t]*psiCC.MI, phiCC[t]*psiCC.CC, phiCC[t]*psiCC.SE, phiCC[t]*psiCC.BR, phiCC[t]*psiCC.AR, 1-phiCC[t],
      phiSE[t]*psiSE.DB, phiSE[t]*psiSE.JB, phiSE[t]*psiSE.MI, phiSE[t]*psiSE.CC, phiSE[t]*psiSE.SE, phiSE[t]*psiSE.BR, phiSE[t]*psiSE.AR, 1-phiSE[t],
      phiBR[t]*psiBR.DB, phiBR[t]*psiBR.JB, phiBR[t]*psiBR.MI, phiBR[t]*psiBR.CC, phiBR[t]*psiBR.SE, phiBR[t]*psiBR.BR, phiBR[t]*psiBR.AR, 1-phiBR[t],
      phiAR[t]*psiAR.DB, phiAR[t]*psiAR.JB, phiAR[t]*psiAR.MI, phiAR[t]*psiAR.CC, phiAR[t]*psiAR.SE, phiAR[t]*psiAR.BR, phiAR[t]*psiAR.AR, 1-phiAR[t],
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
  
  db[1] <- 0
  psiDB[1] <- (db[1]/sum(db[]))
  
  for (i in 2:6){
  
    db[i] ~ dgamma(1, 1)
    psiDB[i] <- (db[i]/sum(db[]))
    
  }
  
  db[7] <- 0
  psiDB[7] <- (db[7]/sum(db[]))
  
  for (i in 1:3){
    jb[i] <- 0
    psiJB[i] <- (jb[i]/sum(jb[]))
  }
  
  for (i in 4:6){  
    jb[i] ~ dgamma(1, 1)
    psiJB[i] <- (jb[i]/sum(jb[]))
  }
  
    jb[7] <- 0
    psiJB[7] <- (jb[7]/sum(jb[]))
    
  for (i in 1:3){
    mi[i] <- 0
    psiMI[i] <- (mi[i]/sum(mi[]))
  }
  
  for (i in 4:6){  
    mi[i] ~ dgamma(1, 1)
    psiMI[i] <- (mi[i]/sum(mi[]))
  }
  
    mi[7] <- 0
    psiMI[7] <- (mi[7]/sum(mi[]))
  
  for (i in 1:3){
    cc[i] <- 0
    psiCC[i] <- (cc[i]/sum(cc[]))
  }
  
  for (i in 4:7){  
    cc[i] ~ dgamma(1, 1)
    psiCC[i] <- (cc[i]/sum(cc[]))
  }
  
  for (i in 1:7){
    se[i] ~ dgamma(1, 1)
    psiSE[i] <- (se[i]/sum(se[]))
  
    br[i] ~ dgamma(1, 1)
    psiBR[i] <- (br[i]/sum(br[]))
  
    ar[i] ~ dgamma(1, 1)
    psiAR[i] <- (ar[i]/sum(ar[]))
    
  } # i
  
  # Define state-transition and reencounter probabilities - note no i index - is no longer individual 
  for (t in 1:(n.occasions-1)){
    psi[1,t,1] <- phiDB[t] * psiDB[1]
    psi[1,t,2] <- phiDB[t] * psiDB[2]
    psi[1,t,3] <- phiDB[t] * psiDB[3]
    psi[1,t,4] <- phiDB[t] * psiDB[4]
    psi[1,t,5] <- phiDB[t] * psiDB[5]
    psi[1,t,6] <- phiDB[t] * psiDB[6]
    psi[1,t,7] <- phiDB[t] * psiDB[7]
    
    psi[2,t,1] <- phiJB[t] * psiJB[1]
    psi[2,t,2] <- phiJB[t] * psiJB[2]
    psi[2,t,3] <- phiJB[t] * psiJB[3]
    psi[2,t,4] <- phiJB[t] * psiJB[4]
    psi[2,t,5] <- phiJB[t] * psiJB[5]
    psi[2,t,6] <- phiJB[t] * psiJB[6]
    psi[2,t,7] <- phiJB[t] * psiJB[7]
    
    psi[3,t,1] <- phiMI[t] * psiMI[1]
    psi[3,t,2] <- phiMI[t] * psiMI[2]
    psi[3,t,3] <- phiMI[t] * psiMI[3]
    psi[3,t,4] <- phiMI[t] * psiMI[4]
    psi[3,t,5] <- phiMI[t] * psiMI[5]
    psi[3,t,6] <- phiMI[t] * psiMI[6]
    psi[3,t,7] <- phiMI[t] * psiMI[7]
    
    psi[4,t,1] <- phiCC[t] * psiCC[1]
    psi[4,t,2] <- phiCC[t] * psiCC[2]
    psi[4,t,3] <- phiCC[t] * psiCC[3]
    psi[4,t,4] <- phiCC[t] * psiCC[4]
    psi[4,t,5] <- phiCC[t] * psiCC[5]
    psi[4,t,6] <- phiCC[t] * psiCC[6]
    psi[4,t,7] <- phiCC[t] * psiCC[7]
    
    psi[5,t,1] <- phiSE[t] * psiSE[1]
    psi[5,t,2] <- phiSE[t] * psiSE[2]
    psi[5,t,3] <- phiSE[t] * psiSE[3]
    psi[5,t,4] <- phiSE[t] * psiSE[4]
    psi[5,t,5] <- phiSE[t] * psiSE[5]
    psi[5,t,6] <- phiSE[t] * psiSE[6]
    psi[5,t,7] <- phiSE[t] * psiSE[7]
    
    psi[6,t,1] <- phiBR[t] * psiBR[1]
    psi[6,t,2] <- phiBR[t] * psiBR[2]
    psi[6,t,3] <- phiBR[t] * psiBR[3]
    psi[6,t,4] <- phiBR[t] * psiBR[4]
    psi[6,t,5] <- phiBR[t] * psiBR[5]
    psi[6,t,6] <- phiBR[t] * psiBR[6]
    psi[6,t,7] <- phiBR[t] * psiBR[7]
    
    psi[7,t,1] <- phiAR[t] * psiAR[1]
    psi[7,t,2] <- phiAR[t] * psiAR[2]
    psi[7,t,3] <- phiAR[t] * psiAR[3]
    psi[7,t,4] <- phiAR[t] * psiAR[4]
    psi[7,t,5] <- phiAR[t] * psiAR[5]
    psi[7,t,6] <- phiAR[t] * psiAR[6]
    psi[7,t,7] <- phiAR[t] * psiAR[7]
    
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
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 2 days)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/debugging/multistate-debugging2-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/debugging/multistate-debugging2-simslist", Sys.Date(), ".rds"))

# parameter identifiability checks
summary <- readRDS("./analysis-output/debugging/multistate-debugging2-summary2021-11-01.rds")
sims.list <- readRDS("./analysis-output/debugging/multistate-debugging2-simslist2021-11-01.rds")

sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/debugging2/",
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

png(filename = "figures/debugging/debugging2/multistate-debugging2.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phiDB.sim <- tibble(site = c("DB-PreB", "DB-PostB1", "DB-PostB2", "DB-NonB"),
                    value = c(mean.phiDB1, mean.phiDB2, mean.phiDB3, mean.phiDB4))

phiJB.sim <- tibble(site = c("JB-PreB", "JB-PostB1", "JB-PostB2", "JB-NonB"),
                    value = c(mean.phiJB1, mean.phiJB2, mean.phiJB3, mean.phiJB4))

phiMI.sim <- tibble(site = c("MI-PreB", "MI-PostB1", "MI-PostB2", "MI-NonB"),
                    value = c(mean.phiMI1, mean.phiMI2, mean.phiMI3, mean.phiMI4))

phiCC.sim <- tibble(site = c("CC-PreB", "CC-PostB1", "CC-PostB2", "CC-NonB"),
                    value = c(mean.phiCC1, mean.phiCC2, mean.phiCC3, mean.phiCC4))

phiSE.sim <- tibble(site = c("SE-PreB", "SE-PostB1", "SE-PostB2", "SE-NonB"),
                    value = c(mean.phiSE1, mean.phiSE2, mean.phiSE3, mean.phiSE4))

phiBR.sim <- tibble(site = c("BR-PreB", "BR-PostB1", "BR-PostB2", "BR-NonB"),
                    value = c(mean.phiBR1, mean.phiBR2, mean.phiBR3, mean.phiBR4))

phiAR.sim <- tibble(site = c("AR-PreB", "AR-PostB1", "AR-PostB2", "AR-NonB"),
                    value = c(mean.phiAR1, mean.phiAR2, mean.phiAR3, mean.phiAR4))

# resighting
pDB.sim <- tibble(site = c("DB-PreB", "DB-PostB1", "DB-PostB2", "DB-NonB"),
                  value = c(mean.pDB1, mean.pDB2, mean.pDB3, mean.pDB4))

pJB.sim <- tibble(site = c("JB-PreB", "JB-PostB1", "JB-PostB2", "JB-NonB"),
                  value = c(mean.pJB1, mean.pJB2, mean.pJB3, mean.pJB4))

pMI.sim <- tibble(site = c("MI-PreB", "MI-PostB1", "MI-PostB2", "MI-NonB"),
                  value = c(mean.pMI1, mean.pMI2, mean.pMI3, mean.pMI4))

pCC.sim <- tibble(site = c("CC-PreB", "CC-PostB1", "CC-PostB2", "CC-NonB"),
                  value = c(mean.pCC1, mean.pCC2, mean.pCC3, mean.pCC4))

pSE.sim <- tibble(site = c("SE-PreB", "SE-PostB1", "SE-PostB2", "SE-NonB"),
                  value = c(mean.pSE1, mean.pSE2, mean.pSE3, mean.pSE4))

pBR.sim <- tibble(site = c("BR-PreB", "BR-PostB1", "BR-PostB2", "BR-NonB"),
                  value = c(mean.pBR1, mean.pBR2, mean.pBR3, mean.pBR4))

pAR.sim <- tibble(site = c("AR-PreB", "AR-PostB1", "AR-PostB2", "AR-NonB"),
                  value = c(mean.pAR1, mean.pAR2, mean.pAR3, mean.pAR4))

# transitions
psiDB.sim <- tibble(site = c("DB-DB", "DB-JB", "DB-MI", "DB-CC", "DB-SE", "DB-BR", "DB-AR"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI, psiDB.CC, psiDB.SE, psiDB.BR, psiDB.AR))

psiJB.sim <- tibble(site = c("JB-DB", "JB-JB", "JB-MI", "JB-CC", "JB-SE", "JB-BR", "JB-AR"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI, psiJB.CC, psiJB.SE, psiJB.BR, psiJB.AR))

psiMI.sim <- tibble(site = c("MI-DB", "MI-JB", "MI-MI", "MI-CC", "MI-SE", "MI-BR", "MI-AR"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI, psiMI.CC, psiMI.SE, psiMI.BR, psiMI.AR))

psiCC.sim <- tibble(site = c("CC-DB", "CC-JB", "CC-MI", "CC-CC", "CC-SE", "CC-BR", "CC-AR"),
                    value = c(psiCC.DB, psiCC.JB, psiCC.MI, psiCC.CC, psiCC.SE, psiCC.BR, psiCC.AR))

psiSE.sim <- tibble(site = c("SE-DB", "SE-JB", "SE-MI", "SE-CC", "SE-SE", "SE-BR", "SE-AR"),
                    value = c(psiSE.DB, psiSE.JB, psiSE.MI, psiSE.CC, psiSE.SE, psiSE.BR, psiSE.AR))

psiBR.sim <- tibble(site = c("BR-DB", "BR-JB", "BR-MI", "BR-CC", "BR-SE", "BR-BR", "BR-AR"),
                    value = c(psiBR.DB, psiBR.JB, psiBR.MI, psiBR.CC, psiBR.SE, psiBR.BR, psiBR.AR))

psiAR.sim <- tibble(site = c("AR-AR", "AR-JB", "AR-MI", "AR-CC", "AR-SE", "AR-BR", "AR-AR"),
                    value = c(psiAR.AR, psiAR.JB, psiAR.MI, psiAR.CC, psiAR.SE, psiAR.BR, psiAR.AR))

# format survival
phiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiDB.1, mean.phiDB.2, mean.phiDB.3, mean.phiDB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiDB.1", "DB-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiDB.2", "DB-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiDB.3", "DB-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.phiDB.4", "DB-NonB"))

phiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiJB.1, mean.phiJB.2, mean.phiJB.3, mean.phiJB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiJB.1", "JB-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiJB.2", "JB-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiJB.3", "JB-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.phiJB.4", "JB-NonB"))

phiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiMI.1, mean.phiMI.2, mean.phiMI.3, mean.phiMI.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiMI.1", "MI-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiMI.2", "MI-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiMI.3", "MI-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.phiMI.4", "MI-NonB"))

phiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiCC.1, mean.phiCC.2, mean.phiCC.3, mean.phiCC.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiCC.1", "CC-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiCC.2", "CC-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiCC.3", "CC-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.phiCC.4", "CC-NonB"))

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
  select(mean.phiAR.1, mean.phiAR.2, mean.phiAR.3, mean.phiAR.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiAR.1", "AR-PreB")) %>%
  mutate(site = str_replace(site, "mean.phiAR.2", "AR-PostB1")) %>%
  mutate(site = str_replace(site, "mean.phiAR.3", "AR-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.phiAR.4", "AR-NonB"))

# plot survival
phiDB.plot <- ggplot() +
  geom_violin(phiDB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiDB.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiDB.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiDB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

phiJB.plot <- ggplot() +
  geom_violin(phiJB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiJB.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiJB.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiJB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

phiMI.plot <- ggplot() +
  geom_violin(phiMI.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiMI.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiMI.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiMI.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

phiCC.plot <- ggplot() +
  geom_violin(phiCC.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phiCC.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phiCC.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phiCC.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
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
pDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pDB.1, mean.pDB.2, mean.pDB.3, mean.pDB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pDB.1", "DB-PreB")) %>%
  mutate(site = str_replace(site, "mean.pDB.2", "DB-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pDB.3", "DB-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.pDB.4", "DB-NonB"))

pJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pJB.1, mean.pJB.2, mean.pJB.3, mean.pJB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pJB.1", "JB-PreB")) %>%
  mutate(site = str_replace(site, "mean.pJB.2", "JB-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pJB.3", "JB-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.pJB.4", "JB-NonB"))

pMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pMI.1, mean.pMI.2, mean.pMI.3, mean.pMI.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pMI.1", "MI-PreB")) %>%
  mutate(site = str_replace(site, "mean.pMI.2", "MI-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pMI.3", "MI-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.pMI.4", "MI-NonB"))

pCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pCC.1, mean.pCC.2, mean.pCC.3, mean.pCC.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pCC.1", "CC-PreB")) %>%
  mutate(site = str_replace(site, "mean.pCC.2", "CC-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pCC.3", "CC-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.pCC.4", "CC-NonB"))

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
  select(mean.pAR.1, mean.pAR.2, mean.pAR.3, mean.pAR.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pAR.1", "AR-PreB")) %>%
  mutate(site = str_replace(site, "mean.pAR.2", "AR-PostB1")) %>%
  mutate(site = str_replace(site, "mean.pAR.3", "AR-PostB2")) %>% 
  mutate(site = str_replace(site, "mean.pAR.4", "AR-NonB"))

# plot resighting
pDB.plot <- ggplot() +
  geom_violin(pDB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pDB.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pDB.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pDB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

pJB.plot <- ggplot() +
  geom_violin(pJB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pJB.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pJB.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pJB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

pMI.plot <- ggplot() +
  geom_violin(pMI.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pMI.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pMI.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pMI.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

pCC.plot <- ggplot() +
  geom_violin(pCC.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(pCC.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(pCC.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(pCC.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
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
  ylab("Survival probability") +
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
  ylab("Survival probability") +
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
  ylab("Survival probability") +
  theme(legend.position = "none")

# format transitions
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiDB.5, psiDB.6, psiDB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiDB.1", "DB-DB")) %>%
  mutate(site = str_replace(site, "psiDB.2", "DB-JB")) %>%
  mutate(site = str_replace(site, "psiDB.3", "DB-MI")) %>% 
  mutate(site = str_replace(site, "psiDB.4", "DB-CC")) %>%
  mutate(site = str_replace(site, "psiDB.5", "DB-SE")) %>%
  mutate(site = str_replace(site, "psiDB.6", "DB-BR")) %>%
  mutate(site = str_replace(site, "psiDB.7", "DB-AR"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiJB.1, psiJB.2, psiJB.3, psiJB.4, psiJB.5, psiJB.6, psiJB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiJB.1", "JB-DB")) %>%
  mutate(site = str_replace(site, "psiJB.2", "JB-JB")) %>%
  mutate(site = str_replace(site, "psiJB.3", "JB-MI")) %>% 
  mutate(site = str_replace(site, "psiJB.4", "JB-CC")) %>%
  mutate(site = str_replace(site, "psiJB.5", "JB-SE")) %>%
  mutate(site = str_replace(site, "psiJB.6", "JB-BR")) %>%
  mutate(site = str_replace(site, "psiJB.7", "JB-AR"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiMI.5, psiMI.6, psiMI.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiMI.1", "MI-DB")) %>%
  mutate(site = str_replace(site, "psiMI.2", "MI-JB")) %>%
  mutate(site = str_replace(site, "psiMI.3", "MI-MI")) %>% 
  mutate(site = str_replace(site, "psiMI.4", "MI-CC")) %>%
  mutate(site = str_replace(site, "psiMI.5", "MI-SE")) %>%
  mutate(site = str_replace(site, "psiMI.6", "MI-BR")) %>%
  mutate(site = str_replace(site, "psiMI.7", "MI-AR"))

psiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiCC.1, psiCC.2, psiCC.3, psiCC.4, psiCC.5, psiCC.6, psiCC.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiCC.1", "CC-DB")) %>%
  mutate(site = str_replace(site, "psiCC.2", "CC-JB")) %>%
  mutate(site = str_replace(site, "psiCC.3", "CC-MI")) %>% 
  mutate(site = str_replace(site, "psiCC.4", "CC-CC")) %>%
  mutate(site = str_replace(site, "psiCC.5", "CC-SE")) %>%
  mutate(site = str_replace(site, "psiCC.6", "CC-BR")) %>%
  mutate(site = str_replace(site, "psiCC.7", "CC-AR"))

psiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiSE.1, psiSE.2, psiSE.3, psiSE.4, psiSE.5, psiSE.6, psiSE.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiSE.1", "SE-DB")) %>%
  mutate(site = str_replace(site, "psiSE.2", "SE-JB")) %>%
  mutate(site = str_replace(site, "psiSE.3", "SE-MI")) %>% 
  mutate(site = str_replace(site, "psiSE.4", "SE-CC")) %>%
  mutate(site = str_replace(site, "psiSE.5", "SE-SE")) %>%
  mutate(site = str_replace(site, "psiSE.6", "SE-BR")) %>%
  mutate(site = str_replace(site, "psiSE.7", "SE-AR"))

psiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiBR.1, psiBR.2, psiBR.3, psiBR.4, psiBR.5, psiBR.6, psiBR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiBR.1", "BR-DB")) %>%
  mutate(site = str_replace(site, "psiBR.2", "BR-JB")) %>%
  mutate(site = str_replace(site, "psiBR.3", "BR-MI")) %>% 
  mutate(site = str_replace(site, "psiBR.4", "BR-CC")) %>%
  mutate(site = str_replace(site, "psiBR.5", "BR-SE")) %>%
  mutate(site = str_replace(site, "psiBR.6", "BR-BR")) %>%
  mutate(site = str_replace(site, "psiBR.7", "BR-AR"))

psiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiAR.1, psiAR.2, psiAR.3, psiAR.4, psiAR.5, psiAR.6, psiAR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiAR.1", "AR-AR")) %>%
  mutate(site = str_replace(site, "psiAR.2", "AR-JB")) %>%
  mutate(site = str_replace(site, "psiAR.3", "AR-MI")) %>% 
  mutate(site = str_replace(site, "psiAR.4", "AR-CC")) %>%
  mutate(site = str_replace(site, "psiAR.5", "AR-SE")) %>%
  mutate(site = str_replace(site, "psiAR.6", "AR-BR")) %>%
  mutate(site = str_replace(site, "psiAR.7", "AR-AR"))

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
png(filename = "figures/debugging/debugging2/phiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiDB.plot)

dev.off()

png(filename = "figures/debugging/debugging2/phiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiJB.plot)

dev.off()

png(filename = "figures/debugging/debugging2/phiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiMI.plot)

dev.off()

png(filename = "figures/debugging/debugging2/phiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiCC.plot)

dev.off()

png(filename = "figures/debugging/debugging2/phiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiSE.plot)

dev.off()

png(filename = "figures/debugging/debugging2/phiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiBR.plot)

dev.off()

png(filename = "figures/debugging/debugging2/phiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiAR.plot)

dev.off()

# resighting
png(filename = "figures/debugging/debugging2/pDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pDB.plot)

dev.off()

png(filename = "figures/debugging/debugging2/pJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pJB.plot)

dev.off()

png(filename = "figures/debugging/debugging2/pMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(pMI.plot)

dev.off()

png(filename = "figures/debugging/debugging2/pCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(pCC.plot)

dev.off()

png(filename = "figures/debugging/debugging2/pSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(pSE.plot)

dev.off()

png(filename = "figures/debugging/debugging2/pBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(pBR.plot)

dev.off()

png(filename = "figures/debugging/debugging2/pAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(pAR.plot)

dev.off()

# transition
png(filename = "figures/debugging/debugging2/psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/debugging/debugging2/psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/debugging/debugging2/psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

png(filename = "figures/debugging/debugging2/psiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiCC.plot)

dev.off()

png(filename = "figures/debugging/debugging2/psiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiSE.plot)

dev.off()

png(filename = "figures/debugging/debugging2/psiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiBR.plot)

dev.off()

png(filename = "figures/debugging/debugging2/psiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiAR.plot)

dev.off()
