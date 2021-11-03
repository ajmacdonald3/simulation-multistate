################################################################################
# MULTISTATE SEASONAL SURVIVAL MODEL
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

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 30
n.years <- 10
n.states <- 6
n.obs <- 6

# DB
mean.phiDB1 <- 0.85
mean.phiDB2 <- 0
mean.phiDB3 <- 0
var.phiDB <- 0.3                       # Temporal variance of survival

mean.pDB1 <- 0
mean.pDB2 <- 0
mean.pDB3 <- 0.7
var.pDB <- 0.5

# Determine annual survival probabilities
logit.phiDB1 <- rnorm(n.years, qlogis(mean.phiDB1), var.phiDB^0.5)
phiDB1 <- plogis(logit.phiDB1)

logit.phiDB2 <- rnorm(n.years, qlogis(mean.phiDB2), var.phiDB^0.5)
phiDB2 <- plogis(logit.phiDB2)

logit.phiDB3 <- rnorm(n.years, qlogis(mean.phiDB3), var.phiDB^0.5)
phiDB3 <- plogis(logit.phiDB3)

phiDB <- c(rbind(phiDB1, phiDB2, phiDB3))
phiDB <- phiDB[1:29]

# Determine annual resighting probabilities
logit.pDB1 <- rnorm(n.years, qlogis(mean.pDB1), var.pDB^0.5)
pDB1 <- plogis(logit.pDB1)

logit.pDB2 <- rnorm(n.years, qlogis(mean.pDB2), var.pDB^0.5)
pDB2 <- plogis(logit.pDB2)

logit.pDB3 <- rnorm(n.years, qlogis(mean.pDB3), var.pDB^0.5)
pDB3 <- plogis(logit.pDB3)

pDB <- c(rbind(pDB1, pDB2, pDB3))
pDB <- pDB[1:29]

# JB
mean.phiJB1 <- 0
mean.phiJB2 <- 0.95
mean.phiJB3 <- 0
var.phiJB <- 0.3                       # Temporal variance of survival

mean.pJB1 <- 0.4
mean.pJB2 <- 0
mean.pJB3 <- 0
var.pJB <- 0.5

# Determine annual survival probabilities
logit.phiJB1 <- rnorm(n.years, qlogis(mean.phiJB1), var.phiJB^0.5)
phiJB1 <- plogis(logit.phiJB1)

logit.phiJB2 <- rnorm(n.years, qlogis(mean.phiJB2), var.phiJB^0.5)
phiJB2 <- plogis(logit.phiJB2)

logit.phiJB3 <- rnorm(n.years, qlogis(mean.phiJB3), var.phiJB^0.5)
phiJB3 <- plogis(logit.phiJB3)

phiJB <- c(rbind(phiJB1, phiJB2, phiJB3))
phiJB <- phiJB[1:29]

# Determine annual resighting probabilities
logit.pJB1 <- rnorm(n.years, qlogis(mean.pJB1), var.pJB^0.5)
pJB1 <- plogis(logit.pJB1)

logit.pJB2 <- rnorm(n.years, qlogis(mean.pJB2), var.pJB^0.5)
pJB2 <- plogis(logit.pJB2)

logit.pJB3 <- rnorm(n.years, qlogis(mean.pJB3), var.pJB^0.5)
pJB3 <- plogis(logit.pJB3)

pJB <- c(rbind(pJB1, pJB2, pJB3))
pJB <- pJB[1:29]

# MI
mean.phiMI1 <- 0
mean.phiMI2 <- 0.9
mean.phiMI3 <- 0
var.phiMI <- 0.3                       # Temporal variance of survival

mean.pMI1 <- 0.5
mean.pMI2 <- 0
mean.pMI3 <- 0
var.pMI <- 0.5

# Determine annual survival probabilities
logit.phiMI1 <- rnorm(n.years, qlogis(mean.phiMI1), var.phiMI^0.5)
phiMI1 <- plogis(logit.phiMI1)

logit.phiMI2 <- rnorm(n.years, qlogis(mean.phiMI2), var.phiMI^0.5)
phiMI2 <- plogis(logit.phiMI2)

logit.phiMI3 <- rnorm(n.years, qlogis(mean.phiMI3), var.phiMI^0.5)
phiMI3 <- plogis(logit.phiMI3)

phiMI <- c(rbind(phiMI1, phiMI2, phiMI3))
phiMI <- phiMI[1:29]

# Determine annual resighting probabilities
logit.pMI1 <- rnorm(n.years, qlogis(mean.pMI1), var.pMI^0.5)
pMI1 <- plogis(logit.pMI1)

logit.pMI2 <- rnorm(n.years, qlogis(mean.pMI2), var.pMI^0.5)
pMI2 <- plogis(logit.pMI2)

logit.pMI3 <- rnorm(n.years, qlogis(mean.pMI3), var.pMI^0.5)
pMI3 <- plogis(logit.pMI3)

pMI <- c(rbind(pMI1, pMI2, pMI3))
pMI <- pMI[1:29]

# SE
mean.phiSE1 <- 0.8
mean.phiSE2 <- 0.93
mean.phiSE3 <- 0.9
var.phiSE <- 0.3                       # Temporal variance of survival

mean.pSE1 <- 0.6
mean.pSE2 <- 0.6
mean.pSE3 <- 0.4
var.pSE <- 0.5

# Determine annual survival probabilities
logit.phiSE1 <- rnorm(n.years, qlogis(mean.phiSE1), var.phiSE^0.5)
phiSE1 <- plogis(logit.phiSE1)

logit.phiSE2 <- rnorm(n.years, qlogis(mean.phiSE2), var.phiSE^0.5)
phiSE2 <- plogis(logit.phiSE2)

logit.phiSE3 <- rnorm(n.years, qlogis(mean.phiSE3), var.phiSE^0.5)
phiSE3 <- plogis(logit.phiSE3)

phiSE <- c(rbind(phiSE1, phiSE2, phiSE3))
phiSE <- phiSE[1:29]

# Determine annual resighting probabilities
logit.pSE1 <- rnorm(n.years, qlogis(mean.pSE1), var.pSE^0.5)
pSE1 <- plogis(logit.pSE1)

logit.pSE2 <- rnorm(n.years, qlogis(mean.pSE2), var.pSE^0.5)
pSE2 <- plogis(logit.pSE2)

logit.pSE3 <- rnorm(n.years, qlogis(mean.pSE3), var.pSE^0.5)
pSE3 <- plogis(logit.pSE3)

pSE <- c(rbind(pSE1, pSE2, pSE3))
pSE <- pSE[1:29]

# BR
mean.phiBR1 <- 0
mean.phiBR2 <- 0
mean.phiBR3 <- 0.92
var.phiBR <- 0.3                       # Temporal variance of survival

mean.pBR1 <- 0
mean.pBR2 <- 0.3
mean.pBR3 <- 0
var.pBR <- 0.5

# Determine annual survival probabilities
logit.phiBR1 <- rnorm(n.years, qlogis(mean.phiBR1), var.phiBR^0.5)
phiBR1 <- plogis(logit.phiBR1)

logit.phiBR2 <- rnorm(n.years, qlogis(mean.phiBR2), var.phiBR^0.5)
phiBR2 <- plogis(logit.phiBR2)

logit.phiBR3 <- rnorm(n.years, qlogis(mean.phiBR3), var.phiBR^0.5)
phiBR3 <- plogis(logit.phiBR3)

phiBR <- c(rbind(phiBR1, phiBR2, phiBR3))
phiBR <- phiBR[1:29]

# Determine annual resighting probabilities
logit.pBR1 <- rnorm(n.years, qlogis(mean.pBR1), var.pBR^0.5)
pBR1 <- plogis(logit.pBR1)

logit.pBR2 <- rnorm(n.years, qlogis(mean.pBR2), var.pBR^0.5)
pBR2 <- plogis(logit.pBR2)

logit.pBR3 <- rnorm(n.years, qlogis(mean.pBR3), var.pBR^0.5)
pBR3 <- plogis(logit.pBR3)

pBR <- c(rbind(pBR1, pBR2, pBR3))
pBR <- pBR[1:29]

# transition probabilities (sum to 1)
psiDB.DB <- rep(0, n.occasions)
psiDB.JB <- rep(c(0.5, 0, 0), length.out = n.occasions)
psiDB.MI <- rep(c(0.4, 0, 0), length.out = n.occasions)
psiDB.SE <- rep(c(0.1, 0, 0), length.out = n.occasions)
psiDB.BR <- rep(0, n.occasions)

psiJB.DB <- rep(0, n.occasions)
psiJB.JB <- rep(0, n.occasions)
psiJB.MI <- rep(0, n.occasions)
psiJB.SE <- rep(c(0, 0.4, 0), length.out = n.occasions)
psiJB.BR <- rep(c(0, 0.6, 0), length.out = n.occasions)

psiMI.DB <- rep(0, n.occasions)
psiMI.JB <- rep(0, n.occasions)
psiMI.MI <- rep(0, n.occasions)
psiMI.SE <- rep(c(0, 0.6, 0), length.out = n.occasions)
psiMI.BR <- rep(c(0, 0.4, 0), length.out = n.occasions)

psiSE.DB <- rep(c(0, 0, 0.8), length.out = n.occasions)
psiSE.JB <- rep(c(0.4, 0, 0), length.out = n.occasions)
psiSE.MI <- rep(c(0.5, 0, 0), length.out = n.occasions)
psiSE.SE <- rep(c(0.1, 0.7, 0.2), length.out = n.occasions)
psiSE.BR <- rep(c(0, 0.3, 0), length.out = n.occasions)

psiBR.DB <- rep(c(0, 0, 0.7), length.out = n.occasions)
psiBR.JB <- rep(0, n.occasions)
psiBR.MI <- rep(0, n.occasions)
psiBR.SE <- rep(c(0, 0, 0.3), length.out = n.occasions)
psiBR.BR <- rep(0, n.occasions)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(100, 0, 0), length.out = n.occasions) # DB
marked[,2] <- rep(c(0, 50, 0), length.out = n.occasions) # JB
marked[,3] <- rep(c(0, 50, 0), length.out = n.occasions) # MI
marked[,4] <- rep(20, n.occasions) # SE
marked[,5] <- rep(c(0, 0, 10), length.out = n.occasions) # BR
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
      phiDB[t]*psiDB.DB[t], phiDB[t]*psiDB.JB[t], phiDB[t]*psiDB.MI[t], phiDB[t]*psiDB.SE[t], phiDB[t]*psiDB.BR[t], 1-phiDB[t],
      phiJB[t]*psiJB.DB[t], phiJB[t]*psiJB.JB[t], phiJB[t]*psiJB.MI[t], phiJB[t]*psiJB.SE[t], phiJB[t]*psiJB.BR[t], 1-phiJB[t],
      phiMI[t]*psiMI.DB[t], phiMI[t]*psiMI.JB[t], phiMI[t]*psiMI.MI[t], phiMI[t]*psiMI.SE[t], phiMI[t]*psiMI.BR[t], 1-phiMI[t],
      phiSE[t]*psiSE.DB[t], phiSE[t]*psiSE.JB[t], phiSE[t]*psiSE.MI[t], phiSE[t]*psiSE.SE[t], phiSE[t]*psiSE.BR[t], 1-phiSE[t],
      phiBR[t]*psiBR.DB[t], phiBR[t]*psiBR.JB[t], phiBR[t]*psiBR.MI[t], phiBR[t]*psiBR.SE[t], phiBR[t]*psiBR.BR[t], 1-phiBR[t],
      0,                    0,                    0,                    0,                    0,                    1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pDB[t], 0,      0,      0,      0,      1-pDB[t],
      0,      pJB[t], 0,      0,      0,      1-pJB[t],
      0,      0,      pMI[t], 0,      0,      1-pMI[t],
      0,      0,      0,      pSE[t], 0,      1-pSE[t],
      0,      0,      0,      0,      pBR[t], 1-pBR[t],
      0,      0,      0,      0,      0,      1),
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
  for (t in 1:10){
    logit(phiDB[t]) <- mu.phiDB + eps.phiDB[t]
    eps.phiDB[t] ~ dnorm(0, tau.phiDB)T(-10,10)
  }
  
  for (t in 1:9){
    logit(pDB[t]) <- mu.pDB + eps.pDB[t]
    eps.pDB[t] ~ dnorm(0, tau.pDB)T(-10,10)
  }
  
    mean.phiDB ~ dunif(0, 1)
    mu.phiDB <- logit(mean.phiDB)
    sigma.phiDB ~ dunif(0, 10)               # Prior for standard deviation
    tau.phiDB <- pow(sigma.phiDB, -2)
    sigma2.phiDB <- pow(sigma.phiDB, 2)
    
    mean.pDB ~ dunif(0, 1)
    mu.pDB <- logit(mean.pDB)
    sigma.pDB ~ dunif(0, 10)               # Prior for standard deviation
    tau.pDB <- pow(sigma.pDB, -2)
    sigma2.pDB <- pow(sigma.pDB, 2)
  
  phitDB[1] <- phiDB[1]
  phitDB[2] <- 0
  phitDB[3] <- 0
  phitDB[4] <- phiDB[2]
  phitDB[5] <- 0
  phitDB[6] <- 0
  phitDB[7] <- phiDB[3]
  phitDB[8] <- 0
  phitDB[9] <- 0
  phitDB[10] <- phiDB[4]
  phitDB[11] <- 0
  phitDB[12] <- 0
  phitDB[13] <- phiDB[5]
  phitDB[14] <- 0
  phitDB[15] <- 0
  phitDB[16] <- phiDB[6]
  phitDB[17] <- 0
  phitDB[18] <- 0
  phitDB[19] <- phiDB[7]
  phitDB[20] <- 0
  phitDB[21] <- 0
  phitDB[22] <- phiDB[8]
  phitDB[23] <- 0
  phitDB[24] <- 0
  phitDB[25] <- phiDB[9]
  phitDB[26] <- 0
  phitDB[27] <- 0
  phitDB[28] <- phiDB[10]
  phitDB[29] <- 0
  
  psightDB[1] <- 0
  psightDB[2] <- 0
  psightDB[3] <- pDB[1]
  psightDB[4] <- 0
  psightDB[5] <- 0
  psightDB[6] <- pDB[2]
  psightDB[7] <- 0
  psightDB[8] <- 0
  psightDB[9] <- pDB[3]
  psightDB[10] <- 0
  psightDB[11] <- 0
  psightDB[12] <- pDB[4]
  psightDB[13] <- 0
  psightDB[14] <- 0
  psightDB[15] <- pDB[5]
  psightDB[16] <- 0
  psightDB[17] <- 0
  psightDB[18] <- pDB[6]
  psightDB[19] <- 0
  psightDB[20] <- 0
  psightDB[21] <- pDB[7]
  psightDB[22] <- 0
  psightDB[23] <- 0
  psightDB[24] <- pDB[8]
  psightDB[25] <- 0
  psightDB[26] <- 0
  psightDB[27] <- pDB[9]
  psightDB[28] <- 0
  psightDB[29] <- 0
  
  # JB
  for (t in 1:10){
    logit(phiJB[t]) <- mu.phiJB + eps.phiJB[t]
    eps.phiJB[t] ~ dnorm(0, tau.phiJB)T(-10,10)
    
    logit(pJB[t]) <- mu.pJB + eps.pJB[t]
    eps.pJB[t] ~ dnorm(0, tau.pJB)T(-10,10)
  }
  
    mean.phiJB ~ dunif(0, 1)
    mu.phiJB <- logit(mean.phiJB)
    sigma.phiJB ~ dunif(0, 10)               # Prior for standard deviation
    tau.phiJB <- pow(sigma.phiJB, -2)
    sigma2.phiJB <- pow(sigma.phiJB, 2)
    
    mean.pJB ~ dunif(0, 1)
    mu.pJB <- logit(mean.pJB)
    sigma.pJB ~ dunif(0, 10)               # Prior for standard deviation
    tau.pJB <- pow(sigma.pJB, -2)
    sigma2.pJB <- pow(sigma.pJB, 2)

  phitJB[1] <- 0
  phitJB[2] <- phiJB[1]
  phitJB[3] <- 0
  phitJB[4] <- 0
  phitJB[5] <- phiJB[2]
  phitJB[6] <- 0
  phitJB[7] <- 0
  phitJB[8] <- phiJB[3]
  phitJB[9] <- 0
  phitJB[10] <- 0
  phitJB[11] <- phiJB[4]
  phitJB[12] <- 0
  phitJB[13] <- 0
  phitJB[14] <- phiJB[5]
  phitJB[15] <- 0
  phitJB[16] <- 0
  phitJB[17] <- phiJB[6]
  phitJB[18] <- 0
  phitJB[19] <- 0
  phitJB[20] <- phiJB[7]
  phitJB[21] <- 0
  phitJB[22] <- 0
  phitJB[23] <- phiJB[8]
  phitJB[24] <- 0
  phitJB[25] <- 0
  phitJB[26] <- phiJB[9]
  phitJB[27] <- 0
  phitJB[28] <- 0
  phitJB[29] <- phiJB[10]
  
  psightJB[1] <- pJB[1]
  psightJB[2] <- 0
  psightJB[3] <- 0
  psightJB[4] <- pJB[2]
  psightJB[5] <- 0
  psightJB[6] <- 0
  psightJB[7] <- pJB[3]
  psightJB[8] <- 0
  psightJB[9] <- 0
  psightJB[10] <- pJB[4]
  psightJB[11] <- 0
  psightJB[12] <- 0
  psightJB[13] <- pJB[5]
  psightJB[14] <- 0
  psightJB[15] <- 0
  psightJB[16] <- pJB[6]
  psightJB[17] <- 0
  psightJB[18] <- 0
  psightJB[19] <- pJB[7]
  psightJB[20] <- 0
  psightJB[21] <- 0
  psightJB[22] <- pJB[8]
  psightJB[23] <- 0
  psightJB[24] <- 0
  psightJB[25] <- pJB[9]
  psightJB[26] <- 0
  psightJB[27] <- 0
  psightJB[28] <- pJB[10]
  psightJB[29] <- 0

  # MI
  for (t in 1:10){
    logit(phiMI[t]) <- mu.phiMI + eps.phiMI[t]
    eps.phiMI[t] ~ dnorm(0, tau.phiMI)T(-10,10)
    
    logit(pMI[t]) <- mu.pMI + eps.pMI[t]
    eps.pMI[t] ~ dnorm(0, tau.pMI)T(-10,10)
  }
  
    mean.phiMI ~ dunif(0, 1)
    mu.phiMI <- logit(mean.phiMI)
    sigma.phiMI ~ dunif(0, 10)               # Prior for standard deviation
    tau.phiMI <- pow(sigma.phiMI, -2)
    sigma2.phiMI <- pow(sigma.phiMI, 2)
    
    mean.pMI ~ dunif(0, 1)
    mu.pMI <- logit(mean.pMI)
    sigma.pMI ~ dunif(0, 10)               # Prior for standard deviation
    tau.pMI <- pow(sigma.pMI, -2)
    sigma2.pMI <- pow(sigma.pMI, 2)

  phitMI[1] <- 0
  phitMI[2] <- phiMI[1]
  phitMI[3] <- 0
  phitMI[4] <- 0
  phitMI[5] <- phiMI[2]
  phitMI[6] <- 0
  phitMI[7] <- 0
  phitMI[8] <- phiMI[3]
  phitMI[9] <- 0
  phitMI[10] <- 0
  phitMI[11] <- phiMI[4]
  phitMI[12] <- 0
  phitMI[13] <- 0
  phitMI[14] <- phiMI[5]
  phitMI[15] <- 0
  phitMI[16] <- 0
  phitMI[17] <- phiMI[6]
  phitMI[18] <- 0
  phitMI[19] <- 0
  phitMI[20] <- phiMI[7]
  phitMI[21] <- 0
  phitMI[22] <- 0
  phitMI[23] <- phiMI[8]
  phitMI[24] <- 0
  phitMI[25] <- 0
  phitMI[26] <- phiMI[9]
  phitMI[27] <- 0
  phitMI[28] <- 0
  phitMI[29] <- phiMI[10]
  
  psightMI[1] <- pMI[1]
  psightMI[2] <- 0
  psightMI[3] <- 0
  psightMI[4] <- pMI[2]
  psightMI[5] <- 0
  psightMI[6] <- 0
  psightMI[7] <- pMI[3]
  psightMI[8] <- 0
  psightMI[9] <- 0
  psightMI[10] <- pMI[4]
  psightMI[11] <- 0
  psightMI[12] <- 0
  psightMI[13] <- pMI[5]
  psightMI[14] <- 0
  psightMI[15] <- 0
  psightMI[16] <- pMI[6]
  psightMI[17] <- 0
  psightMI[18] <- 0
  psightMI[19] <- pMI[7]
  psightMI[20] <- 0
  psightMI[21] <- 0
  psightMI[22] <- pMI[8]
  psightMI[23] <- 0
  psightMI[24] <- 0
  psightMI[25] <- pMI[9]
  psightMI[26] <- 0
  psightMI[27] <- 0
  psightMI[28] <- pMI[10]
  psightMI[29] <- 0

  # SE
  for (t in 1:(n.occasions-1)){
    logit(phiSE[t]) <- mu.phiSE[season[t]] + eps.phiSE[t]
    eps.phiSE[t] ~ dnorm(0, tau.phiSE[season[t]])T(-10,10)
    
    logit(pSE[t]) <- mu.pSE[season[t]] + eps.pSE[t]
    eps.pSE[t] ~ dnorm(0, tau.pSE[season[t]])T(-10,10)
  }
  
  for (s in 1:3){
    mean.phiSE[s] ~ dunif(0, 1)
    mu.phiSE[s] <- logit(mean.phiSE[s])
    sigma.phiSE[s] ~ dunif(0, 10)               # Prior for standard deviation
    tau.phiSE[s] <- pow(sigma.phiSE[s], -2)
    sigma2.phiSE[s] <- pow(sigma.phiSE[s], 2)
    
    mean.pSE[s] ~ dunif(0, 1)
    mu.pSE[s] <- logit(mean.pSE[s])
    sigma.pSE[s] ~ dunif(0, 10)               # Prior for standard deviation
    tau.pSE[s] <- pow(sigma.pSE[s], -2)
    sigma2.pSE[s] <- pow(sigma.pSE[s], 2)
  }
  
  # BR
  for (t in 1:10){
    logit(phiBR[t]) <- mu.phiBR + eps.phiBR[t]
    eps.phiBR[t] ~ dnorm(0, tau.phiBR)T(-15,15)
    
    logit(pBR[t]) <- mu.pBR + eps.pBR[t]
    eps.pBR[t] ~ dnorm(0, tau.pBR)T(-15,15)
  }
  
    mean.phiBR ~ dunif(0, 1)
    mu.phiBR <- logit(mean.phiBR)
    sigma.phiBR ~ dunif(0, 10)               # Prior for standard deviation
    tau.phiBR <- pow(sigma.phiBR, -2)
    sigma2.phiBR <- pow(sigma.phiBR, 2)
    
    mean.pBR ~ dunif(0, 1)
    mu.pBR <- logit(mean.pBR)
    sigma.pBR ~ dunif(0, 10)               # Prior for standard deviation
    tau.pBR <- pow(sigma.pBR, -2)
    sigma2.pBR <- pow(sigma.pBR, 2)
    
    phitBR[1] <- 0
    phitBR[2] <- 0
    phitBR[3] <- phiBR[1]
    phitBR[4] <- 0
    phitBR[5] <- 0
    phitBR[6] <- phiBR[2]
    phitBR[7] <- 0
    phitBR[8] <- 0
    phitBR[9] <- phiBR[3]
    phitBR[10] <- 0
    phitBR[11] <- 0
    phitBR[12] <- phiBR[4]
    phitBR[13] <- 0
    phitBR[14] <- 0
    phitBR[15] <- phiBR[5]
    phitBR[16] <- 0
    phitBR[17] <- 0
    phitBR[18] <- phiBR[6]
    phitBR[19] <- 0
    phitBR[20] <- 0
    phitBR[21] <- phiBR[7]
    phitBR[22] <- 0
    phitBR[23] <- 0
    phitBR[24] <- phiBR[8]
    phitBR[25] <- 0
    phitBR[26] <- 0
    phitBR[27] <- phiBR[9]
    phitBR[28] <- 0
    phitBR[29] <- 0
    
  psightBR[1] <- 0
  psightBR[2] <- pBR[1]
  psightBR[3] <- 0
  psightBR[4] <- 0
  psightBR[5] <- pBR[2]
  psightBR[6] <- 0
  psightBR[7] <- 0
  psightBR[8] <- pBR[3]
  psightBR[9] <- 0
  psightBR[10] <- 0
  psightBR[11] <- pBR[4]
  psightBR[12] <- 0
  psightBR[13] <- 0
  psightBR[14] <- pBR[5]
  psightBR[15] <- 0
  psightBR[16] <- 0
  psightBR[17] <- pBR[6]
  psightBR[18] <- 0
  psightBR[19] <- 0
  psightBR[20] <- pBR[7]
  psightBR[21] <- 0
  psightBR[22] <- 0
  psightBR[23] <- pBR[8]
  psightBR[24] <- 0
  psightBR[25] <- 0
  psightBR[26] <- pBR[9]
  psightBR[27] <- 0
  psightBR[28] <- 0
  psightBR[29] <- pBR[10]

  # Transitions: gamma priors
  for (i in 1:3){
    db[i] ~ dgamma(1, 1)
  }
  for (t in 1:(n.occasions-1)){
    psiDB[1,t] <- 0
    psiDB[2,t] <- (db[1]/sum(db[])) * transDB[t]
    psiDB[3,t] <- (db[2]/sum(db[])) * transDB[t]
    psiDB[4,t] <- (db[3]/sum(db[])) * transDB[t]
    psiDB[5,t] <- 0
  }
  
 for (i in 1:2){
    jb[i] ~ dgamma(1, 1)
 }
 for (t in 1:(n.occasions-1)){
    psiJB[1,t] <- 0
    psiJB[2,t] <- 0
    psiJB[3,t] <- 0
    psiJB[4,t] <- (jb[1]/sum(jb[])) * transJB[t]
    psiJB[5,t] <- (jb[2]/sum(jb[])) * transJB[t]
 }
 
 for (i in 1:2){
    mi[i] ~ dgamma(1, 1)
 }
 for (t in 1:(n.occasions-1)){
    psiMI[1,t] <- 0
    psiMI[2,t] <- 0
    psiMI[3,t] <- 0
    psiMI[4,t] <- (mi[1]/sum(mi[])) * transMI[t]
    psiMI[5,t] <- (mi[2]/sum(mi[])) * transMI[t]
 }
 
 for (i in 1:5){
    se[i] ~ dgamma(1, 1)
 }
 for (t in 1:(n.occasions-1)){
    psiSE[1,t] <- (se[1]/sum(se[])) * transSEDB[t]
    psiSE[2,t] <- (se[2]/sum(se[])) * transSEJB[t]
    psiSE[3,t] <- (se[3]/sum(se[])) * transSEMI[t]
    psiSE[4,t] <- (se[4]/sum(se[])) * transSESE[t]
    psiSE[5,t] <- (se[5]/sum(se[])) * transSEBR[t]
 }
 
 for (i in 1:2){
    br[i] ~ dgamma(1, 1)
 }
 for (t in 1:(n.occasions-1)){
    psiBR[1,t] <- (br[1]/sum(br[])) * transBR[t]
    psiBR[2,t] <- 0
    psiBR[3,t] <- 0
    psiBR[4,t] <- (br[2]/sum(br[])) * transBR[t]
    psiBR[5,t] <- 0
 }
  
  # Define state-transition and reencounter probabilities - note no i index - is no longer individual 
  for (t in 1:(n.occasions-1)){
    psi[1,t,1] <- 0
    psi[1,t,2] <- phitDB[t] * psiDB[2,t]
    psi[1,t,3] <- phitDB[t] * psiDB[3,t]
    psi[1,t,4] <- phitDB[t] * psiDB[4,t]
    psi[1,t,5] <- 0

    psi[2,t,1] <- 0
    psi[2,t,2] <- 0
    psi[2,t,3] <- 0
    psi[2,t,4] <- phitJB[t] * psiJB[4,t]
    psi[2,t,5] <- phitJB[t] * psiJB[5,t]

    psi[3,t,1] <- 0
    psi[3,t,2] <- 0
    psi[3,t,3] <- 0
    psi[3,t,4] <- phitMI[t] * psiMI[4,t]
    psi[3,t,5] <- phitMI[t] * psiMI[5,t]

    psi[4,t,1] <- phiSE[t] * psiSE[1,t]
    psi[4,t,2] <- phiSE[t] * psiSE[2,t]
    psi[4,t,3] <- phiSE[t] * psiSE[3,t]
    psi[4,t,4] <- phiSE[t] * psiSE[4,t]
    psi[4,t,5] <- phiSE[t] * psiSE[5,t]

    psi[5,t,1] <- phitBR[t] * psiBR[1,t]
    psi[5,t,2] <- 0
    psi[5,t,3] <- 0
    psi[5,t,4] <- phitBR[t] * psiBR[4,t]
    psi[5,t,5] <- 0

    po[1,t] <- psightDB[t]
    po[2,t] <- psightJB[t]
    po[3,t] <- psightMI[t]
    po[4,t] <- pSE[t]
    po[5,t] <- psightBR[t]

    
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
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns),
                  season = rep(c(1, 2, 3), length.out = n.occasions-1),
                  transDB = rep(c(1, 0, 0), length.out = n.occasions-1),
                  transJB = rep(c(0, 1, 0), length.out = n.occasions-1),
                  transMI = rep(c(0, 1, 0), length.out = n.occasions-1),
                  transSEDB = rep(c(0, 0, 1), length.out = n.occasions-1),
                  transSEJB = rep(c(1, 0, 0), length.out = n.occasions-1),
                  transSEMI = rep(c(1, 0, 0), length.out = n.occasions-1),
                  transSESE = rep(c(1, 2, 3), length.out = n.occasions-1),
                  transSEBR = rep(c(0, 1, 0), length.out = n.occasions-1),
                  transBR = rep(c(0, 0, 1), length.out = n.occasions-1))

# Initial values 
inits <- function(){list(mean.phiDB = runif(1, 0, 1), mean.phiJB = runif(1, 0, 1), mean.phiMI = runif(1, 0, 1),
                         mean.phiSE = runif(3, 0, 1), mean.phiBR = runif(1, 0, 1),
                         sigma.phiDB = runif(1, 0, 10), sigma.phiJB = runif(1, 0, 10), sigma.phiMI = runif(1, 0, 10),
                         sigma.phiSE = runif(3, 0, 10), sigma.phiBR = runif(1, 0, 10),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1),
                         mean.pSE = runif(3, 0, 1), mean.pBR = runif(1, 0, 1),
                         sigma.pDB = runif(1, 0, 10), sigma.pJB = runif(1, 0, 10), sigma.pMI = runif(1, 0, 10),
                         sigma.pSE = runif(3, 0, 10), sigma.pBR = runif(1, 0, 10))}  


# Parameters monitored
parameters <- c("mean.phiDB", "mean.phiJB", "mean.phiMI", "mean.phiSE", "mean.phiBR",
                "phitDB", "phitJB", "phitMI", "phiSE", "phitBR",
                "mean.pDB", "mean.pJB", "mean.pMI", "mean.pSE", "mean.pBR",
                "psightDB", "psightJB", "psightMI", "pSE", "psightBR",
                "psiDB", "psiJB", "psiMI", "psiSE", "psiBR",
                "fit", "fit.new")

# MCMC settings
ni <- 500000
nt <- 50
nb <- 250000
nc <- 3

# Call JAGS from R (BRT 2 days)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/multistate-seasonal-sim-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/multistate-seasonal-sim-simslist", Sys.Date(), ".rds"))

################################################################################

# plots

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/multistate-seasonal-sim/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          #geom_hline(yintercept = 1, linetype = "dashed") +
          xlab(i))
  
  dev.off()
  
}

# evaluation of fit
mean(sims.list$fit.new > sims.list$fit)

ppcheck <- ggplot(sims.list, aes(x = fit, y = fit.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-multistate-test/reduced-model/multistate-seasonal-sim.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

################################################################################
# Plots

theme_set(theme_bw())

summ <- readRDS("./analysis-output/multistate-seasonal-sim-summary2021-09-06.rds")

sims.list <- readRDS("./analysis-output/multistate-seasonal-sim-simslist2021-09-06.rds")
sims.list <- as.data.frame(sims.list)

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "SE-PreB", "SE-PostB", "SE-NonB",
                             "BR"),
                    value = c(mean.phiDB1, mean.phiJB2, mean.phiMI2, mean.phiSE1,
                              mean.phiSE2, mean.phiSE3, mean.phiBR3))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiDB, mean.phiJB, mean.phiMI, mean.phiSE.1, mean.phiSE.2, mean.phiSE.3, mean.phiBR) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiDB", "DB")) %>%
  mutate(site = str_replace(site, "mean.phiJB", "JB")) %>%
  mutate(site = str_replace(site, "mean.phiMI", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phiSE.1", "SE-PreB")) %>% 
  mutate(site = str_replace(site, "mean.phiSE.2", "SE-PostB")) %>%
  mutate(site = str_replace(site, "mean.phiSE.3", "SE-NonB")) %>%
  mutate(site = str_replace(site, "mean.phiBR", "BR"))

# plot survival
phi.plot <- ggplot() +
  geom_violin(phi.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_boxplot(phi.mod, mapping = aes(x = site, y = estimate, group = site),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phi.mod, mapping = aes(x = site, y = estimate, group = site),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phi.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-multistate-test/reduced-model/phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()
