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

#### SCENARIO 3 ####
# States: DB, JB, MI
# phi: seasonal with temporal random effect
# p: seasonal with temporal random effect
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

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
var.phiDB <- 0.3                       # Temporal variance of survival

mean.pDB1 <- 0.65
mean.pDB2 <- 0.7
mean.pDB3 <- 0.5
mean.pDB4 <- 0.7
var.pDB <- 0.5

# Determine annual survival probabilities
logit.phiDB1 <- rnorm(n.years, qlogis(mean.phiDB1), var.phiDB^0.5)
phiDB1 <- plogis(logit.phiDB1)

logit.phiDB2 <- rnorm(n.years, qlogis(mean.phiDB2), var.phiDB^0.5)
phiDB2 <- plogis(logit.phiDB2)

logit.phiDB3 <- rnorm(n.years, qlogis(mean.phiDB3), var.phiDB^0.5)
phiDB3 <- plogis(logit.phiDB3)

logit.phiDB4 <- rnorm(n.years, qlogis(mean.phiDB4), var.phiDB^0.5)
phiDB4 <- plogis(logit.phiDB4)

phiDB <- c(rbind(phiDB1, phiDB2, phiDB3, phiDB4))
phiDB <- phiDB[1:39]

# Determine annual resighting probabilities
logit.pDB1 <- rnorm(n.years, qlogis(mean.pDB1), var.pDB^0.5)
pDB1 <- plogis(logit.pDB1)

logit.pDB2 <- rnorm(n.years, qlogis(mean.pDB2), var.pDB^0.5)
pDB2 <- plogis(logit.pDB2)

logit.pDB3 <- rnorm(n.years, qlogis(mean.pDB3), var.pDB^0.5)
pDB3 <- plogis(logit.pDB3)

logit.pDB4 <- rnorm(n.years, qlogis(mean.pDB4), var.pDB^0.5)
pDB4 <- plogis(logit.pDB4)

pDB <- c(rbind(pDB1, pDB2, pDB3, pDB4))
pDB <- pDB[1:39]

# JB
mean.phiJB1 <- 0.7
mean.phiJB2 <- 0.9
mean.phiJB3 <- 0.8
mean.phiJB4 <- 0.7
var.phiJB <- 0.3                       # Temporal variance of survival

mean.pJB1 <- 0.4
mean.pJB2 <- 0.6
mean.pJB3 <- 0.3
mean.pJB4 <- 0.7
var.pJB <- 0.5

# Determine annual survival probabilities
logit.phiJB1 <- rnorm(n.years, qlogis(mean.phiJB1), var.phiJB^0.5)
phiJB1 <- plogis(logit.phiJB1)

logit.phiJB2 <- rnorm(n.years, qlogis(mean.phiJB2), var.phiJB^0.5)
phiJB2 <- plogis(logit.phiJB2)

logit.phiJB3 <- rnorm(n.years, qlogis(mean.phiJB3), var.phiJB^0.5)
phiJB3 <- plogis(logit.phiJB3)

logit.phiJB4 <- rnorm(n.years, qlogis(mean.phiJB4), var.phiJB^0.5)
phiJB4 <- plogis(logit.phiJB4)

phiJB <- c(rbind(phiJB1, phiJB2, phiJB3, phiJB4))
phiJB <- phiJB[1:39]

# Determine annual resighting probabilities
logit.pJB1 <- rnorm(n.years, qlogis(mean.pJB1), var.pJB^0.5)
pJB1 <- plogis(logit.pJB1)

logit.pJB2 <- rnorm(n.years, qlogis(mean.pJB2), var.pJB^0.5)
pJB2 <- plogis(logit.pJB2)

logit.pJB3 <- rnorm(n.years, qlogis(mean.pJB3), var.pJB^0.5)
pJB3 <- plogis(logit.pJB3)

logit.pJB4 <- rnorm(n.years, qlogis(mean.pJB4), var.pJB^0.5)
pJB4 <- plogis(logit.pJB4)

pJB <- c(rbind(pJB1, pJB2, pJB3, pJB4))
pJB <- pJB[1:39]

# MI
mean.phiMI1 <- 0.6
mean.phiMI2 <- 0.8
mean.phiMI3 <- 0.7
mean.phiMI4 <- 0.9
var.phiMI <- 0.3                       # Temporal variance of survival

mean.pMI1 <- 0.5
mean.pMI2 <- 0.4
mean.pMI3 <- 0.7
mean.pMI4 <- 0.6
var.pMI <- 0.5

# Determine annual survival probabilities
logit.phiMI1 <- rnorm(n.years, qlogis(mean.phiMI1), var.phiMI^0.5)
phiMI1 <- plogis(logit.phiMI1)

logit.phiMI2 <- rnorm(n.years, qlogis(mean.phiMI2), var.phiMI^0.5)
phiMI2 <- plogis(logit.phiMI2)

logit.phiMI3 <- rnorm(n.years, qlogis(mean.phiMI3), var.phiMI^0.5)
phiMI3 <- plogis(logit.phiMI3)

logit.phiMI4 <- rnorm(n.years, qlogis(mean.phiMI4), var.phiMI^0.5)
phiMI4 <- plogis(logit.phiMI4)

phiMI <- c(rbind(phiMI1, phiMI2, phiMI3, phiMI4))
phiMI <- phiMI[1:39]

# Determine annual resighting probabilities
logit.pMI1 <- rnorm(n.years, qlogis(mean.pMI1), var.pMI^0.5)
pMI1 <- plogis(logit.pMI1)

logit.pMI2 <- rnorm(n.years, qlogis(mean.pMI2), var.pMI^0.5)
pMI2 <- plogis(logit.pMI2)

logit.pMI3 <- rnorm(n.years, qlogis(mean.pMI3), var.pMI^0.5)
pMI3 <- plogis(logit.pMI3)

logit.pMI4 <- rnorm(n.years, qlogis(mean.pMI4), var.pMI^0.5)
pMI4 <- plogis(logit.pMI4)

pMI <- c(rbind(pMI1, pMI2, pMI3, pMI4))
pMI <- pMI[1:39]

# CC
mean.phiCC1 <- 0.75
mean.phiCC2 <- 0.85
mean.phiCC3 <- 0.85
mean.phiCC4 <- 0.8
var.phiCC <- 0.3                       # Temporal variance of survival

mean.pCC1 <- 0.6
mean.pCC2 <- 0.6
mean.pCC3 <- 0.4
mean.pCC4 <- 0.5
var.pCC <- 0.5

# Determine annual survival probabilities
logit.phiCC1 <- rnorm(n.years, qlogis(mean.phiCC1), var.phiCC^0.5)
phiCC1 <- plogis(logit.phiCC1)

logit.phiCC2 <- rnorm(n.years, qlogis(mean.phiCC2), var.phiCC^0.5)
phiCC2 <- plogis(logit.phiCC2)

logit.phiCC3 <- rnorm(n.years, qlogis(mean.phiCC3), var.phiCC^0.5)
phiCC3 <- plogis(logit.phiCC3)

logit.phiCC4 <- rnorm(n.years, qlogis(mean.phiCC4), var.phiCC^0.5)
phiCC4 <- plogis(logit.phiCC4)

phiCC <- c(rbind(phiCC1, phiCC2, phiCC3, phiCC4))
phiCC <- phiCC[1:39]

# Determine annual resighting probabilities
logit.pCC1 <- rnorm(n.years, qlogis(mean.pCC1), var.pCC^0.5)
pCC1 <- plogis(logit.pCC1)

logit.pCC2 <- rnorm(n.years, qlogis(mean.pCC2), var.pCC^0.5)
pCC2 <- plogis(logit.pCC2)

logit.pCC3 <- rnorm(n.years, qlogis(mean.pCC3), var.pCC^0.5)
pCC3 <- plogis(logit.pCC3)

logit.pCC4 <- rnorm(n.years, qlogis(mean.pCC4), var.pCC^0.5)
pCC4 <- plogis(logit.pCC4)

pCC <- c(rbind(pCC1, pCC2, pCC3, pCC4))
pCC <- pCC[1:39]

# SE
mean.phiSE1 <- 0.85
mean.phiSE2 <- 0.9
mean.phiSE3 <- 0.9
mean.phiSE4 <- 0.95
var.phiSE <- 0.3                       # Temporal variance of survival

mean.pSE1 <- 0.6
mean.pSE2 <- 0.6
mean.pSE3 <- 0.4
mean.pSE4 <- 0.3
var.pSE <- 0.5

# Determine annual survival probabilities
logit.phiSE1 <- rnorm(n.years, qlogis(mean.phiSE1), var.phiSE^0.5)
phiSE1 <- plogis(logit.phiSE1)

logit.phiSE2 <- rnorm(n.years, qlogis(mean.phiSE2), var.phiSE^0.5)
phiSE2 <- plogis(logit.phiSE2)

logit.phiSE3 <- rnorm(n.years, qlogis(mean.phiSE3), var.phiSE^0.5)
phiSE3 <- plogis(logit.phiSE3)

logit.phiSE4 <- rnorm(n.years, qlogis(mean.phiSE4), var.phiSE^0.5)
phiSE4 <- plogis(logit.phiSE4)

phiSE <- c(rbind(phiSE1, phiSE2, phiSE3, phiSE4))
phiSE <- phiSE[1:39]

# Determine annual resighting probabilities
logit.pSE1 <- rnorm(n.years, qlogis(mean.pSE1), var.pSE^0.5)
pSE1 <- plogis(logit.pSE1)

logit.pSE2 <- rnorm(n.years, qlogis(mean.pSE2), var.pSE^0.5)
pSE2 <- plogis(logit.pSE2)

logit.pSE3 <- rnorm(n.years, qlogis(mean.pSE3), var.pSE^0.5)
pSE3 <- plogis(logit.pSE3)

logit.pSE4 <- rnorm(n.years, qlogis(mean.pSE4), var.pSE^0.5)
pSE4 <- plogis(logit.pSE4)

pSE <- c(rbind(pSE1, pSE2, pSE3, pSE4))
pSE <- pSE[1:39]

# BR
mean.phiBR1 <- 0.7
mean.phiBR2 <- 0.6
mean.phiBR3 <- 0.85
mean.phiBR4 <- 0.92
var.phiBR <- 0.3                       # Temporal variance of survival

mean.pBR1 <- 0.7
mean.pBR2 <- 0.5
mean.pBR3 <- 0.3
mean.pBR4 <- 0.6
var.pBR <- 0.5

# Determine annual survival probabilities
logit.phiBR1 <- rnorm(n.years, qlogis(mean.phiBR1), var.phiBR^0.5)
phiBR1 <- plogis(logit.phiBR1)

logit.phiBR2 <- rnorm(n.years, qlogis(mean.phiBR2), var.phiBR^0.5)
phiBR2 <- plogis(logit.phiBR2)

logit.phiBR3 <- rnorm(n.years, qlogis(mean.phiBR3), var.phiBR^0.5)
phiBR3 <- plogis(logit.phiBR3)

logit.phiBR4 <- rnorm(n.years, qlogis(mean.phiBR4), var.phiBR^0.5)
phiBR4 <- plogis(logit.phiBR4)

phiBR <- c(rbind(phiBR1, phiBR2, phiBR3, phiBR4))
phiBR <- phiBR[1:39]

# Determine annual resighting probabilities
logit.pBR1 <- rnorm(n.years, qlogis(mean.pBR1), var.pBR^0.5)
pBR1 <- plogis(logit.pBR1)

logit.pBR2 <- rnorm(n.years, qlogis(mean.pBR2), var.pBR^0.5)
pBR2 <- plogis(logit.pBR2)

logit.pBR3 <- rnorm(n.years, qlogis(mean.pBR3), var.pBR^0.5)
pBR3 <- plogis(logit.pBR3)

logit.pBR4 <- rnorm(n.years, qlogis(mean.pBR4), var.pBR^0.5)
pBR4 <- plogis(logit.pBR4)

pBR <- c(rbind(pBR1, pBR2, pBR3, pBR4))
pBR <- pBR[1:39]

# AR
mean.phiAR1 <- 0.8
mean.phiAR2 <- 0.9
mean.phiAR3 <- 0.7
mean.phiAR4 <- 0.93
var.phiAR <- 0.3                       # Temporal variance of survival

mean.pAR1 <- 0.5
mean.pAR2 <- 0.7
mean.pAR3 <- 0.6
mean.pAR4 <- 0.4
var.pAR <- 0.5

# Determine annual survival probabilities
logit.phiAR1 <- rnorm(n.years, qlogis(mean.phiAR1), var.phiAR^0.5)
phiAR1 <- plogis(logit.phiAR1)

logit.phiAR2 <- rnorm(n.years, qlogis(mean.phiAR2), var.phiAR^0.5)
phiAR2 <- plogis(logit.phiAR2)

logit.phiAR3 <- rnorm(n.years, qlogis(mean.phiAR3), var.phiAR^0.5)
phiAR3 <- plogis(logit.phiAR3)

logit.phiAR4 <- rnorm(n.years, qlogis(mean.phiAR4), var.phiAR^0.5)
phiAR4 <- plogis(logit.phiAR4)

phiAR <- c(rbind(phiAR1, phiAR2, phiAR3, phiAR4))
phiAR <- phiAR[1:39]

# Determine annual resighting probabilities
logit.pAR1 <- rnorm(n.years, qlogis(mean.pAR1), var.pAR^0.5)
pAR1 <- plogis(logit.pAR1)

logit.pAR2 <- rnorm(n.years, qlogis(mean.pAR2), var.pAR^0.5)
pAR2 <- plogis(logit.pAR2)

logit.pAR3 <- rnorm(n.years, qlogis(mean.pAR3), var.pAR^0.5)
pAR3 <- plogis(logit.pAR3)

logit.pAR4 <- rnorm(n.years, qlogis(mean.pAR4), var.pAR^0.5)
pAR4 <- plogis(logit.pAR4)

pAR <- c(rbind(pAR1, pAR2, pAR3, pAR4))
pAR <- pAR[1:39]

# transition probabilities
psiDB.DB <- rep(0, n.occasions-1)
psiDB.JB <- rep(0.35, n.occasions-1)
psiDB.MI <- rep(0.30, n.occasions-1)
psiDB.CC <- rep(0.15, n.occasions-1)
psiDB.SE <- rep(0.10, n.occasions-1)
psiDB.BR <- rep(0, n.occasions-1)
psiDB.AR <- rep(0, n.occasions-1)

psiJB.DB <- rep(0, n.occasions-1)
psiJB.JB <- rep(0, n.occasions-1)
psiJB.MI <- rep(0, n.occasions-1)
psiJB.CC <- rep(0.45, n.occasions-1)
psiJB.SE <- rep(0.35, n.occasions-1)
psiJB.BR <- rep(0.20, n.occasions-1)
psiJB.AR <- rep(0, n.occasions-1)

psiMI.DB <- rep(0, n.occasions-1)
psiMI.JB <- rep(0, n.occasions-1)
psiMI.MI <- rep(0, n.occasions-1)
psiMI.CC <- rep(0.35, n.occasions-1)
psiMI.SE <- rep(0.45, n.occasions-1)
psiMI.BR <- rep(0.20, n.occasions-1)
psiMI.AR <- rep(0, n.occasions-1)

psiCC.DB <- rep(0, n.occasions-1)
psiCC.JB <- rep(0, n.occasions-1)
psiCC.MI <- rep(0, n.occasions-1)
psiCC.CC <- rep(0.20, n.occasions-1)
psiCC.SE <- rep(0.20, n.occasions-1)
psiCC.BR <- rep(0.35, n.occasions-1)
psiCC.AR <- rep(0.25, n.occasions-1)

psiSE.DB <- rep(0.20, n.occasions-1)
psiSE.JB <- rep(0.05, n.occasions-1)
psiSE.MI <- rep(0.05, n.occasions-1)
psiSE.CC <- rep(0.10, n.occasions-1)
psiSE.SE <- rep(0.20, n.occasions-1)
psiSE.BR <- rep(0.25, n.occasions-1)
psiSE.AR <- rep(0.15, n.occasions-1)

psiBR.DB <- rep(0.3, n.occasions-1)
psiBR.JB <- rep(0, n.occasions-1)
psiBR.MI <- rep(0, n.occasions-1)
psiBR.CC <- rep(0, n.occasions-1)
psiBR.SE <- rep(0.5, n.occasions-1)
psiBR.BR <- rep(0.1, n.occasions-1)
psiBR.AR <- rep(0.1, n.occasions-1)

psiAR.DB <- rep(0.90, n.occasions-1)
psiAR.JB <- rep(0, n.occasions-1)
psiAR.MI <- rep(0, n.occasions-1)
psiAR.CC <- rep(0, n.occasions-1)
psiAR.SE <- rep(0.10, n.occasions-1)
psiAR.BR <- rep(0, n.occasions-1)
psiAR.AR <- rep(0, n.occasions-1)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(25, 25, 25, 25), length.out = n.occasions) # DB
marked[,2] <- rep(c(10, 10, 10, 10), length.out = n.occasions) # JB
marked[,3] <- rep(c(5, 5, 5, 5), length.out = n.occasions) # MI
marked[,4] <- rep(c(5, 10, 10, 5), length.out = n.occasions) # CC
marked[,5] <- rep(c(5, 5, 5, 5), length.out = n.occasions) # SE
marked[,6] <- rep(c(5, 5, 10, 10), length.out = n.occasions) # BR
marked[,7] <- rep(c(10, 10, 10, 10), length.out = n.occasions) # AR
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
         logit(phiDB[t]) <- mu.phiDB[season[t]] + eps.phiDB[t]
         eps.phiDB[t] ~ dnorm(0, tau.phiDB[season[t]])T(-10,10)

         logit(pDB[t]) <- mu.pDB[season[t]] + eps.pDB[t]
         eps.pDB[t] ~ dnorm(0, tau.pDB[season[t]])T(-10,10)
   }
   
   for (s in 1:4){
         mean.phiDB[s] ~ dunif(0, 1)
         mu.phiDB[s] <- logit(mean.phiDB[s])
         sigma.phiDB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiDB[s] <- pow(sigma.phiDB[s], -2)
         sigma2.phiDB[s] <- pow(sigma.phiDB[s], 2)
      
         mean.pDB[s] ~ dunif(0, 1)
         mu.pDB[s] <- logit(mean.pDB[s])
         sigma.pDB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pDB[s] <- pow(sigma.pDB[s], -2)
         sigma2.pDB[s] <- pow(sigma.pDB[s], 2)
   }
   
   # JB
   for (t in 1:(n.occasions-1)){
         logit(phiJB[t]) <- mu.phiJB[season[t]] + eps.phiJB[t]
         eps.phiJB[t] ~ dnorm(0, tau.phiJB[season[t]])T(-10,10)

         logit(pJB[t]) <- mu.pJB[season[t]] + eps.pJB[t]
         eps.pJB[t] ~ dnorm(0, tau.pJB[season[t]])T(-10,10)
   }

   for (s in 1:4){
         mean.phiJB[s] ~ dunif(0, 1)
         mu.phiJB[s] <- logit(mean.phiJB[s])
         sigma.phiJB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiJB[s] <- pow(sigma.phiJB[s], -2)
         sigma2.phiJB[s] <- pow(sigma.phiJB[s], 2)
      
         mean.pJB[s] ~ dunif(0, 1)
         mu.pJB[s] <- logit(mean.pJB[s])
         sigma.pJB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pJB[s] <- pow(sigma.pJB[s], -2)
         sigma2.pJB[s] <- pow(sigma.pJB[s], 2)
   }

   # MI
   for (t in 1:(n.occasions-1)){
         logit(phiMI[t]) <- mu.phiMI[season[t]] + eps.phiMI[t]
         eps.phiMI[t] ~ dnorm(0, tau.phiMI[season[t]])T(-10,10)

         logit(pMI[t]) <- mu.pMI[season[t]] + eps.pMI[t]
         eps.pMI[t] ~ dnorm(0, tau.pMI[season[t]])T(-10,10)
   }

   for (s in 1:4){
         mean.phiMI[s] ~ dunif(0, 1)
         mu.phiMI[s] <- logit(mean.phiMI[s])
         sigma.phiMI[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiMI[s] <- pow(sigma.phiMI[s], -2)
         sigma2.phiMI[s] <- pow(sigma.phiMI[s], 2)

         mean.pMI[s] ~ dunif(0, 1)
         mu.pMI[s] <- logit(mean.pMI[s])
         sigma.pMI[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pMI[s] <- pow(sigma.pMI[s], -2)
         sigma2.pMI[s] <- pow(sigma.pMI[s], 2)
   }

   # CC
   for (t in 1:(n.occasions-1)){
         logit(phiCC[t]) <- mu.phiCC[season[t]] + eps.phiCC[t]
         eps.phiCC[t] ~ dnorm(0, tau.phiCC[season[t]])T(-10,10)

         logit(pCC[t]) <- mu.pCC[season[t]] + eps.pCC[t]
         eps.pCC[t] ~ dnorm(0, tau.pCC[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiCC[s] ~ dunif(0, 1)
         mu.phiCC[s] <- logit(mean.phiCC[s])
         sigma.phiCC[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiCC[s] <- pow(sigma.phiCC[s], -2)
         sigma2.phiCC[s] <- pow(sigma.phiCC[s], 2)
      
         mean.pCC[s] ~ dunif(0, 1)
         mu.pCC[s] <- logit(mean.pCC[s])
         sigma.pCC[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pCC[s] <- pow(sigma.pCC[s], -2)
         sigma2.pCC[s] <- pow(sigma.pCC[s], 2)
      }

   # SE
   for (t in 1:(n.occasions-1)){
         logit(phiSE[t]) <- mu.phiSE[season[t]] + eps.phiSE[t]
         eps.phiSE[t] ~ dnorm(0, tau.phiSE[season[t]])T(-10,10)
         
         logit(pSE[t]) <- mu.pSE[season[t]] + eps.pSE[t]
         eps.pSE[t] ~ dnorm(0, tau.pSE[season[t]])T(-10,10)
   }

      for (s in 1:4){
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
   for (t in 1:(n.occasions-1)){
         logit(phiBR[t]) <- mu.phiBR[season[t]] + eps.phiBR[t]
         eps.phiBR[t] ~ dnorm(0, tau.phiBR[season[t]])T(-10,10)

         logit(pBR[t]) <- mu.pBR[season[t]] + eps.pBR[t]
         eps.pBR[t] ~ dnorm(0, tau.pBR[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiBR[s] ~ dunif(0, 1)
         mu.phiBR[s] <- logit(mean.phiBR[s])
         sigma.phiBR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiBR[s] <- pow(sigma.phiBR[s], -2)
         sigma2.phiBR[s] <- pow(sigma.phiBR[s], 2)

         mean.pBR[s] ~ dunif(0, 1)
         mu.pBR[s] <- logit(mean.pBR[s])
         sigma.pBR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pBR[s] <- pow(sigma.pBR[s], -2)
         sigma2.pBR[s] <- pow(sigma.pBR[s], 2)
      }

      # AR
      for (t in 1:(n.occasions-1)){
         logit(phiAR[t]) <- mu.phiAR[season[t]] + eps.phiAR[t]
         eps.phiAR[t] ~ dnorm(0, tau.phiAR[season[t]])T(-10,10)

         logit(pAR[t]) <- mu.pAR[season[t]] + eps.pAR[t]
         eps.pAR[t] ~ dnorm(0, tau.pAR[season[t]])T(-10,10)
   }

   for (s in 1:4){
         mean.phiAR[s] ~ dunif(0, 1)
         mu.phiAR[s] <- logit(mean.phiAR[s])
         sigma.phiAR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiAR[s] <- pow(sigma.phiAR[s], -2)
         sigma2.phiAR[s] <- pow(sigma.phiAR[s], 2)

         mean.pAR[s] ~ dunif(0, 1)
         mu.pAR[s] <- logit(mean.pAR[s])
         sigma.pAR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pAR[s] <- pow(sigma.pAR[s], -2)
         sigma2.pAR[s] <- pow(sigma.pAR[s], 2)
   }

# Transitions: gamma priors
   for (i in 1:4){
     db[i] ~ dgamma(1, 1)
   } # i
   
   for (t in 1:(n.occasions-1)){
      psiDB[1,t] <- 0
      psiDB[2,t] <- (db[1]/sum(db[]))*season[t]
      psiDB[3,t] <- (db[2]/sum(db[]))*season[t]
      psiDB[4,t] <- (db[3]/sum(db[]))*season[t]
      psiDB[5,t] <- (db[4]/sum(db[]))*season[t]
      psiDB[6,t] <- 0
      psiDB[7,t] <- 0
      }
   
    for (i in 1:3){
      jb[i] ~ dgamma(1, 1)
    }
    
    for (t in 1:(n.occasions-1)){
      psiJB[1,t] <- 0
      psiJB[2,t] <- 0
      psiJB[3,t] <- 0
      psiJB[4,t] <- (jb[1]/sum(jb[]))*season[t]
      psiJB[5,t] <- (jb[2]/sum(jb[]))*season[t]
      psiJB[6,t] <- (jb[3]/sum(jb[]))*season[t]
      psiJB[7,t] <- 0
    }
    
    for (i in 1:3){
      mi[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){
      psiMI[1,t] <- 0
      psiMI[2,t] <- 0
      psiMI[3,t] <- 0
      psiMI[4,t] <- (mi[1]/sum(mi[]))*season[t]
      psiMI[5,t] <- (mi[2]/sum(mi[]))*season[t]
      psiMI[6,t] <- (mi[3]/sum(mi[]))*season[t]
      psiMI[7,t] <- 0
    }  
    
    for (i in 1:4){
      cc[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){
      psiCC[1,t] <- 0
      psiCC[2,t] <- 0
      psiCC[3,t] <- 0
      psiCC[4,t] <- (cc[1]/sum(cc[]))*season[t]
      psiCC[5,t] <- (cc[2]/sum(cc[]))*season[t]
      psiCC[6,t] <- (cc[3]/sum(cc[]))*season[t]
      psiCC[7,t] <- (cc[4]/sum(cc[]))*season[t]
    }
      
    for (i in 1:7){
      se[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){  
      psiSE[1,t] <- (se[1]/sum(se[]))*season[t]
      psiSE[2,t] <- (se[2]/sum(se[]))*season[t]
      psiSE[3,t] <- (se[3]/sum(se[]))*season[t]
      psiSE[4,t] <- (se[4]/sum(se[]))*season[t]
      psiSE[5,t] <- (se[5]/sum(se[]))*season[t]
      psiSE[6,t] <- (se[6]/sum(se[]))*season[t]
      psiSE[7,t] <- (se[7]/sum(se[]))*season[t]
    }  
    
    for (i in 1:4){
      br[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){
      psiBR[1,t] <- (br[1]/sum(br[]))*season[t]
      psiBR[2,t] <- 0
      psiBR[3,t] <- 0
      psiBR[4,t] <- 0
      psiBR[5,t] <- (br[2]/sum(br[]))*season[t]
      psiBR[6,t] <- (br[3]/sum(br[]))*season[t]
      psiBR[7,t] <- (br[4]/sum(br[]))*season[t]
    }  
    
    for (i in 1:2){
      ar[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){  
      psiAR[1,t] <- (ar[1]/sum(ar[]))*season[t]
      psiAR[2,t] <- 0
      psiAR[3,t] <- 0
      psiAR[4,t] <- 0
      psiAR[5,t] <- (ar[2]/sum(ar[]))*season[t]
      psiAR[6,t] <- 0
      psiAR[7,t] <- 0
  } # t

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- 0
      psi[1,t,2] <- phiDB[t] * psiDB[2,t]
      psi[1,t,3] <- phiDB[t] * psiDB[3,t]
      psi[1,t,4] <- phiDB[t] * psiDB[4,t]
      psi[1,t,5] <- phiDB[t] * psiDB[5,t]
      psi[1,t,6] <- 0
      psi[1,t,7] <- 0
      
      psi[2,t,1] <- 0
      psi[2,t,2] <- 0
      psi[2,t,3] <- 0
      psi[2,t,4] <- phiJB[t] * psiJB[4,t]
      psi[2,t,5] <- phiJB[t] * psiJB[5,t]
      psi[2,t,6] <- phiJB[t] * psiJB[6,t]
      psi[2,t,7] <- 0

      psi[3,t,1] <- 0
      psi[3,t,2] <- 0
      psi[3,t,3] <- 0
      psi[3,t,4] <- phiMI[t] * psiMI[4,t]
      psi[3,t,5] <- phiMI[t] * psiMI[5,t]
      psi[3,t,6] <- phiMI[t] * psiMI[6,t]
      psi[3,t,7] <- 0
      
      psi[4,t,1] <- 0
      psi[4,t,2] <- 0
      psi[4,t,3] <- 0
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
      psi[6,t,2] <- 0
      psi[6,t,3] <- 0
      psi[6,t,4] <- 0
      psi[6,t,5] <- phiBR[t] * psiBR[5,t]
      psi[6,t,6] <- phiBR[t] * psiBR[6,t]
      psi[6,t,7] <- phiBR[t] * psiBR[7,t]
      
      psi[7,t,1] <- phiAR[t] * psiAR[1,t]
      psi[7,t,2] <- 0
      psi[7,t,3] <- 0
      psi[7,t,4] <- 0
      psi[7,t,5] <- phiAR[t] * psiAR[5,t]
      psi[7,t,6] <- 0
      psi[7,t,7] <- 0
      
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
                         sigma.phiDB = runif(4, 0, 10), sigma.phiJB = runif(4, 0, 10), sigma.phiMI = runif(4, 0, 10),
                         sigma.phiCC = runif(4, 0, 10), sigma.phiSE = runif(4, 0, 10), sigma.phiBR = runif(4, 0, 10), sigma.phiAR = runif(4, 0, 10),
                         mean.pDB = runif(4, 0, 1), mean.pJB = runif(4, 0, 1), mean.pMI = runif(4, 0, 1),
                         mean.pCC = runif(4, 0, 1), mean.pSE = runif(4, 0, 1), mean.pBR = runif(4, 0, 1), mean.pAR = runif(4, 0, 1),
                         sigma.pDB = runif(4, 0, 10), sigma.pJB = runif(4, 0, 10), sigma.pMI = runif(4, 0, 10),
                         sigma.pCC = runif(4, 0, 10), sigma.pSE = runif(4, 0, 10), sigma.pBR = runif(4, 0, 10), sigma.pAR = runif(4, 0, 10))}  


# Parameters monitored
parameters <- c("mean.phiDB", "mean.phiJB", "mean.phiMI", "mean.phiCC", "mean.phiSE", "mean.phiBR", "mean.phiAR",
                "mean.pDB", "mean.pJB", "mean.pMI", "mean.pCC", "mean.pSE", "mean.pBR", "mean.pAR",
                "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
                "fit", "fit.new")

# MCMC settings
ni <- 100
nt <- 1
nb <- 70
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/seasonal-multistate-scenario-1-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/seasonal-multistate-scenario-1-simslist", Sys.Date(), ".rds"))
