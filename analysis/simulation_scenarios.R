################################################################################
# SEASONAL SURVIVAL SIMULATION SCENARIOS
#
# gradually increasing model complexity to simulate realistic biological
# scenario
################################################################################

library(rjags)
library(jagsUI)

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

# ms.init.z <- function(ch, f){
#   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
#   states <- max(ch, na.rm = TRUE)
#   known.states <- 1:(states-1)
#   v <- which(ch==states)
#   ch[-v] <- NA
#   ch[v] <- sample(known.states, length(v), replace = TRUE)
#   return(ch)
# }

# modified function for initial z values
mod.init.z <- function(ch, notseen){
  z.known <- ch # get known states
  z.known[z.known==notseen] <- NA
  z.init <- as.data.frame(t(z.known)) # transpose
  z.init <- tidyr::fill(z.init, names(z.init)) # fill all subsequent unknown states with previous known state
  z.init <- t(z.init) # transpose back
  z.init[!is.na(z.known)] <- NA # replace known values with NA
  #z.init[z.init > 1] <- 1
  return(z.init)
}

################################################################################

#### SCENARIO 1 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: all constant
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
n.states <- 6
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93
phiUN <- 0.92

pDB <- 0.6
pJB <- 0.4
pMI <- 0.5
pAR <- 0.5

psiDB.DB <- 0.3
psiDB.JB <- 0.2
psiDB.MI <- 0.2
psiDB.AR <- 0.1

psiJB.DB <- 0.3
psiJB.JB <- 0.3
psiJB.MI <- 0.1
psiJB.AR <- 0.2

psiMI.DB <- 0.1
psiMI.JB <- 0.1
psiMI.MI <- 0.2
psiMI.AR <- 0.5

psiAR.DB <- 0.4
psiAR.JB <- 0.1
psiAR.MI <- 0.1
psiAR.AR <- 0.2

psiUN.DB <- 0.2
psiUN.JB <- 0.2
psiUN.MI <- 0.2
psiUN.AR <- 0.2

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions) # lots banded in DB
marked[,2] <- rep(0, n.occasions) # none banded in JB
marked[,3] <- rep(20, n.occasions) # some banded in MI
marked[,4] <- rep(20, n.occasions) # some banded in AR
marked[,5] <- rep(0, n.occasions) # none banded while in unobserved state
marked[,6] <- rep(0, n.occasions) # none banded while dead

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
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.AR, phiDB*(1-psiDB.DB-psiDB.JB-psiDB.MI-psiDB.AR), 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.AR, phiJB*(1-psiJB.DB-psiJB.JB-psiJB.MI-psiJB.AR), 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.AR, phiMI*(1-psiMI.DB-psiMI.JB-psiMI.MI-psiMI.AR), 1-phiMI,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.AR, phiAR*(1-psiAR.DB-psiAR.JB-psiAR.MI-psiAR.AR), 1-phiAR,
      phiUN*psiUN.DB, phiUN*psiUN.JB, phiUN*psiUN.MI, phiUN*psiUN.AR, phiUN*(1-psiUN.DB-psiUN.JB-psiUN.MI-psiUN.AR), 1-phiUN,
      0,              0,              0,              0,              0,                                             1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pDB, 0,   0,   0,   1-pDB,
      0,   pJB, 0,   0,   1-pJB,
      0,   0,   pMI, 0,   1-pMI,
      0,   0,   0,   pAR, 1-pAR,
      0,   0,   0,   0,   1,   
      0,   0,   0,   0,   1),
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
   phiDB ~ dunif(0, 1)
   phiJB ~ dunif(0, 1)
   phiMI ~ dunif(0, 1)
   phiAR ~ dunif(0, 1)
   phiUN ~ dunif(0, 1)
   
   pDB ~ dunif(0, 1)
   pJB ~ dunif(0, 1)
   pMI ~ dunif(0, 1)
   pAR ~ dunif(0, 1)

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
      for (i in 1:4){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
         psiUN[i] <- exp(lpsiUN[i]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      }
         
    # Calculate the last transition probability
      psiDB[5] <- 1-psiDB[1]-psiDB[2]-psiDB[3]-psiDB[4]
      psiJB[5] <- 1-psiJB[1]-psiJB[2]-psiJB[3]-psiJB[4]
      psiMI[5] <- 1-psiMI[1]-psiMI[2]-psiMI[3]-psiMI[4]
      psiAR[5] <- 1-psiAR[1]-psiAR[2]-psiAR[3]-psiAR[4]
      psiUN[5] <- 1-psiUN[1]-psiUN[2]-psiUN[3]-psiUN[4]

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiDB * psiDB[1]
      ps[1,i,t,2] <- phiDB * psiDB[2]
      ps[1,i,t,3] <- phiDB * psiDB[3]
      ps[1,i,t,4] <- phiDB * psiDB[4]
      ps[1,i,t,5] <- phiDB * psiDB[5]
      ps[1,i,t,6] <- 1-phiDB
      
      ps[2,i,t,1] <- phiJB * psiJB[1]
      ps[2,i,t,2] <- phiJB * psiJB[2]
      ps[2,i,t,3] <- phiJB * psiJB[3]
      ps[2,i,t,4] <- phiJB * psiJB[4]
      ps[2,i,t,5] <- phiJB * psiJB[5]
      ps[2,i,t,6] <- 1-phiJB
      
      ps[3,i,t,1] <- phiMI * psiMI[1]
      ps[3,i,t,2] <- phiMI * psiMI[2]
      ps[3,i,t,3] <- phiMI * psiMI[3]
      ps[3,i,t,4] <- phiMI * psiMI[4]
      ps[3,i,t,5] <- phiMI * psiMI[5]
      ps[3,i,t,6] <- 1-phiMI
      
      ps[4,i,t,1] <- phiAR * psiAR[1]
      ps[4,i,t,2] <- phiAR * psiAR[2]
      ps[4,i,t,3] <- phiAR * psiAR[3]
      ps[4,i,t,4] <- phiAR * psiAR[4]
      ps[4,i,t,5] <- phiAR * psiAR[5]
      ps[4,i,t,6] <- 1-phiAR
      
      ps[5,i,t,1] <- phiUN * psiUN[1]
      ps[5,i,t,2] <- phiUN * psiUN[2]
      ps[5,i,t,3] <- phiUN * psiUN[3]
      ps[5,i,t,4] <- phiUN * psiUN[4]
      ps[5,i,t,5] <- phiUN * psiUN[5]
      ps[5,i,t,6] <- 1-phiUN     
      
      ps[6,i,t,1] <- 0
      ps[6,i,t,2] <- 0
      ps[6,i,t,3] <- 0
      ps[6,i,t,4] <- 0
      ps[6,i,t,5] <- 0
      ps[6,i,t,6] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pDB
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[1,i,t,4] <- 0
      po[1,i,t,5] <- 1-pDB
      
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pJB
      po[2,i,t,3] <- 0
      po[2,i,t,4] <- 0
      po[2,i,t,5] <- 1-pJB
      
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- pMI
      po[3,i,t,4] <- 0
      po[3,i,t,5] <- 1-pMI
      
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 0
      po[4,i,t,4] <- pAR
      po[4,i,t,5] <- 1-pAR
      
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

# Initial values 
inits <- function(){list(phiDB = runif(1, 0, 1), phiJB = runif(1, 0, 1), phiMI = runif(1, 0, 1), phiAR = runif(1, 0, 1), phiUN = runif(1, 0, 1),
                         lpsiDB = rnorm(4), lpsiJB = rnorm(4), lpsiMI = rnorm(4), lpsiAR = rnorm(4), lpsiUN = rnorm(4),
                         pDB = runif(1, 0, 1) , pJB = runif(1, 0, 1), pMI = runif(1, 0, 1), pAR = runif(1, 0, 1),
                         z = mod.init.z(rCH, 5))}  


# Parameters monitored
parameters <- c("phiDB", "phiJB", "phiMI", "phiAR", "phiUN",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "pDB", "pJB", "pMI", "pAR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 2 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: all time-dependent
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
n.states <- 6
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93
phiUN <- 0.92

pDB <- c(0.6, 0.4, 0.5, 0.7, 0.6, 0.6, 0.3, 0.8, 0.7, 0.4, 0.5, 0.6, 0.2, 0.6, 0.8)
pJB <- c(0.4, 0.5, 0.3, 0.6, 0.5, 0.2, 0.4, 0.5, 0.3, 0.7, 0.5, 0.4, 0.6, 0.5, 0.5)
pMI <- c(0.5, 0.6, 0.5, 0.4, 0.3, 0.7, 0.6, 0.4, 0.4, 0.6, 0.7, 0.5, 0.3, 0.2, 0.7)
pAR <- c(0.5, 0.4, 0.7, 0.5, 0.2, 0.6, 0.4, 0.5, 0.7, 0.7, 0.3, 0.2, 0.6, 0.5, 0.6)

psiDB.DB <- 0.3
psiDB.JB <- 0.2
psiDB.MI <- 0.2
psiDB.AR <- 0.1

psiJB.DB <- 0.3
psiJB.JB <- 0.3
psiJB.MI <- 0.1
psiJB.AR <- 0.2

psiMI.DB <- 0.1
psiMI.JB <- 0.1
psiMI.MI <- 0.2
psiMI.AR <- 0.5

psiAR.DB <- 0.4
psiAR.JB <- 0.1
psiAR.MI <- 0.1
psiAR.AR <- 0.2

psiUN.DB <- 0.2
psiUN.JB <- 0.2
psiUN.MI <- 0.2
psiUN.AR <- 0.2

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions) # lots banded in DB
marked[,2] <- rep(0, n.occasions) # none banded in JB
marked[,3] <- rep(20, n.occasions) # some banded in MI
marked[,4] <- rep(20, n.occasions) # some banded in AR
marked[,5] <- rep(0, n.occasions) # none banded while in unobserved state
marked[,6] <- rep(0, n.occasions) # none banded while dead

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.AR, phiDB*(1-psiDB.DB-psiDB.JB-psiDB.MI-psiDB.AR), 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.AR, phiJB*(1-psiJB.DB-psiJB.JB-psiJB.MI-psiJB.AR), 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.AR, phiMI*(1-psiMI.DB-psiMI.JB-psiMI.MI-psiMI.AR), 1-phiMI,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.AR, phiAR*(1-psiAR.DB-psiAR.JB-psiAR.MI-psiAR.AR), 1-phiAR,
      phiUN*psiUN.DB, phiUN*psiUN.JB, phiUN*psiUN.MI, phiUN*psiUN.AR, phiUN*(1-psiUN.DB-psiUN.JB-psiUN.MI-psiUN.AR), 1-phiUN,
      0,              0,              0,              0,              0,                                             1),
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
   
   pDB[t] ~ dunif(0, 1)
   pJB[t] ~ dunif(0, 1)
   pMI[t] ~ dunif(0, 1)
   pAR[t] ~ dunif(0, 1)
   }
   
   for (u in 1:5){
   mean.phi[u] ~ dunif(0, 1)
   }

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
      for (i in 1:4){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
         psiUN[i] <- exp(lpsiUN[i]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      }
         
    # Calculate the last transition probability
      psiDB[5] <- 1-psiDB[1]-psiDB[2]-psiDB[3]-psiDB[4]
      psiJB[5] <- 1-psiJB[1]-psiJB[2]-psiJB[3]-psiJB[4]
      psiMI[5] <- 1-psiMI[1]-psiMI[2]-psiMI[3]-psiMI[4]
      psiAR[5] <- 1-psiAR[1]-psiAR[2]-psiAR[3]-psiAR[4]
      psiUN[5] <- 1-psiUN[1]-psiUN[2]-psiUN[3]-psiUN[4]

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiDB[t] * psiDB[1]
      ps[1,i,t,2] <- phiDB[t] * psiDB[2]
      ps[1,i,t,3] <- phiDB[t] * psiDB[3]
      ps[1,i,t,4] <- phiDB[t] * psiDB[4]
      ps[1,i,t,5] <- phiDB[t] * psiDB[5]
      ps[1,i,t,6] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- phiJB[t] * psiJB[1]
      ps[2,i,t,2] <- phiJB[t] * psiJB[2]
      ps[2,i,t,3] <- phiJB[t] * psiJB[3]
      ps[2,i,t,4] <- phiJB[t] * psiJB[4]
      ps[2,i,t,5] <- phiJB[t] * psiJB[5]
      ps[2,i,t,6] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- phiMI[t] * psiMI[1]
      ps[3,i,t,2] <- phiMI[t] * psiMI[2]
      ps[3,i,t,3] <- phiMI[t] * psiMI[3]
      ps[3,i,t,4] <- phiMI[t] * psiMI[4]
      ps[3,i,t,5] <- phiMI[t] * psiMI[5]
      ps[3,i,t,6] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[1]
      ps[4,i,t,2] <- phiAR[t] * psiAR[2]
      ps[4,i,t,3] <- phiAR[t] * psiAR[3]
      ps[4,i,t,4] <- phiAR[t] * psiAR[4]
      ps[4,i,t,5] <- phiAR[t] * psiAR[5]
      ps[4,i,t,6] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- phiUN[t] * psiUN[1]
      ps[5,i,t,2] <- phiUN[t] * psiUN[2]
      ps[5,i,t,3] <- phiUN[t] * psiUN[3]
      ps[5,i,t,4] <- phiUN[t] * psiUN[4]
      ps[5,i,t,5] <- phiUN[t] * psiUN[5]
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

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = rnorm(4), lpsiJB = rnorm(4), lpsiMI = rnorm(4), lpsiAR = rnorm(4), lpsiUN = rnorm(4),
                         pDB = runif(n.occasions-1, 0, 1) , pJB = runif(n.occasions-1, 0, 1), pMI = runif(n.occasions-1, 0, 1), pAR = runif(n.occasions-1, 0, 1),
                         z = mod.init.z(rCH, 5))}  


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "pDB", "pJB", "pMI", "pAR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 3 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: temporal random effect
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
n.states <- 6
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93
phiUN <- 0.92

pDB <- c(0.6, 0.4, 0.5, 0.7, 0.6, 0.6, 0.3, 0.8, 0.7, 0.4, 0.5, 0.6, 0.2, 0.6, 0.8)
pJB <- c(0.4, 0.5, 0.3, 0.6, 0.5, 0.2, 0.4, 0.5, 0.3, 0.7, 0.5, 0.4, 0.6, 0.5, 0.5)
pMI <- c(0.5, 0.6, 0.5, 0.4, 0.3, 0.7, 0.6, 0.4, 0.4, 0.6, 0.7, 0.5, 0.3, 0.2, 0.7)
pAR <- c(0.5, 0.4, 0.7, 0.5, 0.2, 0.6, 0.4, 0.5, 0.7, 0.7, 0.3, 0.2, 0.6, 0.5, 0.6)

psiDB.DB <- 0.3
psiDB.JB <- 0.2
psiDB.MI <- 0.2
psiDB.AR <- 0.1

psiJB.DB <- 0.3
psiJB.JB <- 0.3
psiJB.MI <- 0.1
psiJB.AR <- 0.2

psiMI.DB <- 0.1
psiMI.JB <- 0.1
psiMI.MI <- 0.2
psiMI.AR <- 0.5

psiAR.DB <- 0.4
psiAR.JB <- 0.1
psiAR.MI <- 0.1
psiAR.AR <- 0.2

psiUN.DB <- 0.2
psiUN.JB <- 0.2
psiUN.MI <- 0.2
psiUN.AR <- 0.2

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions) # lots banded in DB
marked[,2] <- rep(0, n.occasions) # none banded in JB
marked[,3] <- rep(20, n.occasions) # some banded in MI
marked[,4] <- rep(20, n.occasions) # some banded in AR
marked[,5] <- rep(0, n.occasions) # none banded while in unobserved state
marked[,6] <- rep(0, n.occasions) # none banded while dead

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.AR, phiDB*(1-psiDB.DB-psiDB.JB-psiDB.MI-psiDB.AR), 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.AR, phiJB*(1-psiJB.DB-psiJB.JB-psiJB.MI-psiJB.AR), 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.AR, phiMI*(1-psiMI.DB-psiMI.JB-psiMI.MI-psiMI.AR), 1-phiMI,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.AR, phiAR*(1-psiAR.DB-psiAR.JB-psiAR.MI-psiAR.AR), 1-phiAR,
      phiUN*psiUN.DB, phiUN*psiUN.JB, phiUN*psiUN.MI, phiUN*psiUN.AR, phiUN*(1-psiUN.DB-psiUN.JB-psiUN.MI-psiUN.AR), 1-phiUN,
      0,              0,              0,              0,              0,                                             1),
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
   
   logit(pDB[t]) <- muDB + epsilonDB[t]
   epsilonDB[t] ~ dnorm(0, tauDB)T(-15, 15)
   logit(pJB[t]) <- muJB + epsilonJB[t]
   epsilonJB[t] ~ dnorm(0, tauJB)T(-15, 15)
   logit(pMI[t]) <- muMI + epsilonMI[t]
   epsilonMI[t] ~ dnorm(0, tauMI)T(-15, 15)
   logit(pAR[t]) <- muAR + epsilonAR[t]
   epsilonAR[t] ~ dnorm(0, tauAR)T(-15, 15)
   }
   
   for (u in 1:5){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   muDB <- log(mean.pDB / (1-mean.pDB))     # Logit transformation
   mean.pDB ~ dunif(0, 1)                   # Prior for mean survival
   tauDB <- pow(sigmaDB, -2)
   sigmaDB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2DB <- pow(sigmaDB, 2)              # Temporal variance
   
   muJB <- log(mean.pJB / (1-mean.pJB))     # Logit transformation
   mean.pJB ~ dunif(0, 1)                   # Prior for mean survival
   tauJB <- pow(sigmaJB, -2)
   sigmaJB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2JB <- pow(sigmaJB, 2)              # Temporal variance
   
   muMI <- log(mean.pMI / (1-mean.pMI))     # Logit transformation
   mean.pMI ~ dunif(0, 1)                   # Prior for mean survival
   tauMI <- pow(sigmaMI, -2)
   sigmaMI ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2MI <- pow(sigmaMI, 2)              # Temporal variance
   
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
      for (i in 1:4){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
         psiUN[i] <- exp(lpsiUN[i]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      }
         
    # Calculate the last transition probability
      psiDB[5] <- 1-psiDB[1]-psiDB[2]-psiDB[3]-psiDB[4]
      psiJB[5] <- 1-psiJB[1]-psiJB[2]-psiJB[3]-psiJB[4]
      psiMI[5] <- 1-psiMI[1]-psiMI[2]-psiMI[3]-psiMI[4]
      psiAR[5] <- 1-psiAR[1]-psiAR[2]-psiAR[3]-psiAR[4]
      psiUN[5] <- 1-psiUN[1]-psiUN[2]-psiUN[3]-psiUN[4]

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiDB[t] * psiDB[1]
      ps[1,i,t,2] <- phiDB[t] * psiDB[2]
      ps[1,i,t,3] <- phiDB[t] * psiDB[3]
      ps[1,i,t,4] <- phiDB[t] * psiDB[4]
      ps[1,i,t,5] <- phiDB[t] * psiDB[5]
      ps[1,i,t,6] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- phiJB[t] * psiJB[1]
      ps[2,i,t,2] <- phiJB[t] * psiJB[2]
      ps[2,i,t,3] <- phiJB[t] * psiJB[3]
      ps[2,i,t,4] <- phiJB[t] * psiJB[4]
      ps[2,i,t,5] <- phiJB[t] * psiJB[5]
      ps[2,i,t,6] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- phiMI[t] * psiMI[1]
      ps[3,i,t,2] <- phiMI[t] * psiMI[2]
      ps[3,i,t,3] <- phiMI[t] * psiMI[3]
      ps[3,i,t,4] <- phiMI[t] * psiMI[4]
      ps[3,i,t,5] <- phiMI[t] * psiMI[5]
      ps[3,i,t,6] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[1]
      ps[4,i,t,2] <- phiAR[t] * psiAR[2]
      ps[4,i,t,3] <- phiAR[t] * psiAR[3]
      ps[4,i,t,4] <- phiAR[t] * psiAR[4]
      ps[4,i,t,5] <- phiAR[t] * psiAR[5]
      ps[4,i,t,6] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- phiUN[t] * psiUN[1]
      ps[5,i,t,2] <- phiUN[t] * psiUN[2]
      ps[5,i,t,3] <- phiUN[t] * psiUN[3]
      ps[5,i,t,4] <- phiUN[t] * psiUN[4]
      ps[5,i,t,5] <- phiUN[t] * psiUN[5]
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

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = rnorm(4), lpsiJB = rnorm(4), lpsiMI = rnorm(4), lpsiAR = rnorm(4), lpsiUN = rnorm(4),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 1), sigmaJB = runif(1, 0, 1), sigmaMI = runif(1, 0, 1), sigmaAR = runif(1, 0, 1),
                         z = mod.init.z(rCH, 5))} 


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "pDB", "mean.pDB", "sigma2DB",
                "pJB", "mean.pJB", "sigma2JB",
                "pMI", "mean.pMI", "sigma2MI",
                "pAR", "mean.pAR", "sigma2AR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 4 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: temporal random effect, some JB periods set to 0
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
n.states <- 6
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93
phiUN <- 0.92

pDB <- c(0.6, 0.4, 0.5, 0.7, 0.6, 0.6, 0.3, 0.8, 0.7, 0.4, 0.5, 0.6, 0.2, 0.6)
pJB <- c(0, 0, 0.3, 0.6, 0.5, 0.2, 0.4, 0.5, 0.3, 0.7, 0.5, 0.4, 0.6, 0.5)
pMI <- c(0.5, 0.6, 0.5, 0.4, 0.3, 0.7, 0.6, 0.4, 0.4, 0.6, 0.7, 0.5, 0.3, 0.2)
pAR <- c(0.5, 0.4, 0.7, 0.5, 0.2, 0.6, 0.4, 0.5, 0.7, 0.7, 0.3, 0.2, 0.6, 0.5)

psiDB.DB <- 0.3
psiDB.JB <- 0.2
psiDB.MI <- 0.2
psiDB.AR <- 0.1

psiJB.DB <- 0.3
psiJB.JB <- 0.3
psiJB.MI <- 0.1
psiJB.AR <- 0.2

psiMI.DB <- 0.1
psiMI.JB <- 0.1
psiMI.MI <- 0.2
psiMI.AR <- 0.5

psiAR.DB <- 0.4
psiAR.JB <- 0.1
psiAR.MI <- 0.1
psiAR.AR <- 0.2

psiUN.DB <- 0.2
psiUN.JB <- 0.2
psiUN.MI <- 0.2
psiUN.AR <- 0.2

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions) # lots banded in DB
marked[,2] <- rep(0, n.occasions) # none banded in JB
marked[,3] <- rep(20, n.occasions) # some banded in MI
marked[,4] <- rep(20, n.occasions) # some banded in AR
marked[,5] <- rep(0, n.occasions) # none banded while in unobserved state
marked[,6] <- rep(0, n.occasions) # none banded while dead

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.AR, phiDB*(1-psiDB.DB-psiDB.JB-psiDB.MI-psiDB.AR), 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.AR, phiJB*(1-psiJB.DB-psiJB.JB-psiJB.MI-psiJB.AR), 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.AR, phiMI*(1-psiMI.DB-psiMI.JB-psiMI.MI-psiMI.AR), 1-phiMI,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.AR, phiAR*(1-psiAR.DB-psiAR.JB-psiAR.MI-psiAR.AR), 1-phiAR,
      phiUN*psiUN.DB, phiUN*psiUN.JB, phiUN*psiUN.MI, phiUN*psiUN.AR, phiUN*(1-psiUN.DB-psiUN.JB-psiUN.MI-psiUN.AR), 1-phiUN,
      0,              0,              0,              0,              0,                                             1),
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
   
   logit(pDB[t]) <- muDB + epsilonDB[t]
   epsilonDB[t] ~ dnorm(0, tauDB)T(-15, 15)
   logit(pMI[t]) <- muMI + epsilonMI[t]
   epsilonMI[t] ~ dnorm(0, tauMI)T(-15, 15)
   logit(pAR[t]) <- muAR + epsilonAR[t]
   epsilonAR[t] ~ dnorm(0, tauAR)T(-15, 15)
   }
   
   pJB[1] <- 0
   pJB[2] <- 0
   
   for (t in 3:(n.occasions-1)){
   logit(pJB[t]) <- muJB + epsilonJB[t]
   epsilonJB[t] ~ dnorm(0, tauJB)T(-15, 15)
   }
   
   for (u in 1:5){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   muDB <- log(mean.pDB / (1-mean.pDB))     # Logit transformation
   mean.pDB ~ dunif(0, 1)                   # Prior for mean survival
   tauDB <- pow(sigmaDB, -2)
   sigmaDB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2DB <- pow(sigmaDB, 2)              # Temporal variance
   
   muJB <- log(mean.pJB / (1-mean.pJB))     # Logit transformation
   mean.pJB ~ dunif(0, 1)                   # Prior for mean survival
   tauJB <- pow(sigmaJB, -2)
   sigmaJB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2JB <- pow(sigmaJB, 2)              # Temporal variance
   
   muMI <- log(mean.pMI / (1-mean.pMI))     # Logit transformation
   mean.pMI ~ dunif(0, 1)                   # Prior for mean survival
   tauMI <- pow(sigmaMI, -2)
   sigmaMI ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2MI <- pow(sigmaMI, 2)              # Temporal variance
   
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
      for (i in 1:4){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
         psiUN[i] <- exp(lpsiUN[i]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      }
         
    # Calculate the last transition probability
      psiDB[5] <- 1-psiDB[1]-psiDB[2]-psiDB[3]-psiDB[4]
      psiJB[5] <- 1-psiJB[1]-psiJB[2]-psiJB[3]-psiJB[4]
      psiMI[5] <- 1-psiMI[1]-psiMI[2]-psiMI[3]-psiMI[4]
      psiAR[5] <- 1-psiAR[1]-psiAR[2]-psiAR[3]-psiAR[4]
      psiUN[5] <- 1-psiUN[1]-psiUN[2]-psiUN[3]-psiUN[4]

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiDB[t] * psiDB[1]
      ps[1,i,t,2] <- phiDB[t] * psiDB[2]
      ps[1,i,t,3] <- phiDB[t] * psiDB[3]
      ps[1,i,t,4] <- phiDB[t] * psiDB[4]
      ps[1,i,t,5] <- phiDB[t] * psiDB[5]
      ps[1,i,t,6] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- phiJB[t] * psiJB[1]
      ps[2,i,t,2] <- phiJB[t] * psiJB[2]
      ps[2,i,t,3] <- phiJB[t] * psiJB[3]
      ps[2,i,t,4] <- phiJB[t] * psiJB[4]
      ps[2,i,t,5] <- phiJB[t] * psiJB[5]
      ps[2,i,t,6] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- phiMI[t] * psiMI[1]
      ps[3,i,t,2] <- phiMI[t] * psiMI[2]
      ps[3,i,t,3] <- phiMI[t] * psiMI[3]
      ps[3,i,t,4] <- phiMI[t] * psiMI[4]
      ps[3,i,t,5] <- phiMI[t] * psiMI[5]
      ps[3,i,t,6] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[1]
      ps[4,i,t,2] <- phiAR[t] * psiAR[2]
      ps[4,i,t,3] <- phiAR[t] * psiAR[3]
      ps[4,i,t,4] <- phiAR[t] * psiAR[4]
      ps[4,i,t,5] <- phiAR[t] * psiAR[5]
      ps[4,i,t,6] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- phiUN[t] * psiUN[1]
      ps[5,i,t,2] <- phiUN[t] * psiUN[2]
      ps[5,i,t,3] <- phiUN[t] * psiUN[3]
      ps[5,i,t,4] <- phiUN[t] * psiUN[4]
      ps[5,i,t,5] <- phiUN[t] * psiUN[5]
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

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = rnorm(4), lpsiJB = rnorm(4), lpsiMI = rnorm(4), lpsiAR = rnorm(4), lpsiUN = rnorm(4),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 1), sigmaJB = runif(1, 0, 1), sigmaMI = runif(1, 0, 1), sigmaAR = runif(1, 0, 1),
                         z = mod.init.z(rCH, 5))} 


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "pDB", "mean.pDB", "sigma2DB",
                "pJB", "mean.pJB", "sigma2JB",
                "pMI", "mean.pMI", "sigma2MI",
                "pAR", "mean.pAR", "sigma2AR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 5 ####

# States: DB, JB, MI, AR
# phi: all constant
# p: temporal random effect, some JB periods set to 0
# psi: following annual cycle

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
n.states <- 5
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93

pDB <- c(0.6, 0.4, 0.5, 0.7, 0.6, 0.6, 0.3, 0.8, 0.7, 0.4, 0.5, 0.6, 0.2, 0.6)
pJB <- c(0, 0, 0.3, 0.6, 0.5, 0.2, 0.4, 0.5, 0.3, 0.7, 0.5, 0.4, 0.6, 0.5)
pMI <- c(0.5, 0.6, 0.5, 0.4, 0.3, 0.7, 0.6, 0.4, 0.4, 0.6, 0.7, 0.5, 0.3, 0.2)
pAR <- c(0.5, 0.4, 0.7, 0.5, 0.2, 0.6, 0.4, 0.5, 0.7, 0.7, 0.3, 0.2, 0.6, 0.5)

psiDB.DB <- 0
psiDB.JB <- 0.5
psiDB.MI <- 0.5
psiDB.AR <- 0

psiJB.DB <- 0
psiJB.JB <- 0
psiJB.MI <- 0
psiJB.AR <- 1

psiMI.DB <- 0
psiMI.JB <- 0
psiMI.MI <- 0
psiMI.AR <- 1

psiAR.DB <- 1
psiAR.JB <- 0
psiAR.MI <- 0
psiAR.AR <- 0

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(100, 0, 0), n.occasions/3)
marked[,2] <- rep(0, n.occasions)
marked[,3] <- rep(c(0, 20, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 20), n.occasions/3)
marked[,5] <- rep(0, n.occasions)

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
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*(1-psiDB.DB-psiDB.JB-psiDB.MI), 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*(1-psiJB.DB-psiJB.JB-psiJB.MI), 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*(1-psiMI.DB-psiMI.JB-psiMI.MI), 1-phiMI,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*(1-psiAR.DB-psiAR.JB-psiAR.MI), 1-phiAR,
      0,              0,              0,              0,                                    1),
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
      0,      0,      0,      0,      1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = NA)
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

   logit(pDB[t]) <- muDB + epsilonDB[t]
   epsilonDB[t] ~ dnorm(0, tauDB)T(-15, 15)
   logit(pMI[t]) <- muMI + epsilonMI[t]
   epsilonMI[t] ~ dnorm(0, tauMI)T(-15, 15)
   logit(pAR[t]) <- muAR + epsilonAR[t]
   epsilonAR[t] ~ dnorm(0, tauAR)T(-15, 15)
   }
   
   pJB[1] <- 0
   pJB[2] <- 0
   
   for (t in 3:(n.occasions-1)){
   logit(pJB[t]) <- muJB + epsilonJB[t]
   epsilonJB[t] ~ dnorm(0, tauJB)T(-15, 15)
   }
   
   for (u in 1:4){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   muDB <- log(mean.pDB / (1-mean.pDB))     # Logit transformation
   mean.pDB ~ dunif(0, 1)                   # Prior for mean survival
   tauDB <- pow(sigmaDB, -2)
   sigmaDB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2DB <- pow(sigmaDB, 2)              # Temporal variance
   
   muJB <- log(mean.pJB / (1-mean.pJB))     # Logit transformation
   mean.pJB ~ dunif(0, 1)                   # Prior for mean survival
   tauJB <- pow(sigmaJB, -2)
   sigmaJB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2JB <- pow(sigmaJB, 2)              # Temporal variance
   
   muMI <- log(mean.pMI / (1-mean.pMI))     # Logit transformation
   mean.pMI ~ dunif(0, 1)                   # Prior for mean survival
   tauMI <- pow(sigmaMI, -2)
   sigmaMI ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2MI <- pow(sigmaMI, 2)              # Temporal variance
   
   muAR <- log(mean.pAR / (1-mean.pAR))     # Logit transformation
   mean.pAR ~ dunif(0, 1)                   # Prior for mean survival
   tauAR <- pow(sigmaAR, -2)
   sigmaAR ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2AR <- pow(sigmaAR, 2)              # Temporal variance

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probs
      for (i in 1:3){
         lpsiDB[i] ~ dnorm(0, 0.001)
         lpsiJB[i] ~ dnorm(0, 0.001)
         lpsiMI[i] ~ dnorm(0, 0.001)
         lpsiAR[i] ~ dnorm(0, 0.001)
      }
         
      # Constrain the transitions such that their sum is < 1
      for (i in 1:3){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]))
      }
         
    # Calculate the last transition probability
      psiDB[4] <- 1-psiDB[1]-psiDB[2]-psiDB[3]
      psiJB[4] <- 1-psiJB[1]-psiJB[2]-psiJB[3]
      psiMI[4] <- 1-psiMI[1]-psiMI[2]-psiMI[3]
      psiAR[4] <- 1-psiAR[1]-psiAR[2]-psiAR[3]

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phiDB[t] * psiDB[2]
      ps[1,i,t,3] <- phiDB[t] * psiDB[3]
      ps[1,i,t,4] <- 0
      ps[1,i,t,5] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 0
      ps[2,i,t,4] <- phiJB[t] * psiJB[4]
      ps[2,i,t,5] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phiMI[t] * psiMI[4]
      ps[3,i,t,5] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[1]
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- 0
      ps[5,i,t,2] <- 0
      ps[5,i,t,3] <- 0
      ps[5,i,t,4] <- 0
      ps[5,i,t,5] <- 1

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

# modify function for initial z values again to follow annual cycle
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

# Initial values 
inits <- function(){list(mean.phi = runif(4, 0, 1),
                         lpsiDB = rnorm(3), lpsiJB = rnorm(3), lpsiMI = rnorm(3), lpsiAR = rnorm(3),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 1), sigmaJB = runif(1, 0, 1), sigmaMI = runif(1, 0, 1), sigmaAR = runif(1, 0, 1),
                         z = mod.init.z2(rCH, 5))} 


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR",
                "pDB", "mean.pDB", "sigma2DB",
                "pJB", "mean.pJB", "sigma2JB",
                "pMI", "mean.pMI", "sigma2MI",
                "pAR", "mean.pAR", "sigma2AR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 6 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: temporal random effect, some JB periods set to 0
# psi: following annual cycle (but constant)

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
n.states <- 6
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93
phiUN <- 0.92

pDB <- c(0.6, 0.4, 0.5, 0.7, 0.6, 0.6, 0.3, 0.8, 0.7, 0.4, 0.5, 0.6, 0.2, 0.6)
pJB <- c(0, 0, 0.3, 0.6, 0.5, 0.2, 0.4, 0.5, 0.3, 0.7, 0.5, 0.4, 0.6, 0.5)
pMI <- c(0.5, 0.6, 0.5, 0.4, 0.3, 0.7, 0.6, 0.4, 0.4, 0.6, 0.7, 0.5, 0.3, 0.2)
pAR <- c(0.5, 0.4, 0.7, 0.5, 0.2, 0.6, 0.4, 0.5, 0.7, 0.7, 0.3, 0.2, 0.6, 0.5)

psiDB.DB <- 0
psiDB.JB <- 0.3
psiDB.MI <- 0.3
psiDB.AR <- 0

psiJB.DB <- 0
psiJB.JB <- 0
psiJB.MI <- 0
psiJB.AR <- 0.7

psiMI.DB <- 0
psiMI.JB <- 0
psiMI.MI <- 0
psiMI.AR <- 0.6

psiAR.DB <- 0.8
psiAR.JB <- 0
psiAR.MI <- 0
psiAR.AR <- 0

psiUN.DB <- rep(c(0, 0, 0.8), n.occasions/3)
psiUN.JB <- rep(c(0.4, 0, 0), n.occasions/3)
psiUN.MI <- rep(c(0.4, 0, 0), n.occasions/3)
psiUN.AR <- rep(c(0, 0.7, 0), n.occasions/3)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(100, 0, 0), n.occasions/3)
marked[,2] <- rep(0, n.occasions)
marked[,3] <- rep(c(0, 20, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 20), n.occasions/3)
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
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.AR, phiDB*(1-psiDB.JB-psiDB.MI),                   1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.AR, phiJB*(1-psiJB.AR),                            1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.AR, phiMI*(1-psiMI.AR),                            1-phiMI,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.AR, phiAR*(1-psiAR.DB),                            1-phiAR,
      phiUN*psiUN.DB[t], phiUN*psiUN.JB[t], phiUN*psiUN.MI[t], phiUN*psiUN.AR[t], phiUN*(1-psiUN.DB[t]-psiUN.JB[t]-psiUN.MI[t]-psiUN.AR[t]), 1-phiUN,
      0,              0,              0,              0,              0,                                             1),
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

   logit(pDB[t]) <- muDB + epsilonDB[t]
   epsilonDB[t] ~ dnorm(0, tauDB)T(-15, 15)
   logit(pMI[t]) <- muMI + epsilonMI[t]
   epsilonMI[t] ~ dnorm(0, tauMI)T(-15, 15)
   logit(pAR[t]) <- muAR + epsilonAR[t]
   epsilonAR[t] ~ dnorm(0, tauAR)T(-15, 15)
   }
   
   pJB[1] <- 0
   pJB[2] <- 0
   
   for (t in 3:(n.occasions-1)){
   logit(pJB[t]) <- muJB + epsilonJB[t]
   epsilonJB[t] ~ dnorm(0, tauJB)T(-15, 15)
   }
   
   for (u in 1:5){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   muDB <- log(mean.pDB / (1-mean.pDB))     # Logit transformation
   mean.pDB ~ dunif(0, 1)                   # Prior for mean survival
   tauDB <- pow(sigmaDB, -2)
   sigmaDB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2DB <- pow(sigmaDB, 2)              # Temporal variance
   
   muJB <- log(mean.pJB / (1-mean.pJB))     # Logit transformation
   mean.pJB ~ dunif(0, 1)                   # Prior for mean survival
   tauJB <- pow(sigmaJB, -2)
   sigmaJB ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2JB <- pow(sigmaJB, 2)              # Temporal variance
   
   muMI <- log(mean.pMI / (1-mean.pMI))     # Logit transformation
   mean.pMI ~ dunif(0, 1)                   # Prior for mean survival
   tauMI <- pow(sigmaMI, -2)
   sigmaMI ~ dunif(0, 10)                   # Prior on standard deviation
   sigma2MI <- pow(sigmaMI, 2)              # Temporal variance
   
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
      for (i in 1:4){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
         psiUN[i] <- exp(lpsiUN[i]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      }
         
    # Calculate the last transition probability
      psiDB[5] <- 1-psiDB[1]-psiDB[2]-psiDB[3]-psiDB[4]
      psiJB[5] <- 1-psiJB[1]-psiJB[2]-psiJB[3]-psiJB[4]
      psiMI[5] <- 1-psiMI[1]-psiMI[2]-psiMI[3]-psiMI[4]
      psiAR[5] <- 1-psiAR[1]-psiAR[2]-psiAR[3]-psiAR[4]
      psiUN[5] <- 1-psiUN[1]-psiUN[2]-psiUN[3]-psiUN[4]
      
    #   # Normal priors on logit of all but one transition probs
    #   for (i in 1:4){
    #   for (t in 1:(n.occasions-1)){
    #      lpsiUN[i,t] ~ dnorm(0, 0.001)
    #   }
    #   }
    #      
    #   # Constrain the transitions such that their sum is < 1
    #   for (i in 1:4){
    #   for (t in 1:(n.occasions-1)){
    #      psiUN[i,t] <- exp(lpsiUN[i,t]) / (1 + exp(lpsiUN[1,t]) + exp(lpsiUN[2,t]) + exp(lpsiUN[3,t]) + exp(lpsiUN[4,t]))
    #   }
    #   }
    #      
    # # Calculate the last transition probability
    # for (t in 1:(n.occasions-1)){
    #   psiUN[5,t] <- 1-psiUN[1,t]-psiUN[2,t]-psiUN[3,t]-psiUN[4,t]
    # }

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phiDB[t] * psiDB[2]
      ps[1,i,t,3] <- phiDB[t] * psiDB[3]
      ps[1,i,t,4] <- 0
      ps[1,i,t,5] <- phiDB[t] * psiDB[5]
      ps[1,i,t,6] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 0
      ps[2,i,t,4] <- phiJB[t] * psiJB[4]
      ps[2,i,t,5] <- phiJB[t] * psiJB[5]
      ps[2,i,t,6] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phiMI[t] * psiMI[4]
      ps[3,i,t,5] <- phiMI[t] * psiMI[5]
      ps[3,i,t,6] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[1]
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- phiAR[t] * psiAR[5]
      ps[4,i,t,6] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- phiUN[t] * psiUN[1]
      ps[5,i,t,2] <- phiUN[t] * psiUN[2]
      ps[5,i,t,3] <- phiUN[t] * psiUN[3]
      ps[5,i,t,4] <- phiUN[t] * psiUN[4]
      ps[5,i,t,5] <- phiUN[t] * psiUN[5]
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

# modify function for initial z values again to follow annual cycle
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

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = rnorm(4), lpsiJB = rnorm(4), lpsiMI = rnorm(4), lpsiAR = rnorm(4), lpsiUN = rnorm(4),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 1), sigmaJB = runif(1, 0, 1), sigmaMI = runif(1, 0, 1), sigmaAR = runif(1, 0, 1),
                         z = mod.init.z2(rCH, 5))} 


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "pDB", "mean.pDB", "sigma2DB",
                "pJB", "mean.pJB", "sigma2JB",
                "pMI", "mean.pMI", "sigma2MI",
                "pAR", "mean.pAR", "sigma2AR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 7 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: temporal random effect, following fieldwork schedule
# psi: following annual cycle but constant

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
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
marked[,1] <- rep(c(100, 0, 0), n.occasions/3)
marked[,2] <- rep(0, n.occasions)
marked[,3] <- rep(c(0, 20, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 20), n.occasions/3)
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
      for (i in 1:4){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
         psiUN[i] <- exp(lpsiUN[i]) / (1 + exp(lpsiUN[1]) + exp(lpsiUN[2]) + exp(lpsiUN[3]) + exp(lpsiUN[4]))
      }
         
    # Calculate the last transition probability
      psiDB[5] <- 1-psiDB[1]-psiDB[2]-psiDB[3]-psiDB[4]
      psiJB[5] <- 1-psiJB[1]-psiJB[2]-psiJB[3]-psiJB[4]
      psiMI[5] <- 1-psiMI[1]-psiMI[2]-psiMI[3]-psiMI[4]
      psiAR[5] <- 1-psiAR[1]-psiAR[2]-psiAR[3]-psiAR[4]
      psiUN[5] <- 1-psiUN[1]-psiUN[2]-psiUN[3]-psiUN[4]
      
    #   # Normal priors on logit of all but one transition probs
    #   for (i in 1:4){
    #   for (t in 1:(n.occasions-1)){
    #      lpsiUN[i,t] ~ dnorm(0, 0.001)
    #   }
    #   }
    #      
    #   # Constrain the transitions such that their sum is < 1
    #   for (i in 1:4){
    #   for (t in 1:(n.occasions-1)){
    #      psiUN[i,t] <- exp(lpsiUN[i,t]) / (1 + exp(lpsiUN[1,t]) + exp(lpsiUN[2,t]) + exp(lpsiUN[3,t]) + exp(lpsiUN[4,t]))
    #   }
    #   }
    #      
    # # Calculate the last transition probability
    # for (t in 1:(n.occasions-1)){
    #   psiUN[5,t] <- 1-psiUN[1,t]-psiUN[2,t]-psiUN[3,t]-psiUN[4,t]
    # }

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phiDB[t] * psiDB[2]
      ps[1,i,t,3] <- phiDB[t] * psiDB[3]
      ps[1,i,t,4] <- 0
      ps[1,i,t,5] <- phiDB[t] * psiDB[5]
      ps[1,i,t,6] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 0
      ps[2,i,t,4] <- phiJB[t] * psiJB[4]
      ps[2,i,t,5] <- phiJB[t] * psiJB[5]
      ps[2,i,t,6] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phiMI[t] * psiMI[4]
      ps[3,i,t,5] <- phiMI[t] * psiMI[5]
      ps[3,i,t,6] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[1]
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- phiAR[t] * psiAR[5]
      ps[4,i,t,6] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- phiUN[t] * psiUN[1]
      ps[5,i,t,2] <- phiUN[t] * psiUN[2]
      ps[5,i,t,3] <- phiUN[t] * psiUN[3]
      ps[5,i,t,4] <- phiUN[t] * psiUN[4]
      ps[5,i,t,5] <- phiUN[t] * psiUN[5]
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

# modify function for initial z values again to follow annual cycle
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

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = rnorm(4), lpsiJB = rnorm(4), lpsiMI = rnorm(4), lpsiAR = rnorm(4), lpsiUN = rnorm(4),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 1), sigmaJB = runif(1, 0, 1), sigmaMI = runif(1, 0, 1), sigmaAR = runif(1, 0, 1),
                         z = mod.init.z2(rCH, 5))} 


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "pDB","mean.pDB", "sigma2DB", "ptDB",
                "pJB", "mean.pJB", "sigma2JB", "ptJB",
                "pMI", "mean.pMI", "sigma2MI", "ptMI",
                "pAR", "mean.pAR", "sigma2AR", "ptAR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 8 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: temporal random effect, following fieldwork schedule
# psi: following annual cycle and time-dependent

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
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
marked[,1] <- rep(c(100, 0, 0), n.occasions/3)
marked[,2] <- rep(0, n.occasions)
marked[,3] <- rep(c(0, 20, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 20), n.occasions/3)
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
      for (t in 1:(n.occasions-1)){
         lpsiDB[i,t] ~ dnorm(0, 0.001)
         lpsiJB[i,t] ~ dnorm(0, 0.001)
         lpsiMI[i,t] ~ dnorm(0, 0.001)
         lpsiAR[i,t] ~ dnorm(0, 0.001)
         lpsiUN[i,t] ~ dnorm(0, 0.001)
      }
      }

      # Constrain the transitions such that their sum is < 1
      for (i in 1:4){
      for (t in 1:(n.occasions-1)){
         psiDB[i,t] <- exp(lpsiDB[i,t]) / (1 + exp(lpsiDB[1,t]) + exp(lpsiDB[2,t]) + exp(lpsiDB[3,t]) + exp(lpsiDB[4,t]))
         psiJB[i,t] <- exp(lpsiJB[i,t]) / (1 + exp(lpsiJB[1,t]) + exp(lpsiJB[2,t]) + exp(lpsiJB[3,t]) + exp(lpsiJB[4,t]))
         psiMI[i,t] <- exp(lpsiMI[i,t]) / (1 + exp(lpsiMI[1,t]) + exp(lpsiMI[2,t]) + exp(lpsiMI[3,t]) + exp(lpsiMI[4,t]))
         psiAR[i,t] <- exp(lpsiAR[i,t]) / (1 + exp(lpsiAR[1,t]) + exp(lpsiAR[2,t]) + exp(lpsiAR[3,t]) + exp(lpsiAR[4,t]))
         psiUN[i,t] <- exp(lpsiUN[i,t]) / (1 + exp(lpsiUN[1,t]) + exp(lpsiUN[2,t]) + exp(lpsiUN[3,t]) + exp(lpsiUN[4,t]))
      }
      }

    # Calculate the last transition probability
    for (t in 1:(n.occasions-1)){
      psiDB[5,t] <- 1-psiDB[1,t]-psiDB[2,t]-psiDB[3,t]-psiDB[4,t]
      psiJB[5,t] <- 1-psiJB[1,t]-psiJB[2,t]-psiJB[3,t]-psiJB[4,t]
      psiMI[5,t] <- 1-psiMI[1,t]-psiMI[2,t]-psiMI[3,t]-psiMI[4,t]
      psiAR[5,t] <- 1-psiAR[1,t]-psiAR[2,t]-psiAR[3,t]-psiAR[4,t]
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

# modify function for initial z values again to follow annual cycle
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

# Function to create initial values for psi
init.psi <- function(cols){
  psi.array <- array(NA,dim=c(cols,n.occasions-1))
  for (i in 1:cols){
      for (t in 1:(n.occasions-1)){
        psi.array[i,t] <- rnorm(1)
      } #t
  } #i
  return(psi.array)
}

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = init.psi(4), lpsiJB = init.psi(4), lpsiMI = init.psi(4), lpsiAR = init.psi(4), lpsiUN = init.psi(4),
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
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

################################################################################

#### SCENARIO 9 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: temporal random effect, following fieldwork schedule
# psi: following annual cycle and time-dependent with impossible transitions set to 0

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 15
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
marked[,1] <- rep(c(100, 0, 0), n.occasions/3)
marked[,2] <- rep(0, n.occasions)
marked[,3] <- rep(c(0, 20, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 20), n.occasions/3)
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
      for (t in 1:(n.occasions-1)){
         lpsiDB[i,t] ~ dnorm(0, 0.001)
         lpsiJB[i,t] ~ dnorm(0, 0.001)
         lpsiMI[i,t] ~ dnorm(0, 0.001)
         lpsiAR[i,t] ~ dnorm(0, 0.001)
         lpsiUN[i,t] ~ dnorm(0, 0.001)
      }
      }

      # Constrain the transitions such that their sum is < 1
      for (t in 1:(n.occasions-1)){
         psiDB[1,t] <- 0
         psiDB[2,t] <- exp(lpsiDB[2,t]) / (1 + exp(lpsiDB[1,t]) + exp(lpsiDB[2,t]) + exp(lpsiDB[3,t]) + exp(lpsiDB[4,t]))
         psiDB[3,t] <- exp(lpsiDB[3,t]) / (1 + exp(lpsiDB[1,t]) + exp(lpsiDB[2,t]) + exp(lpsiDB[3,t]) + exp(lpsiDB[4,t]))
         psiDB[4,t] <- 0
      }
      
      for (t in 1:(n.occasions-1)){
         psiJB[1,t] <- 0
         psiJB[2,t] <- 0
         psiJB[3,t] <- 0
         psiJB[4,t] <- exp(lpsiJB[4,t]) / (1 + exp(lpsiJB[1,t]) + exp(lpsiJB[2,t]) + exp(lpsiJB[3,t]) + exp(lpsiJB[4,t]))
      }
      
      for (t in 1:(n.occasions-1)){
         psiMI[1,t] <- 0
         psiMI[2,t] <- 0
         psiMI[3,t] <- 0
         psiMI[4,t] <- exp(lpsiMI[4,t]) / (1 + exp(lpsiMI[1,t]) + exp(lpsiMI[2,t]) + exp(lpsiMI[3,t]) + exp(lpsiMI[4,t]))
      }
      
      for (t in 1:(n.occasions-1)){
         psiAR[1,t] <- exp(lpsiAR[1,t]) / (1 + exp(lpsiAR[1,t]) + exp(lpsiAR[2,t]) + exp(lpsiAR[3,t]) + exp(lpsiAR[4,t]))
         psiAR[2,t] <- 0
         psiAR[3,t] <- 0
         psiAR[4,t] <- 0
      }
      
      for (i in 1:4){   
      for (t in 1:(n.occasions-1)){   
         psiUN[i,t] <- exp(lpsiUN[i,t]) / (1 + exp(lpsiUN[1,t]) + exp(lpsiUN[2,t]) + exp(lpsiUN[3,t]) + exp(lpsiUN[4,t]))
      }
      }
      

    # Calculate the last transition probability
    for (t in 1:(n.occasions-1)){
      psiDB[5,t] <- 1-psiDB[1,t]-psiDB[2,t]-psiDB[3,t]-psiDB[4,t]
      psiJB[5,t] <- 1-psiJB[1,t]-psiJB[2,t]-psiJB[3,t]-psiJB[4,t]
      psiMI[5,t] <- 1-psiMI[1,t]-psiMI[2,t]-psiMI[3,t]-psiMI[4,t]
      psiAR[5,t] <- 1-psiAR[1,t]-psiAR[2,t]-psiAR[3,t]-psiAR[4,t]
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

# modify function for initial z values again to follow annual cycle
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

# Function to create initial values for psi
init.psi <- function(cols){
  psi.array <- array(NA,dim=c(cols,n.occasions-1))
  for (i in 1:cols){
    for (t in 1:(n.occasions-1)){
      psi.array[i,t] <- rnorm(1)
    } #t
  } #i
  return(psi.array)
}

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1),
                         lpsiDB = init.psi(4), lpsiJB = init.psi(4), lpsiMI = init.psi(4), lpsiAR = init.psi(4), lpsiUN = init.psi(4),
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
                           n.thin = nt, save.all.iter = TRUE, Rhat.limit = 1.1,
                           max.iter = 500000, parallel = TRUE)

print(rekn.ms.update, digits = 3)

rekn.ms.summary <- rekn.ms.update$summary

saveRDS(rekn.ms.update$summary, file = paste0("./analysis-output/ms-simulation-summary", Sys.Date(), ".rds"))

saveRDS(rekn.ms.update$sims.list, file = paste0("./analysis-output/ms-simulation-simslist", Sys.Date(), ".rds"))

################################################################################