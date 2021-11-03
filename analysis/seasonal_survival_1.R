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
# 
# ms.init.z <- function(ch, f){
#   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
#   states <- max(ch, na.rm = TRUE)
#   v <- which(ch==states)
#   ch[-v] <- NA
#   ch[v] <- states
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



# 9.6. Estimation of movement among three sites
# 9.6.1. Model description
# 9.6.2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations
n.occasions <- 15
n.states <- 6
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93
phiUN <- 0.92

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
  
pDB <- rep(c(0, 0, 0.6), n.occasions/3)
pJB <- c(0.5, 0, 0, 0.3, 0, 0, 0.1, 0, 0, 0.5, 0, 0, 0, 0, 0)
pMI <- rep(c(0.4, 0, 0), n.occasions/3)
pAR <- rep(c(0, 0.5, 0), n.occasions/3)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(500, 0, 0), n.occasions/3)
marked[,2] <- rep(0, n.occasions)
marked[,3] <- rep(c(0, 50, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 50), n.occasions/3)
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

# 9.6.3. Analysis of the model
# Specify model in BUGS language
source("analysis/JAGS_seasonal_survival_1c.R")

#### old model code ####
sink("seasonal_survival_1.jags")
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
   
   DB.n <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)
   DB.y <- c(3, 6, 9, 12, 15)
   
   for (t in DB.n){
      pDB[t] <- 0
   }
   
   for (t in DB.y){
      pDB[t] <- mean.pDB
   }
   
   mean.pDB <- dunif(0, 1)
   
   JB.n <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15)
   JB.y <- c(1, 4, 7, 10, 13)
   
   for (t in JB.n){
      pJB[t] <- 0
   }
   
   for (t in JB.y){
      pJB[t] <- dunif(0, 1)
   }
   
   MI.n <- c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15)
   MI.y <- c(1, 4, 7, 10, 13)
   
   for (t in MI.n){
      pMI[t] <- 0
   }
   
   for (t in MI.y){
      pMI[t] <- mean.pMI
   }
   
   mean.pMI <- dunif(0, 1)
   
   AR.n <- c(1, 3, 4, 6, 7, 9, 10, 12, 13, 15)
   AR.y <- c(2, 5, 8, 11, 14)
   
   for (t in AR.n){
      pAR[t] <- 0
   }
   
   for (t in AR.y){
      pAR[t] <- mean.pAR
   }
   
   mean.pAR <- dunif(0, 1)

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probs
      for (t in n.occasions){
      for (i in 1:4){
         lpsiDB[t,i] <- mean.lpsiDB
         lpsiJB[t,i] <- mean.lpsiJB
         lpsiMI[t,i] <- mean.lpsiMI
         lpsiAR[t,i] <- mean.lpsiAR
         lpsiUN[t,i] ~ dnorm(0, 0.001)
      }
         mean.lpsiDB ~ dnorm(0, 0.001)
         mean.lpsiJB ~ dnorm(0, 0.001)
         mean.lpsiMI ~ dnorm(0, 0.001)
         mean.lpsiAR ~ dnorm(0, 0.001)
         
      # Constrain the transitions such that their sum is < 1
      for (i in 1:4){
         psiDB[t,i] <- exp(lpsiDB[t,i]) / (1 + exp(lpsiDB[t,1]) + exp(lpsiDB[t,2]) + exp(lpsiDB[t,3]) + exp(lpsiDB[t,4]))
         psiJB[t,i] <- exp(lpsiJB[t,i]) / (1 + exp(lpsiJB[t,1]) + exp(lpsiJB[t,2]) + exp(lpsiJB[t,3]) + exp(lpsiJB[t,4]))
         psiMI[t,i] <- exp(lpsiMI[t,i]) / (1 + exp(lpsiMI[t,1]) + exp(lpsiMI[t,2]) + exp(lpsiMI[t,3]) + exp(lpsiMI[t,4]))
         psiAR[t,i] <- exp(lpsiAR[t,i]) / (1 + exp(lpsiAR[t,1]) + exp(lpsiAR[t,2]) + exp(lpsiAR[t,3]) + exp(lpsiAR[t,4]))
         psiUN[t,i] <- exp(lpsiUN[t,i]) / (1 + exp(lpsiUN[t,1]) + exp(lpsiUN[t,2]) + exp(lpsiUN[t,3]) + exp(lpsiUN[t,4]))
      }
         
      # Calculate the last transition probability
      psiDB[t,5] <- 1-psiDB[t,1]-psiDB[t,2]-psiDB[t,3]-psiDB[t,4]
      psiJB[t,5] <- 1-psiJB[t,1]-psiJB[t,2]-psiJB[t,3]-psiJB[t,4]
      psiMI[t,5] <- 1-psiMI[t,1]-psiMI[t,2]-psiMI[t,3]-psiMI[t,4]
      psiAR[t,5] <- 1-psiAR[t,1]-psiAR[t,2]-psiAR[t,3]-psiAR[t,4]
      psiUN[t,5] <- 1-psiUN[t,1]-psiUN[t,2]-psiUN[t,3]-psiUN[t,4]
  }

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phiDB * psiDB[t,2]
      ps[1,i,t,3] <- phiDB * psiDB[t,3]
      ps[1,i,t,4] <- 0
      ps[1,i,t,5] <- phiDB * psiDB[t,5]
      ps[1,i,t,6] <- 1-phiDB
      
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 0
      ps[2,i,t,4] <- phiJB * psiJB[t,4]
      ps[2,i,t,5] <- phiJB * psiJB[t,5]
      ps[2,i,t,6] <- 1-phiJB
      
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phiMI * psiMI[t,4]
      ps[3,i,t,5] <- phiMI * psiMI[t,5]
      ps[3,i,t,6] <- 1-phiMI
      
      ps[4,i,t,1] <- phiAR * psiAR[t,1]
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- phiAR * psiAR[t,5]
      ps[4,i,t,5] <- 1-phiAR
      
      ps[5,i,t,1] <- phiUN * psiUN[t,1]
      ps[5,i,t,2] <- phiUN * psiUN[t,2]
      ps[5,i,t,3] <- phiUN * psiUN[t,3]
      ps[5,i,t,4] <- phiUN * psiUN[t,4]
      ps[5,i,t,5] <- phiUN * psiUN[t,5]
      ps[5,i,t,6] <- 1-phiUN     
      
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
#### end old model code ####

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 5))

# Initial values 
inits <- function(){list(mean.phiDB = runif(1, 0, 1), mean.phiJB = runif(1, 0, 1), mean.phiMI = runif(1, 0, 1), mean.phiAR = runif(1, 0, 1), mean.phiUN = runif(1, 0, 1),
                         mean.lpsiDB = rnorm(4*n.occasions), mean.lpsiJB = rnorm(4*n.occasions), mean.lpsiMI = rnorm(4*n.occasions), mean.lpsiAR = rnorm(4*n.occasions), lpsiUN = rnorm(4*n.occasions),
                         mean.pDB = runif(1, 0, 1) , pJB = runif(n.occasions, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         z = mod.init.z(rCH, 5))}  


# Parameters monitored
parameters <- c("mean.phiDB", "mean.phiJB", "mean.phiMI", "mean.phiAR", "mean.phiUN",
                "psiDB", "psiJB", "psiMI", "psiAR", "psiUN",
                "mean.pDB", "pJB", "mean.pMI", "mean.pAR")

# MCMC settings
ni <- 100
nt <- 1
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "seasonal_survival_1.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rekn.ms, digits = 3)

