library(rjags)
library(jagsUI)

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

ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}

# 9.6. Estimation of movement among three sites
# 9.6.1. Model description
# 9.6.2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.85
phiB <- 0.75
phiC <- 0.65
psiAB <- 0.3
psiAC <- 0.2
psiBA <- 0.5
psiBC <- 0.1
psiCA <- 0.6
psiCB <- 0.1
pA <- 0.7
pB <- 0.4
#pC <- 0.5
n.occasions <- 6
n.states <- 4
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(50, n.occasions)  
marked[,2] <- rep(50, n.occasions)
marked[,3] <- rep(0, n.occasions)
marked[,4] <- rep(0, n.occasions)

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
      phiA*(1-psiAB-psiAC), phiA*psiAB,           phiA*psiAC,           1-phiA,
      phiB*psiBA,           phiB*(1-psiBA-psiBC), phiB*psiBC,           1-phiB,
      phiC*psiCA,           phiC*psiCB,           phiC*(1-psiCA-psiCB), 1-phiC,
      0,                    0,                    0,                    1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  1-pA,
      0,  pB, 1-pB,
      0,  0,  1,
      0,  0,  1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = 3)
CH <- sim$CH

# Compute vector with occasions of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 3

# 9.6.3. Analysis of the model
# Specify model in BUGS language
sink("ms3-multinomlogit.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# phiC: survival probability when unobservable
# psiAB: movement probability from site A to site B
# psiAC: movement probability from site A to unobservable
# psiBA: movement probability from site B to site A
# psiBC: movement probability from site B to unobservable
# psiCA: movement probability from unobservable to site A
# psiCB: movement probability from unobservable to site B
# pA: recapture probability at site A
# pB: recapture probability at site B
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 alive but unobservable
# 4 dead
# Observations (O):
# 1 seen at A 
# 2 seen at B
# 3 not seen
# -------------------------------------------------

# Priors and constraints
   # Survival and recapture: uniform
   phiA ~ dunif(0, 1)
   phiB ~ dunif(0, 1)
   phiC ~ dunif(0, 1)
   pA ~ dunif(0, 1)
   pB ~ dunif(0, 1)

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probas
      for (i in 1:2){
         lpsiA[i] ~ dnorm(0, 0.001)
         lpsiB[i] ~ dnorm(0, 0.001)
         lpsiC[i] ~ dnorm(0, 0.001)
         }
      # Constrain the transitions such that their sum is < 1
      for (i in 1:2){
         psiA[i] <- exp(lpsiA[i]) / (1 + exp(lpsiA[1]) + exp(lpsiA[2]))
         psiB[i] <- exp(lpsiB[i]) / (1 + exp(lpsiB[1]) + exp(lpsiB[2]))
         psiC[i] <- exp(lpsiC[i]) / (1 + exp(lpsiC[1]) + exp(lpsiC[2]))
         }
      # Calculate the last transition probability
      psiA[3] <- 1-psiA[1]-psiA[2]
      psiB[3] <- 1-psiB[1]-psiB[2]
      psiC[3] <- 1-psiC[1]-psiC[2]

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA * psiA[1]
      ps[1,i,t,2] <- phiA * psiA[2]
      ps[1,i,t,3] <- phiA * psiA[3]
      ps[1,i,t,4] <- 1-phiA
      ps[2,i,t,1] <- phiB * psiB[1]
      ps[2,i,t,2] <- phiB * psiB[2]
      ps[2,i,t,3] <- phiB * psiB[3]
      ps[2,i,t,4] <- 1-phiB
      ps[3,i,t,1] <- phiC * psiC[1]
      ps[3,i,t,2] <- phiC * psiC[2]
      ps[3,i,t,3] <- phiC * psiC[3]
      ps[3,i,t,4] <- 1-phiC
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-pA
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pB
      po[2,i,t,3] <- 1-pB
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 1
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
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))

# Initial values 
inits <- function(){list(phiA = runif(1, 0, 1), phiB = runif(1, 0, 1), phiC = runif(1, 0, 1), lpsiA = rnorm(2), lpsiB = rnorm(2), lpsiC = rnorm(2), pA = runif(1, 0, 1) , pB = runif(1, 0, 1), z = ms.init.z(rCH, f))}  


# Parameters monitored
parameters <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC", "pA", "pB")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3

# Call JAGS from R (BRT 56 min)
ms3 <- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(ms3, digits = 3)



# Estimation of movement among three sites and unobserved state
# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.85
phiB <- 0.75
phiC <- 0.65
phiU <- 0.55

psiAB <- 0.3
psiAC <- 0.2
psiAU <- 0.1
psiBA <- 0.5
psiBC <- 0.1
psiBU <- 0.1
psiCA <- 0.6
psiCB <- 0.1
psiCU <- 0.1
psiUA <- 0.8
psiUB <- 0.7
psiUC <- 0.6

pA <- 0.6
pB <- 0.4
pC <- 0.5

n.occasions <- 30
n.states <- 5
n.obs <- 4

marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(0, 500, 0), n.occasions/3) # site A  
marked[,2] <- rep(c(0, 0, 50), n.occasions/3) # site B
marked[,3] <- rep(c(100, 0, 0), n.occasions/3) # site C
marked[,4] <- rep(0, n.occasions)  # absent
marked[,5] <- rep(0, n.occasions)  # dead

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
      phiA*(1-psiAB-psiAC-psiAU), phiA*psiAB,                 phiA*psiAC,                 phiA*psiAU,                 1-phiA,
      phiB*psiBA,                 phiB*(1-psiBA-psiBC-psiBU), phiB*psiBC,                 phiB*psiBU,                 1-phiB,
      phiC*psiCA,                 phiC*psiCB,                 phiC*(1-psiCA-psiCB-psiCU), phiC*psiCU,                 1-phiC,
      phiU*psiUA,                 phiU*psiUB,                 phiU*psiUC,                 phiU*(1-psiUA-psiUB-psiUC), 1-phiU,
      0,                          0,                          0,                          0,                          1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  0,  1-pA,
      0,  pB, 0,  1-pB,
      0,  0,  pC, 1-pC,
      0,  0,  0,  1,
      0,  0,  0,  1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = 4)
CH <- sim$CH

# Compute vector with occasions of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3, seen alive in C, 4 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 4

# 2. Accounting for immediate trap response
# 2.2. Generation of simulated data 
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phi <- 0.55
pss <- 0.75
pns <- 0.3
n.occasions <- 10  
n.states <- 3
n.obs <- 2
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)	# Alive, seen
marked[,2] <- rep(0, n.occasions)	# Alive, not seen
marked[,3] <- rep(0, n.occasions)	# Dead 

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
      phi*pss, phi*(1-pss), 1-phi,
      phi*pns, phi*(1-pns), 1-phi,
      0,       0,           1       ), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      1, 0,
      0, 1,
      0, 1  ), nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = 2)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen, 2 = not seen 
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 2