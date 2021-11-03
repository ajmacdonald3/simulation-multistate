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

ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}

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

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.9
phiB <- 0.7
phiC <- 0.7
phiD <- 0.8

psiAB <- 0.4
psiAC <- 0.2
psiAD <- 0.1
psiBA <- 0.1
psiBC <- 0.1
psiBD <- 0.3
psiCA <- 0.5
psiCB <- 0.1
psiCD <- 0.2
psiDA <- 0.5
psiDB <- 0.1
psiDC <- 0.1

pA <- 0.7
pB <- 0.6
pC <- 0.6
pD <- 0.5

n.occasions <- 6
n.states <- 5
n.obs <- 5

marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(50, n.occasions)
marked[,2] <- rep(50, n.occasions)
marked[,3] <- rep(50, n.occasions)
marked[,4] <- rep(50, n.occasions)
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
      phiA*(1-psiAB-psiAC-psiAD), phiA*psiAB,                 phiA*psiAC,                 phiA*psiAD,                 1-phiA,
      phiB*psiBA,                 phiB*(1-psiBA-psiBC-psiBD), phiB*psiBC,                 phiB*psiBD,                 1-phiB,
      phiC*psiCA,                 phiC*psiCB,                 phiC*(1-psiCA-psiCB-psiCD), phiC*psiCD,                 1-phiC,
      phiD*psiDA,                 phiD*psiDB,                 phiD*psiDC,                 phiD*(1-psiDA-psiDB-psiDC), 1-phiD,
      0,                          0,                          0,                          0,                          1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  0,  0,  1-pA,
      0,  pB, 0,  0,  1-pB,
      0,  0,  pC, 0,  1-pC,
      0,  0,  0,  pD, 1-pD,
      0,  0,  0,  0,  1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasions of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3 = seen alive in C, 4 = seen alive in D, 5 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 5

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
   phiD ~ dunif(0, 1)
   pA ~ dunif(0, 1)
   pB ~ dunif(0, 1)
   pC ~ dunif(0, 1)
   pD ~ dunif(0, 1)

   # Transitions: gamma priors
   for (i in 1:4){
      a[i] ~ dgamma(1, 1)
      psiA[i] <- a[i]/sum(a[])
      b[i] ~ dgamma(1, 1)
      psiB[i] <- b[i]/sum(b[])
      c[i] ~ dgamma(1, 1)
      psiC[i] <- c[i]/sum(c[])
      d[i] ~ dgamma(1, 1)
      psiD[i] <- d[i]/sum(d[])
      }

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA * psiA[1]
      ps[1,i,t,2] <- phiA * psiA[2]
      ps[1,i,t,3] <- phiA * psiA[3]
      ps[1,i,t,4] <- phiA * psiA[4]
      ps[1,i,t,5] <- 1-phiA
      ps[2,i,t,1] <- phiB * psiB[1]
      ps[2,i,t,2] <- phiB * psiB[2]
      ps[2,i,t,3] <- phiB * psiB[3]
      ps[2,i,t,4] <- phiB * psiB[4]
      ps[2,i,t,5] <- 1-phiB
      ps[3,i,t,1] <- phiC * psiC[1]
      ps[3,i,t,2] <- phiC * psiC[2]
      ps[3,i,t,3] <- phiC * psiC[3]
      ps[3,i,t,4] <- phiC * psiC[4]
      ps[3,i,t,5] <- 1-phiC
      ps[4,i,t,1] <- phiD * psiD[1]
      ps[4,i,t,2] <- phiD * psiD[2]
      ps[4,i,t,3] <- phiD * psiD[3]
      ps[4,i,t,4] <- phiD * psiD[4]
      ps[4,i,t,5] <- 1-phiD
      ps[5,i,t,1] <- 0
      ps[5,i,t,2] <- 0
      ps[5,i,t,3] <- 0
      ps[5,i,t,4] <- 0
      ps[5,i,t,5] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[1,i,t,4] <- 0
      po[1,i,t,5] <- 1-pA
      po[2,i,t,1] <- pB
      po[2,i,t,2] <- 0
      po[2,i,t,3] <- 0
      po[2,i,t,4] <- 0
      po[2,i,t,5] <- 1-pB
      po[3,i,t,1] <- pC
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 0
      po[3,i,t,4] <- 0
      po[3,i,t,5] <- 1-pC
      po[4,i,t,1] <- pD
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 0
      po[4,i,t,4] <- 0
      po[4,i,t,5] <- 1-pD
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

# Initial values 
inits <- function(){list(phiA = runif(1, 0, 1), phiB = runif(1, 0, 1), phiC = runif(1, 0, 1), phiD = runif(1, 0, 1),
                         #lpsiA = rnorm(3), lpsiB = rnorm(3), lpsiC = rnorm(3), lpsiD = rnorm(3),
                         pA = runif(1, 0, 1) , pB = runif(1, 0, 1), pC = runif(1, 0, 1), pD = runif(1, 0, 1),
                         z = ms.init.z(rCH, f))}  


# Parameters monitored
parameters <- c("phiA", "phiB", "phiC", "phiD",
                "psiA", "psiB", "psiC", "psiD",
                "pA", "pB", "pC", "pD")

# MCMC settings
ni <- 100
nt <- 2
nb <- 20
nc <- 3

# Call JAGS from R (BRT 56 min)
ms3 <- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(ms3, digits = 3)

