################################################################################
# SINGLE STATE
#
################################################################################
library(rjags)
library(jagsUI)

set.seed(42)

# 7.9. Simulate data
# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break		# If dead, move to next individual 
      # Bernoulli trial: is individual recaptured? 
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Define parameter values
n.occasions <- 12                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
phi <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p <- c(0.4, 0.65, 0.4, 0.45, 0.55, 0.68, 0.66, 0.28, 0.55, 0.45, 0.35)

# Define matrices with survival and recapture probabilities 
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# 7.10. Fitting the CJS to data in the m-array format: the multinomial likelihood
# 7.10.1. Introduction
# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

# 7.10.2. Time-dependent models
# Specify model in BUGS language
sink("cjs-mnl.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   phi[t] ~ dunif(0, 1)         # Priors for survival
   p[t] ~ dunif(0, 1)           # Priors for recapture
   }
# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]                # Probability of non-recapture
   pr[t,t] <- phi[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } #t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])
}
",fill = TRUE)
sink()

# Create the m-array from the capture-histories
marr <- marray(CH)

# Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2], rel = rowSums(marr))

# Initial values
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0, 1), p = runif(dim(marr)[2]-1, 0, 1))}  

# Parameters monitored
parameters <- c("phi", "p", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs <- jags(jags.data, inits, parameters, "cjs-mnl.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(cjs, digits = 3) 

# Evaluation of fit
plot(cjs$BUGSoutput$sims.list$fit, cjs$BUGSoutput$sims.list$fit.new, xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,  
     ylim = c(5, 25), xlim = c(5, 25), bty ="n") 
abline(0, 1, col = "black", lwd = 2)
mean(cjs$BUGSoutput$sims.list$fit.new > cjs$BUGSoutput$sims.list$fit)

################################################################################
# SINGLE STATE TEMPORAL RANDOM EFFECT
#
################################################################################

# 7.4.2. Random time effects
# Define parameter values
n.occasions <- 20                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
var.phi <- 1                       # Temporal variance of survival
p <- rep(0.4, n.occasions-1)

# Determine annual survival probabilities
logit.phi <- rnorm(n.occasions-1, qlogis(mean.phi), var.phi^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-temp-raneff.jags")
cat("
model {

# Priors and constraints
for (t in 1:(n.occasions-1)){
      logit(phi[t]) <- mu + epsilon[t]
      p[t] <- mean.p
      } #t
      
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   }

#mu ~ dnorm(0, 0.001)                    # Prior for logit of mean survival
#mean.phi <- 1 / (1+exp(-mu))            # Logit transformation
mean.phi ~ dunif(0, 1)                   # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi))       # Logit transformation
sigma ~ dunif(0, 10)                     # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)                  # Temporal variance
mean.p ~ dunif(0, 1)                     # Prior for mean recapture

# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   }
# Define the cell probabilities of the m-array
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]                # Probability of non-recapture
   pr[t,t] <- phi[t]*p[t]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j] <- 0
      } #j
   } #t
# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } #t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
   } #t
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])
}
",fill = TRUE)
sink()

# Create the m-array from the capture-histories
marr <- marray(CH)

# Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2], rel = rowSums(marr))

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 10), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "phi", "mean.p", "sigma2", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs <- jags(jags.data, inits, parameters, "cjs-temp-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(cjs, digits = 3) 

sims.list <- cjs$sims.list
sims.list <- as.data.frame(sims.list)

# evaluation of fit
mean(sims.list$fit.new > sims.list$fit)

ppcheck <- ggplot(sims.list, aes(x = fit, y = fit.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

ppcheck

phi.sim <- tibble(time = 1:19,
                  value = phi)

# format resighting
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(-mean.phi, -mean.p, -deviance, -sigma2, -fit, -fit.new) %>% 
  pivot_longer(cols = 1:19, names_to = "time", values_to = "estimate") %>% 
  mutate(time = str_sub(time, 5))

# plot survival
phi.plot <- ggplot() +
  geom_violin(phi.mod, mapping = aes(x = time, y = estimate, group = time, fill = time), alpha = 0.6) +
  geom_point(phi.sim, mapping = aes(x = time, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Time") +
  ylab("Survival probability") +
  theme(legend.position = "none")

par.sim <- tibble(par = c("mean.phi", "mean.p"),
                  value = c(mean.phi, unique(p)))

par.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi, mean.p) %>% 
  pivot_longer(cols = 1:2, names_to = "par", values_to = "estimate")

par.plot <- ggplot() +
  geom_violin(par.mod, mapping = aes(x = par, y = estimate, group = par, fill = par), alpha = 0.6) +
  geom_point(par.sim, mapping = aes(x = par, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Parameter") +
  ylab("Parameter estimate") +
  theme(legend.position = "none")

################################################################################
# MULTISTATE
#
################################################################################

set.seed(42)

# 9.2. Estimation of movement between two sites
# 9.2.1. Model description
# 9.2.2. Generation of simulated data
# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
phiA <- 0.8
phiB <- 0.7
psiAB <- 0.3
psiBA <- 0.5
pA <- 0.7
pB <- 0.4
n.occasions <- 6
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(100, n.occasions)  
marked[,2] <- rep(60, n.occasions)
marked[,3] <- rep(0, n.occasions)

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
      phiA*(1-psiAB), phiA*psiAB,     1-phiA,
      phiB*psiBA,     phiB*(1-psiBA), 1-phiB,
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  1-pA,
      0,  pB, 1-pB,
      0,  0,  1       ), nrow = n.states, byrow = TRUE)
  } #t
} #i

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

# Execute function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

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

marr <- marray(CH)

# run model with state-space formulation
# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3


# 9.2.3. Analysis of the model
# Specify model in BUGS language
sink("ms.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# psiAB: movement probability from site A to site B
# psiBA: movement probability from site B to site A
# pA: recapture probability at site A
# pB: recapture probability at site B
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 dead
# Observations (O):  
# 1 seen at A 
# 2 seen at B
# 3 not seen
# -------------------------------------------------

# Priors and constraints
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   }
for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
   }

# Define state-transition and observation matrices
for (i in 1:nind){  
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phiA[t] * (1-psiAB[t])
      ps[1,i,t,2] <- phiA[t] * psiAB[t]
      ps[1,i,t,3] <- 1-phiA[t]
      ps[2,i,t,1] <- phiB[t] * psiBA[t]
      ps[2,i,t,2] <- phiB[t] * (1-psiBA[t])
      ps[2,i,t,3] <- 1-phiB[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- pA[t]
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-pA[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pB[t]
      po[2,i,t,3] <- 1-pB[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
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

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
  # notseen: label for ‘not seen’
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(state)
}

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1], z = known.state.ms(rCH, 3))

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2, 0, 1), mean.p = runif(2, 0, 1), z = ms.init.z(rCH, f))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 8 min)
ms <- jags(jags.data, inits, parameters, "ms.jags",
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(ms, digits = 3)

# # run model with m-array formulation
# sink("ms-mnl.jags")
# cat("
# model {
# # Priors and constraints
# for (t in 1:(n.occasions-1)){
#    phiA[t] <- mean.phi[1]
#    phiB[t] <- mean.phi[2]
#    psiAB[t] <- mean.psi[1]
#    psiBA[t] <- mean.psi[2]
#    pA[t] <- mean.p[1]
#    pB[t] <- mean.p[2]
# }
#    
# for (u in 1:2){
#    mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
#    mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
#    mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
# }
# 
# # Define state-transition and observation matrices
#   for (t in 1:(n.occasions-1)){
#       ps[1,t,1] <- phiA[t] * (1-psiAB[t])
#       ps[1,t,2] <- phiA[t] * psiAB[t]
#       ps[1,t,3] <- 1-phiA[t]
#       ps[2,t,1] <- phiB[t] * psiBA[t]
#       ps[2,t,2] <- phiB[t] * (1-psiBA[t])
#       ps[2,t,3] <- 1-phiB[t]
#       ps[3,t,1] <- 0
#       ps[3,t,2] <- 0
#       ps[3,t,3] <- 1
#       
#       # Define probabilities of O(t) given S(t)
#       po[1,t,1] <- pA[t]
#       po[1,t,2] <- 0
#       po[1,t,3] <- 1-pA[t]
#       po[2,t,1] <- 0
#       po[2,t,2] <- pB[t]
#       po[2,t,3] <- 1-pB[t]
#       po[3,t,1] <- 0
#       po[3,t,2] <- 0
#       po[3,t,3] <- 1
#       
#       # for (s in 2:n.states){
#       #   for (u in 1:n.states){
#       #     po[s,t,u] <- po[1,t,u]
#       #   } # u
#       # } # s
#       
#       for (s in 1:n.states){
#         for (u in 1:n.states){
#           qo[s,t,u] <- 1-po[s,t,u]
#         } # u
#       } # s
#     } #t
#    
# # Define the multinomial likelihood
# for (t in 1:(n.occasions-1)*n.states){
#    marr[t,1:(n.occasions*n.states-(n.states-1))] ~ dmulti(pr[t, ], rel[t])
#    }
# 
#     # Define the cell probabilities of the m-array
#     # Define matrix Q: product of probabilities of survival and non-capture 
#     for (t in 1:(n.occasions-2)){
#       Q[(t-1)*n.states+(1:n.states), (t-1)*n.states+(1:n.states)] <- ones                         
#       for (j in (t+1):(n.occasions-1)){
#         Q[(t-1)*n.states+(1:n.states), (j-1)*n.states+(1:n.states)] <-
#         Q[(t-1)*n.states+(1:n.states), (j-2)*n.states+(1:n.states)] %*% (ps[,t,] * qo[,t,])     
#       }
#     }
#     Q[(n.occasions-2)*n.states+(1:n.states), (n.occasions-2)*n.states+(1:n.states)] <- ones
#     
#     # Define the cell probabilities of the multistate m-array  
#     # The main diagonal
#     for (t in 1:(n.occasions-2)){
#       pr[(t-1)*n.states+(1:n.states),(t-1)*n.states+(1:n.states)] <-
#       Q[(t-1)*n.states+(1:n.states), (t-1)*n.states+(1:n.states)] %*% (ps[,t,] * po[,t,])
#       # Above main diagonal
#       for (j in (t+1):(n.occasions-1)){
#         pr[(t-1)*n.states+(1:n.states), (j-1)*n.states+(1:n.states)] <-
#         Q[(t-1)*n.states+(1:n.states), (j-1)*n.states+(1:n.states)] %*% (ps[,j,] * po[,j,])
#       }
#     }
#     pr[(n.occasions-2)*n.states+(1:n.states), (n.occasions-2)*n.states+(1:n.states)] <- 
#     ps[,n.occasions-1,] * po[,n.occasions-1,]       
#     
#     # Below main diagonal
#     for (t in 2:(n.occasions-1)){
#       for (j in 1:(t-1)){
#         pr[(t-1)*n.states+(1:n.states),(j-1)*n.states+(1:n.states)] <- zero
#       } #j
#     } #t
#     
#     # Last column: probability of non-recapture
#     for (t in 1:((n.occasions-1)*n.states)){
#       pr[t,(n.occasions*n.states-(n.states-1))] <- 1-sum(pr[t,1:((n.occasions-   
#       1)*n.states)])
#     } #t
# 
# # # Assess model fit using Freeman-Tukey statistic
# # # Compute fit statistics for observed data
# # for (t in 1:(n.occasions-1)){
# #    for (j in 1:n.occasions){
# #       expmarr[t,j] <- rel[t]*pr[t,j]
# #       E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
# #       } #j
# #    } #t
# # # Generate replicate data and compute fit stats from them
# # for (t in 1:(n.occasions-1)){
# #    marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
# #    for (j in 1:n.occasions){
# #       E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
# #       } #j
# #    } #t
# # fit <- sum(E.org[,])
# # fit.new <- sum(E.new[,])
# }
# ",fill = TRUE)
# sink()

# model code 2
sink("ms-mnl.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   psiAB[t] <- mean.psi[1]
   psiBA[t] <- mean.psi[2]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
}
   
for (u in 1:2){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.psi[u] ~ dunif(0, 1)    # Priors for mean transitions
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
}
# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- phiA[t] * (1-psiAB[t])
      psi[1,t,2] <- phiA[t] * psiAB[t]
      psi[2,t,1] <- phiB[t] * psiBA[t]
      psi[2,t,2] <- phiB[t] * (1-psiBA[t])

      po[1,t] <- pA[t]
      po[2,t] <- pB[t]


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
}
",fill = TRUE)
sink()

# Create the m-array from the capture-histories
marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

# Bundle data
jags.data <- list(marr = marr, n.occasions = dim(CH)[2], rel = rowSums(marr),
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns))

# Initial values
inits <- function(){list(mean.phi = runif(2, 0, 1), mean.psi = runif(2, 0, 1), mean.p = runif(2, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.psi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 8 min)
ms <- jags(jags.data, inits, parameters, "ms-mnl.jags",
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(ms, digits = 3)

################################################################################
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
pC <- 0.5
n.occasions <- 6
n.states <- 4
n.obs <- 4
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(50, n.occasions)  
marked[,2] <- rep(50, n.occasions)
marked[,3] <- rep(50, n.occasions)
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
      pA, 0,  0,  1-pA,
      0,  pB, 0,  1-pB,
      0,  0,  pC, 1-pC,
      0,  0,  0,  1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute vector with occasions of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# 9.6.3. Analysis of the model
# Specify model in BUGS language
sink("ms3-multinomlogit.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# phiC: survival probability at site C
# psiAB: movement probability from site A to site B
# psiAC: movement probability from site A to site C
# psiBA: movement probability from site B to site A
# psiBC: movement probability from site B to site C
# psiCA: movement probability from site C to site A
# psiCB: movement probability from site C to site B
# pA: recapture probability at site A
# pB: recapture probability at site B
# pC: recapture probability at site C 
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 alive at C
# 4 dead
# Observations (O):
# 1 seen at A 
# 2 seen at B
# 3 seen at C
# 4 not seen
# -------------------------------------------------

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   phiC[t] <- mean.phi[3]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   pC[t] <- mean.p[3]
   }
   
   for (u in 1:3){
   mean.phi[u] ~ dunif(0, 1)
   mean.p[u] ~ dunif(0, 1)
   }

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probas
      for (t in 1:(n.occasions-1)){
      for (i in 1:2){
         lpsiA[t,i] <- mean.lpsi[1,i]
         lpsiB[t,i] <- mean.lpsi[2,i]
         lpsiC[t,i] <- mean.lpsi[3,i]
        }
      
      # Constrain the transitions such that their sum is < 1
      for (i in 1:2){
         psiA[t,i] <- exp(lpsiA[t,i]) / (1 + exp(lpsiA[t,1]) + exp(lpsiA[t,2]))
         psiB[t,i] <- exp(lpsiB[t,i]) / (1 + exp(lpsiB[t,1]) + exp(lpsiB[t,2]))
         psiC[t,i] <- exp(lpsiC[t,i]) / (1 + exp(lpsiC[t,1]) + exp(lpsiC[t,2]))
         }
      # Calculate the last transition probability
      psiA[t,3] <- 1-psiA[t,1]-psiA[t,2]
      psiB[t,3] <- 1-psiB[t,1]-psiB[t,2]
      psiC[t,3] <- 1-psiC[t,1]-psiC[t,2]
    }

  for (u in 1:3){
  for (i in 1:2){
    mean.lpsi[u,i] ~ dnorm(0, 0.001)
    }
  }

# Define state-transition and observation matrices 	
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,t,1] <- phiA[t] * psiA[t,1]
      ps[1,t,2] <- phiA[t] * psiA[t,2]
      ps[1,t,3] <- phiA[t] * psiA[t,3]
      ps[1,t,4] <- 1-phiA[t]
      ps[2,t,1] <- phiB[t] * psiB[t,1]
      ps[2,t,2] <- phiB[t] * psiB[t,2]
      ps[2,t,3] <- phiB[t] * psiB[t,3]
      ps[2,t,4] <- 1-phiB[t]
      ps[3,t,1] <- phiC[t] * psiC[t,1]
      ps[3,t,2] <- phiC[t] * psiC[t,2]
      ps[3,t,3] <- phiC[t] * psiC[t,3]
      ps[3,t,4] <- 1-phiC[t]

      # Define probabilities of O(t)
      po[1,t] <- pA[t]
      po[2,t] <- pB[t]
      po[3,t] <- pC[t]
      
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
      U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-2)*ns+(1:ns)] %*% ps[,t,] %*% dq[,t,]
      }
   }
U[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
# Diagonal
for (t in 1:(n.occasions-2)){
   pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] %*% ps[,t,] %*% dp[,t,]
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- U[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] %*% ps[,j,] %*% dp[,j,]
      }
   }
pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ps[,n.occasions-1,] %*% dp[,n.occasions-1,]

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
}
",fill = TRUE)
sink()

# model code 2
sink("ms3-multinomlogit.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phiA: survival probability at site A
# phiB: survival probability at site B
# phiC: survival probability at site C
# psiAB: movement probability from site A to site B
# psiAC: movement probability from site A to site C
# psiBA: movement probability from site B to site A
# psiBC: movement probability from site B to site C
# psiCA: movement probability from site C to site A
# psiCB: movement probability from site C to site B
# pA: recapture probability at site A
# pB: recapture probability at site B
# pC: recapture probability at site C 
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 alive at C
# 4 dead
# Observations (O):
# 1 seen at A 
# 2 seen at B
# 3 seen at C
# 4 not seen
# -------------------------------------------------

# Priors and constraints
   # Survival and recapture: uniform
   phiA ~ dunif(0, 1)
   phiB ~ dunif(0, 1)
   phiC ~ dunif(0, 1)
   pA ~ dunif(0, 1)
   pB ~ dunif(0, 1)
   pC ~ dunif(0, 1)

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
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,t,1] <- phiA * psiA[1]
      ps[1,t,2] <- phiA * psiA[2]
      ps[1,t,3] <- phiA * psiA[3]
      ps[1,t,4] <- 1-phiA
      ps[2,t,1] <- phiB * psiB[1]
      ps[2,t,2] <- phiB * psiB[2]
      ps[2,t,3] <- phiB * psiB[3]
      ps[2,t,4] <- 1-phiB
      ps[3,t,1] <- phiC * psiC[1]
      ps[3,t,2] <- phiC * psiC[2]
      ps[3,t,3] <- phiC * psiC[3]
      ps[3,t,4] <- 1-phiC
      ps[4,t,1] <- 0
      ps[4,t,2] <- 0
      ps[4,t,3] <- 0
      ps[4,t,4] <- 1

# Define probabilities of O(t) given S(t)
      po[1,t,1] <- pA
      po[1,t,2] <- 0
      po[1,t,3] <- 0
      po[1,t,4] <- 1-pA
      
      for (s in 2:ns){
        for (u in 1:ns){
          po[s,t,u] <- po[1,t,u]
        } # u
      } # s
      
      for (s in 1:ns){
        for (u in 1:ns){
          qo[s,t,u] <- 1-po[s,t,u]
        } # u
      } # s
    } # t
    
    # Define the multinomial likelihood
    for (t in 1:((n.occasions-1)*ns)){
      marr[t,1:(n.occasions*ns-(ns-1))] ~ dmulti(pr[t, ], rel[t])
    }
    
    # Define the cell probabilities of the m-array
    # Define matrix Q: product of probabilities of survival and non-capture 
    for (t in 1:(n.occasions-2)){
      Q[(t-1)*ns+(1:ns), (t-1)*ns+(1:ns)] <- ones                         
      for (j in (t+1):(n.occasions-1)){
        Q[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- Q[(t-1)*ns+(1:ns), (j-
        2)*ns+(1:ns)] %*% (ps[,t,] * qo[,t,])     
      }
    }
    Q[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- ones
    
    # Define the cell probabilities of the multistate m-array  
    # The main diagonal
    for (t in 1:(n.occasions-2)){
      pr[(t-1)*ns+(1:ns),(t-1)*ns+(1:ns)] <- Q[(t-1)*ns+(1:ns), (t-
      1)*ns+(1:ns)] %*% (ps[,t,] * po[,t,])
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr[(t-1)*ns+(1:ns), (j-1)*ns+(1:ns)] <- Q[(t-1)*ns+(1:ns), (j-
        1)*ns+(1:ns)] %*% (ps[,j,] * po[,j,])
      }
    }
    pr[(n.occasions-2)*ns+(1:ns), (n.occasions-2)*ns+(1:ns)] <- 
    ps[,n.occasions-1,] * po[,n.occasions-1,]       
    
    # Below main diagonal
    for (t in 2:(n.occasions-1)){
      for (j in 1:(t-1)){
        pr[(t-1)*ns+(1:ns),(j-1)*ns+(1:ns)] <- zero
      } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:((n.occasions-1)*ns)){
      pr[t,(n.occasions*ns-(ns-1))] <- 1-sum(pr[t,1:((n.occasions-   
      1)*ns)])
    } #t
        
    }
    ",fill = TRUE)
sink()
      
# Create the m-array from the capture-histories
marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

# Bundle data
jags.data <- list(marr = marr, n.occasions = ncol(CH), rel = rowSums(marr),
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(3, 0, 1),
                         mean.lpsi = matrix(rnorm(1), nrow = 3, ncol = 2),
                         mean.p = runif(3, 0, 1))}  


# Parameters monitored
parameters <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC", "pA", "pB", "pC")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3

# Call JAGS from R (BRT 56 min)
ms3 <- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
