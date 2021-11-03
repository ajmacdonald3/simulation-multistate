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
for (t in 1:(n.occasions-1)){
   phiA[t] <- mean.phi[1]
   phiB[t] <- mean.phi[2]
   phiC[t] <- mean.phi[3]
   pA[t] <- mean.p[1]
   pB[t] <- mean.p[2]
   pC[t] <- mean.p[3]
}
   
for (u in 1:3){
   mean.phi[u] ~ dunif(0, 1)    # Priors for mean state-spec. survival
   mean.p[u] ~ dunif(0, 1)      # Priors for mean state-spec. recapture
}

# Transitions: gamma priors
   for (i in 1:3){
      a[i] ~ dgamma(1, 1)
      psiA[i] <- a[i]/sum(a[])
      b[i] ~ dgamma(1, 1)
      psiB[i] <- b[i]/sum(b[])
      c[i] ~ dgamma(1, 1)
      psiC[i] <- c[i]/sum(c[])
      }

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- phiA[t] * psiA[1]
      psi[1,t,2] <- phiA[t] * psiA[2]
      psi[1,t,3] <- phiA[t] * psiA[3]
      #psi[1,t,4] <- 1-phiA[t]
      
      psi[2,t,1] <- phiB[t] * psiB[1]
      psi[2,t,2] <- phiB[t] * psiB[2]
      psi[2,t,3] <- phiB[t] * psiB[3]
      #psi[2,t,4] <- 1-phiB[t] 

      psi[3,t,1] <- phiC[t] * psiC[1]
      psi[3,t,2] <- phiC[t] * psiC[2]
      psi[3,t,3] <- phiC[t] * psiC[3]
      #psi[3,t,4] <- 1-phiC[t]
      
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

# Create the m-array from the capture-histories
marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

# Bundle data
jags.data <- list(marr = marr, n.occasions = ncol(CH), rel = rowSums(marr),
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(3, 0, 1), mean.p = runif(3, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phi", "psiA", "psiB", "psiC", "mean.p", "fit", "fit.new")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3

# Call JAGS from R (BRT 56 min)
ms3 <- jags(jags.data, inits, parameters, "ms3-multinomlogit.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(ms3, digits = 3)

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################
library(tidyverse)
library(cowplot)
library(viridis)

# parameter identifiability checks
sims.list <- ms3$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/plover-example/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
          xlab(i))
  
  dev.off()
  
}

# Evaluation of fit
plot(ms3$sims.list$fit, ms3$sims.list$fit.new, xlab = "Discrepancy actual data",
     ylab = "Discrepancy replicate data", las = 1,  
     bty ="n") 
abline(0, 1, col = "black", lwd = 2)
mean(ms3$sims.list$fit.new > ms3$sims.list$fit)

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("A", "B", "C"),
                  value = c(phiA, phiB, phiC))

# resighting
p.sim <- tibble(site = c("A", "B", "C"),
                value = c(pA, pB, pC))

# transition
psi.sim <- tibble(transition = c("A-B", "A-C", "B-A", "B-C", "C-A", "C-B"),
                  value = c(psiAB, psiAC, psiBA, psiBC, psiCA, psiCB))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3) %>% 
  pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "A")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "B")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "C"))

# plot survival
phi.plot <- ggplot() +
  geom_violin(phi.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_point(phi.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Survival probability") +
  theme(legend.position = "none")

# format resighting
p.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.p.1, mean.p.2, mean.p.3) %>% 
  pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.p.1", "A")) %>%
  mutate(site = str_replace(site, "mean.p.2", "B")) %>%
  mutate(site = str_replace(site, "mean.p.3", "C")) 

# plot resighting
p.plot <- ggplot() +
  geom_violin(p.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_point(p.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

# format transitions
psi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiA.2, psiA.3, psiB.1, psiB.3, psiC.1, psiC.2) %>% 
  pivot_longer(cols = 1:6, names_to = "transition", values_to = "estimate") %>% 
  mutate(transition = str_replace(transition, "psiA.2", "A-B")) %>%
  mutate(transition = str_replace(transition, "psiA.3", "A-C")) %>%
  mutate(transition = str_replace(transition, "psiB.1", "B-A")) %>%
  mutate(transition = str_replace(transition, "psiB.3", "B-C")) %>%
  mutate(transition = str_replace(transition, "psiC.1", "C-A")) %>%
  mutate(transition = str_replace(transition, "psiC.2", "C-B"))

# plot transitions
psi.plot <- ggplot() +
  geom_violin(psi.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_point(psi.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all in grid
sim.plots <- plot_grid(phi.plot, p.plot, labels = c("A", "B"), ncol = 2)

png(filename = "figures/plover-example.png", width = 8, height = 8,
    units = "in", res = 600)

plot_grid(sim.plots, psi.plot, labels = c("", "C"), nrow = 2)

dev.off()

################################################################################
# LADY SLIPPER EXAMPLE
#
################################################################################

# 9.7. Real data example: the showy lady’s slipper
CH <- as.matrix(read.table("orchids.txt", sep=" ", header = F))
n.occasions <- dim(CH)[2]

# Create the m-array from the capture-histories
marr <- marray(CH, unobs = 1)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

# Specify model in BUGS language
sink("ladyslipper.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# s: survival probability
# psiV: transitions from vegetative
# psiF: transitions from flowering
# psiD: transitions from dormant
# -------------------------------------------------
# States (S):
# 1 vegetative
# 2 flowering
# 3 dormant
# 4 dead
# Observations (O):  
# 1 seen vegetative 
# 2 seen flowering
# 3 not seen
# -------------------------------------------------

# Priors and constraints
   # Survival: uniform
   for (t in 1:(n.occasions-1)){  
      s[t] ~ dunif(0, 1)
      }
   # Transitions: gamma priors
   for (i in 1:3){
      a[i] ~ dgamma(1, 1)
      psiD[i] <- a[i]/sum(a[])
      b[i] ~ dgamma(1, 1)
      psiV[i] <- b[i]/sum(b[])
      c[i] ~ dgamma(1, 1)
      psiF[i] <- c[i]/sum(c[])
      }

# Define state-transition and observation matrices 	
   # Define probabilities of state S(t+1) given S(t)
   for (t in 1:(n.occasions-1)){
      ps[1,t,1] <- s[t] * psiV[1]
      ps[1,t,2] <- s[t] * psiV[2]
      ps[1,t,3] <- s[t] * psiV[3]
      #ps[1,i,t,4] <- 1-s[t]
      ps[2,t,1] <- s[t] * psiF[1]
      ps[2,t,2] <- s[t] * psiF[2]
      ps[2,t,3] <- s[t] * psiF[3]
      #ps[2,i,t,4] <- 1-s[t]
      ps[3,t,1] <- s[t] * psiD[1]
      ps[3,t,2] <- s[t] * psiD[2]
      ps[3,t,3] <- s[t] * psiD[3]
      #ps[3,i,t,4] <- 1-s[t]
      #ps[4,i,t,1] <- 0
      #ps[4,i,t,2] <- 0
      #ps[4,i,t,3] <- 0
      #ps[4,i,t,4] <- 1

      # Define probabilities of O(t) given S(t)
      po[1,t] <- 1
      #po[1,i,t,2] <- 0
      #po[1,i,t,3] <- 0
      #po[2,i,t,1] <- 0
      po[2,t] <- 1
      #po[2,i,t,3] <- 0
      #po[3,i,t,1] <- 0
      #po[3,i,t,2] <- 0
      po[3,t] <- 1
      #po[4,i,t,1] <- 0
      #po[4,i,t,2] <- 0
      #po[4,i,t,3] <- 1
      } #t
   
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

# Bundle data
jags.data <- list(marr = marr, n.occasions = ncol(CH), rel = rowSums(marr),
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns))

# Initial values 
inits <- function(){list(s = runif((dim(rCH)[2]-1), 0, 1))}  


# Parameters monitored
parameters <- c("s", "psiV", "psiF", "psiD")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R (BRT 56 min)
ls <- jags(jags.data, inits, parameters, "ladyslipper.jags",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(ls, digits = 3)