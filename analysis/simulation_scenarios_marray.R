################################################################################
# MULTISTATE SURVIVAL SIMULATION SCENARIOS
#
# using m-array format
#
# gradually increasing model complexity to simulate realistic biological
# scenario
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

#### SCENARIO 1 ####

# States: DB, JB, MI, AR
# phi: all constant
# p: all constant
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 30
n.states <- 5
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93

pDB <- 0.6
pJB <- 0.4
pMI <- 0.5
pAR <- 0.5

psiDB.DB <- 0.2
psiDB.JB <- 0.3
psiDB.MI <- 0.4

psiJB.DB <- 0.1
psiJB.JB <- 0.4
psiJB.MI <- 0.2

psiMI.DB <- 0.5
psiMI.JB <- 0.1
psiMI.MI <- 0.1

psiAR.DB <- 0.4
psiAR.JB <- 0.2
psiAR.MI <- 0.2

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(200, n.occasions)
marked[,2] <- rep(50, n.occasions)
marked[,3] <- rep(80, n.occasions)
marked[,4] <- rep(100, n.occasions)
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
      pDB, 0,      0,      0,      1-pDB,
      0,      pJB, 0,      0,      1-pJB,
      0,      0,      pMI, 0,      1-pMI,
      0,      0,      0,      pAR, 1-pAR,
      0,      0,      0,      0,      1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = NA)
CH <- sim$CH

marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phi[1]
   phiJB[t] <- mean.phi[2]
   phiMI[t] <- mean.phi[3]
   phiAR[t] <- mean.phi[4]
   }

   for (u in 1:4){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   for (t in 1:(n.occasions-1)){
   pDB[t] <- mean.p[1]
   pJB[t] <- mean.p[2]
   pMI[t] <- mean.p[3]
   pAR[t] <- mean.p[4]
   }
   
   for (u in 1:4){
   mean.p[u] ~ dunif(0, 1)
   }

# Transitions: gamma priors
   
   for (i in 1:4){
      db[i] ~ dgamma(1, 1)
      psiDB[i] <- db[i]/sum(db[])
      
      jb[i] ~ dgamma(1, 1)
      psiJB[i] <- jb[i]/sum(jb[])
      
      mi[i] ~ dgamma(1, 1)
      psiMI[i] <- mi[i]/sum(mi[])
      
      ar[i] ~ dgamma(1, 1)
      psiAR[i] <- ar[i]/sum(ar[])
   }
   

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- phiDB[t] * psiDB[1]
      psi[1,t,2] <- phiDB[t] * psiDB[2]
      psi[1,t,3] <- phiDB[t] * psiDB[3]
      psi[1,t,4] <- phiDB[t] * psiDB[4]
      
      psi[2,t,1] <- phiJB[t] * psiJB[1]
      psi[2,t,2] <- phiJB[t] * psiJB[2]
      psi[2,t,3] <- phiJB[t] * psiJB[3]
      psi[2,t,4] <- phiJB[t] * psiJB[4]

      psi[3,t,1] <- phiMI[t] * psiMI[1]
      psi[3,t,2] <- phiMI[t] * psiMI[2]
      psi[3,t,3] <- phiMI[t] * psiMI[3]
      psi[3,t,4] <- phiMI[t] * psiMI[4]
      
      psi[4,t,1] <- phiAR[t] * psiAR[1]
      psi[4,t,2] <- phiAR[t] * psiAR[2]
      psi[4,t,3] <- phiAR[t] * psiAR[3]
      psi[4,t,4] <- phiAR[t] * psiAR[4]
      
      po[1,t] <- pDB[t]
      po[2,t] <- pJB[t]
      po[3,t] <- pMI[t]
      po[4,t] <- pAR[t]

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
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(4, 0, 1), mean.p = runif(4, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiAR", "mean.p", "fit", "fit.new")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 70000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/scenario-1-marray-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/scenario-1-marray-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/marray-scenario-1/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
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

png(filename = "figures/scenario-1.marray.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "AR"),
                  value = c(phiDB, phiJB, phiMI, phiAR))

# resighting
p.sim <- tibble(site = c("DB", "JB", "MI", "AR"),
                value = c(pDB, pJB, pMI, pAR))

# transition
psi.sim <- tibble(site1 = c("DB", "DB", "DB", "DB", "JB", "JB", "JB", "JB", "MI", "MI", "MI", "MI", "AR", "AR", "AR", "AR"),
                  site2 = c("DB", "JB", "MI", "AR", "DB", "JB", "MI", "AR", "DB", "JB", "MI", "AR", "DB", "JB", "MI", "AR"),
                  value = c(psiDB.DB, psiDB.JB, psiDB.MI, 1-psiDB.DB-psiDB.JB-psiDB.MI, psiJB.DB, psiJB.JB, psiJB.MI, 1-psiJB.DB-psiJB.JB-psiJB.MI,
                            psiMI.DB, psiMI.JB, psiMI.MI, 1-psiMI.DB-psiMI.JB-psiMI.MI, psiAR.DB, psiAR.JB, psiAR.MI, 1-psiAR.DB-psiAR.JB-psiAR.MI))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phi.4", "AR"))

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
  select(mean.p.1, mean.p.2, mean.p.3, mean.p.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.p.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.p.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.p.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.p.4", "AR"))

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
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiJB.1, psiJB.2, psiJB.3, psiJB.4,
         psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiAR.1, psiAR.2, psiAR.3, psiAR.4) %>% 
  pivot_longer(cols = 1:16, names_to = "site1", values_to = "estimate") %>%
  mutate(site2 = site1) %>% 
  mutate(site2 = str_replace(site2, "psiDB.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiDB.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiDB.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiDB.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiJB.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiJB.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiJB.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiJB.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiMI.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiMI.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiMI.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiMI.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiAR.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiAR.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiAR.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiAR.4", "AR")) %>% 
  mutate(site1 = str_replace(site1, "psiDB.1", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.2", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.3", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.4", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.1", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.2", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.3", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.4", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiMI.1", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.2", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.3", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.4", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiAR.1", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.2", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.3", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.4", "AR"))

# plot transitions
psi.plot <- ggplot() +
  geom_violin(psi.mod,
              mapping = aes(x = site2, y = estimate, group = site2,
                            fill = site2), alpha = 0.6) +
  geom_point(psi.sim,
             mapping = aes(x = site2, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(. ~ site1, ncol = 2) +
  #coord_flip() +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/scenario-1.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/scenario-1.marray.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(p.plot)

dev.off()

# transition
png(filename = "figures/scenario-1.marray.psi.png", width = 8, height = 8,
    units = "in", res = 600)

print(psi.plot)

dev.off()

#### SCENARIO 2 ####

# States: DB, JB, MI, AR
# phi: all constant
# p: temporal random effect
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 30
n.states <- 5
n.obs <- 5

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiAR <- 0.93

mean.pDB <- 0.6
mean.pJB <- 0.4
mean.pMI <- 0.5
mean.pAR <- 0.5

var.p <- 0.5

logit.pDB <- rnorm(n.occasions-1, qlogis(mean.pDB), var.p^0.5)
pDB <- plogis(logit.pDB)
logit.pJB <- rnorm(n.occasions-1, qlogis(mean.pJB), var.p^0.5)
pJB <- plogis(logit.pJB)
logit.pMI <- rnorm(n.occasions-1, qlogis(mean.pMI), var.p^0.5)
pMI <- plogis(logit.pMI)
logit.pAR <- rnorm(n.occasions-1, qlogis(mean.pAR), var.p^0.5)
pAR <- plogis(logit.pAR)

psiDB.DB <- 0.2
psiDB.JB <- 0.3
psiDB.MI <- 0.4

psiJB.DB <- 0.1
psiJB.JB <- 0.4
psiJB.MI <- 0.2

psiMI.DB <- 0.5
psiMI.JB <- 0.1
psiMI.MI <- 0.1

psiAR.DB <- 0.4
psiAR.JB <- 0.2
psiAR.MI <- 0.2

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(200, n.occasions)
marked[,2] <- rep(50, n.occasions)
marked[,3] <- rep(80, n.occasions)
marked[,4] <- rep(100, n.occasions)
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

marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phi[1]
   phiJB[t] <- mean.phi[2]
   phiMI[t] <- mean.phi[3]
   phiAR[t] <- mean.phi[4]
   }

   for (u in 1:4){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   for (t in 1:(n.occasions-1)){
   logit(pDB[t]) <- muDB + epsilonDB[t]
   epsilonDB[t] ~ dnorm(0, tauDB)T(-15, 15)
   logit(pJB[t]) <- muJB + epsilonJB[t]
   epsilonJB[t] ~ dnorm(0, tauJB)T(-15, 15)
   logit(pMI[t]) <- muMI + epsilonMI[t]
   epsilonMI[t] ~ dnorm(0, tauMI)T(-15, 15)
   logit(pAR[t]) <- muAR + epsilonAR[t]
   epsilonAR[t] ~ dnorm(0, tauAR)T(-15, 15)
   }
   
   muDB <- log(mean.pDB / (1-mean.pDB))     # Logit transformation
   mean.pDB ~ dunif(0, 1)                   # Prior for mean survival
   tauDB <- pow(sigmaDB, -2)
   sigmaDB ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2DB <- pow(sigmaDB, 2)              # Temporal variance
   
   muJB <- log(mean.pJB / (1-mean.pJB))     # Logit transformation
   mean.pJB ~ dunif(0, 1)                   # Prior for mean survival
   tauJB <- pow(sigmaJB, -2)
   sigmaJB ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2JB <- pow(sigmaJB, 2)              # Temporal variance
   
   muMI <- log(mean.pMI / (1-mean.pMI))     # Logit transformation
   mean.pMI ~ dunif(0, 1)                   # Prior for mean survival
   tauMI <- pow(sigmaMI, -2)
   sigmaMI ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2MI <- pow(sigmaMI, 2)              # Temporal variance
   
   muAR <- log(mean.pAR / (1-mean.pAR))     # Logit transformation
   mean.pAR ~ dunif(0, 1)                   # Prior for mean survival
   tauAR <- pow(sigmaAR, -2)
   sigmaAR ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2AR <- pow(sigmaAR, 2)              # Temporal variance

# Transitions: gamma priors
   
   for (i in 1:4){
      db[i] ~ dgamma(1, 1)
      psiDB[i] <- db[i]/sum(db[])
      
      jb[i] ~ dgamma(1, 1)
      psiJB[i] <- jb[i]/sum(jb[])
      
      mi[i] ~ dgamma(1, 1)
      psiMI[i] <- mi[i]/sum(mi[])
      
      ar[i] ~ dgamma(1, 1)
      psiAR[i] <- ar[i]/sum(ar[])
   }
   

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- phiDB[t] * psiDB[1]
      psi[1,t,2] <- phiDB[t] * psiDB[2]
      psi[1,t,3] <- phiDB[t] * psiDB[3]
      psi[1,t,4] <- phiDB[t] * psiDB[4]
      
      psi[2,t,1] <- phiJB[t] * psiJB[1]
      psi[2,t,2] <- phiJB[t] * psiJB[2]
      psi[2,t,3] <- phiJB[t] * psiJB[3]
      psi[2,t,4] <- phiJB[t] * psiJB[4]

      psi[3,t,1] <- phiMI[t] * psiMI[1]
      psi[3,t,2] <- phiMI[t] * psiMI[2]
      psi[3,t,3] <- phiMI[t] * psiMI[3]
      psi[3,t,4] <- phiMI[t] * psiMI[4]
      
      psi[4,t,1] <- phiAR[t] * psiAR[1]
      psi[4,t,2] <- phiAR[t] * psiAR[2]
      psi[4,t,3] <- phiAR[t] * psiAR[3]
      psi[4,t,4] <- phiAR[t] * psiAR[4]
      
      po[1,t] <- pDB[t]
      po[2,t] <- pJB[t]
      po[3,t] <- pMI[t]
      po[4,t] <- pAR[t]

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
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(4, 0, 1),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 5), sigmaJB = runif(1, 0, 5), sigmaMI = runif(1, 0, 5), sigmaAR = runif(1, 0, 5))}  


# Parameters monitored
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiAR",
                "pDB", "pJB", "pMI", "pAR", "mean.pDB", "mean.pJB", "mean.pMI", "mean.pAR",
                "fit", "fit.new")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 70000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/scenario-2-marray-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/scenario-2-marray-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/marray-scenario-2/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
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

png(filename = "figures/scenario-2.marray.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "AR"),
                  value = c(phiDB, phiJB, phiMI, phiAR))

# resighting
p.sim <- tibble(site = c(rep("DB", 29), rep("JB", 29), rep("MI", 29), rep("AR", 29)),
                time = rep(1:29, 4),
                value = c(pDB, pJB, pMI, pAR))

mean.p.sim <- tibble(site = c("DB", "JB", "MI", "AR"),
                     value = c(mean.pDB, mean.pJB, mean.pMI, mean.pAR))

# transition
psi.sim <- tibble(site1 = c("DB", "DB", "DB", "DB", "JB", "JB", "JB", "JB", "MI", "MI", "MI", "MI", "AR", "AR", "AR", "AR"),
                  site2 = c("DB", "JB", "MI", "AR", "DB", "JB", "MI", "AR", "DB", "JB", "MI", "AR", "DB", "JB", "MI", "AR"),
                  value = c(psiDB.DB, psiDB.JB, psiDB.MI, 1-psiDB.DB-psiDB.JB-psiDB.MI, psiJB.DB, psiJB.JB, psiJB.MI, 1-psiJB.DB-psiJB.JB-psiJB.MI,
                            psiMI.DB, psiMI.JB, psiMI.MI, 1-psiMI.DB-psiMI.JB-psiMI.MI, psiAR.DB, psiAR.JB, psiAR.MI, 1-psiAR.DB-psiAR.JB-psiAR.MI))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phi.4", "AR"))

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
  select(-psiDB.1, -psiDB.2, -psiDB.3, -psiDB.4, -psiJB.1, -psiJB.2, -psiJB.3, -psiJB.4,
         -psiMI.1, -psiMI.2, -psiMI.3, -psiMI.4, -psiAR.1, -psiAR.2, -psiAR.3, -psiAR.4,
         -mean.phi.1, -mean.phi.2, -mean.phi.3, -mean.phi.4, -deviance,
         -mean.pDB, -mean.pJB, -mean.pMI, -mean.pAR) %>% 
  pivot_longer(cols = 1:116, names_to = "site", values_to = "estimate") %>% 
  mutate(time = str_sub(site, 5)) %>% 
  mutate(site = str_sub(site, 2, 3))

# plot resighting
p.plot <- ggplot() +
  geom_violin(p.mod, mapping = aes(x = time, y = estimate, group = time, fill = time), alpha = 0.6) +
  geom_point(p.sim, mapping = aes(x = time, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(. ~ site, nrow = 4) +
  #coord_flip() +
  xlab("Site") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

# format mean resighting
mean.p.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pDB, mean.pJB, mean.pMI, mean.pAR) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pDB", "DB")) %>%
  mutate(site = str_replace(site, "mean.pJB", "JB")) %>%
  mutate(site = str_replace(site, "mean.pMI", "MI")) %>% 
  mutate(site = str_replace(site, "mean.pAR", "AR"))

# plot mean resighting
mean.p.plot <- ggplot() +
  geom_violin(mean.p.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
  geom_point(mean.p.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Site") +
  ylab("Mean resighting probability") +
  theme(legend.position = "none")

# format transitions
psi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiJB.1, psiJB.2, psiJB.3, psiJB.4,
         psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiAR.1, psiAR.2, psiAR.3, psiAR.4) %>% 
  pivot_longer(cols = 1:16, names_to = "site1", values_to = "estimate") %>%
  mutate(site2 = site1) %>% 
  mutate(site2 = str_replace(site2, "psiDB.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiDB.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiDB.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiDB.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiJB.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiJB.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiJB.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiJB.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiMI.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiMI.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiMI.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiMI.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiAR.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiAR.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiAR.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiAR.4", "AR")) %>% 
  mutate(site1 = str_replace(site1, "psiDB.1", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.2", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.3", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.4", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.1", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.2", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.3", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.4", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiMI.1", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.2", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.3", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.4", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiAR.1", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.2", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.3", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.4", "AR"))

# plot transitions
psi.plot <- ggplot() +
  geom_violin(psi.mod,
              mapping = aes(x = site2, y = estimate, group = site2,
                            fill = site2), alpha = 0.6) +
  geom_point(psi.sim,
             mapping = aes(x = site2, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(. ~ site1, ncol = 2) +
  #coord_flip() +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/scenario-2.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/scenario-2.marray.p.png", width = 8, height = 4,
    units = "in", res = 600)

print(p.plot)

dev.off()

png(filename = "figures/scenario-2.marray.mean.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.p.plot)

dev.off()

# transition
png(filename = "figures/scenario-2.marray.psi.png", width = 8, height = 8,
    units = "in", res = 600)

print(psi.plot)

dev.off()

#### SCENARIO 3 ####

# States: DB, JB, MI, AR, UN
# phi: all constant
# p: all constant
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 30
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

psiDB.DB <- 0.2
psiDB.JB <- 0.3
psiDB.MI <- 0.3
psiDB.AR <- 0.1

psiJB.DB <- 0.1
psiJB.JB <- 0.4
psiJB.MI <- 0.2
psiJB.AR <- 0.2

psiMI.DB <- 0.5
psiMI.JB <- 0.1
psiMI.MI <- 0.1
psiMI.AR <- 0.2

psiAR.DB <- 0.4
psiAR.JB <- 0.2
psiAR.MI <- 0.2
psiAR.AR <- 0.1

psiUN.DB <- 0.3
psiUN.JB <- 0.2
psiUN.MI <- 0.2
psiUN.AR <- 0.2

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(200, n.occasions)
marked[,2] <- rep(50, n.occasions)
marked[,3] <- rep(80, n.occasions)
marked[,4] <- rep(100, n.occasions)
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
      pDB, 0,      0,      0,      1-pDB,
      0,      pJB, 0,      0,      1-pJB,
      0,      0,      pMI, 0,      1-pMI,
      0,      0,      0,      pAR, 1-pAR,
      0,      0,      0,      0,   1,
      0,      0,      0,      0,   1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = 5)
CH <- sim$CH

marr <- marray(CH, unobs = 1)

# Calculate the number of states
unobs <- 1
ns <- length(unique(as.numeric(CH))) - 1 + unobs

sink("simulation_scenarios.jags")
cat("
model {

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
   
   for (t in 1:(n.occasions-1)){
   pDB[t] <- mean.p[1]
   pJB[t] <- mean.p[2]
   pMI[t] <- mean.p[3]
   pAR[t] <- mean.p[4]
   }
   
   for (u in 1:4){
   mean.p[u] ~ dunif(0, 1)
   }

# Transitions: gamma priors
   
   for (i in 1:5){
      db[i] ~ dgamma(1, 1)
      psiDB[i] <- db[i]/sum(db[])
      
      jb[i] ~ dgamma(1, 1)
      psiJB[i] <- jb[i]/sum(jb[])
      
      mi[i] ~ dgamma(1, 1)
      psiMI[i] <- mi[i]/sum(mi[])
      
      ar[i] ~ dgamma(1, 1)
      psiAR[i] <- ar[i]/sum(ar[])
      
      un[i] ~ dgamma(1, 1)
      psiUN[i] <- un[i]/sum(un[])
   }
   

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- phiDB[t] * psiDB[1]
      psi[1,t,2] <- phiDB[t] * psiDB[2]
      psi[1,t,3] <- phiDB[t] * psiDB[3]
      psi[1,t,4] <- phiDB[t] * psiDB[4]
      psi[1,t,5] <- phiDB[t] * psiDB[5]
      
      psi[2,t,1] <- phiJB[t] * psiJB[1]
      psi[2,t,2] <- phiJB[t] * psiJB[2]
      psi[2,t,3] <- phiJB[t] * psiJB[3]
      psi[2,t,4] <- phiJB[t] * psiJB[4]
      psi[2,t,5] <- phiJB[t] * psiJB[5]

      psi[3,t,1] <- phiMI[t] * psiMI[1]
      psi[3,t,2] <- phiMI[t] * psiMI[2]
      psi[3,t,3] <- phiMI[t] * psiMI[3]
      psi[3,t,4] <- phiMI[t] * psiMI[4]
      psi[3,t,5] <- phiMI[t] * psiMI[5]
      
      psi[4,t,1] <- phiAR[t] * psiAR[1]
      psi[4,t,2] <- phiAR[t] * psiAR[2]
      psi[4,t,3] <- phiAR[t] * psiAR[3]
      psi[4,t,4] <- phiAR[t] * psiAR[4]
      psi[4,t,5] <- phiAR[t] * psiAR[5]
      
      psi[5,t,1] <- phiUN[t] * psiUN[1]
      psi[5,t,2] <- phiUN[t] * psiUN[2]
      psi[5,t,3] <- phiUN[t] * psiUN[3]
      psi[5,t,4] <- phiUN[t] * psiUN[4]
      psi[5,t,5] <- phiUN[t] * psiUN[5]
      
      po[1,t] <- pDB[t]
      po[2,t] <- pJB[t]
      po[3,t] <- pMI[t]
      po[4,t] <- pAR[t]
      po[5,t] <- 0

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
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(5, 0, 1), mean.p = runif(4, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiAR", "psiUN", "mean.p", "fit", "fit.new")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 70000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/scenario-3-marray-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/scenario-3-marray-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/marray-scenario-3/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
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

png(filename = "figures/scenario-3.marray.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "AR", "UN"),
                  value = c(phiDB, phiJB, phiMI, phiAR, phiUN))

# resighting
p.sim <- tibble(site = c("DB", "JB", "MI", "AR"),
                value = c(pDB, pJB, pMI, pAR))

# transition
psi.sim <- tibble(site1 = c("DB", "DB", "DB", "DB", "DB", "JB", "JB", "JB", "JB", "JB", "MI", "MI", "MI", "MI", "MI", "AR", "AR", "AR", "AR", "AR", "UN", "UN", "UN", "UN", "UN"),
                  site2 = c("DB", "JB", "MI", "AR", "UN", "DB", "JB", "MI", "AR", "UN", "DB", "JB", "MI", "AR", "UN", "DB", "JB", "MI", "AR", "UN", "DB", "JB", "MI", "AR", "UN"),
                  value = c(psiDB.DB, psiDB.JB, psiDB.MI, psiDB.AR, 1-psiDB.DB-psiDB.JB-psiDB.MI-psiDB.AR, psiJB.DB, psiJB.JB, psiJB.MI, psiJB.AR, 1-psiJB.DB-psiJB.JB-psiJB.MI-psiJB.AR,
                            psiMI.DB, psiMI.JB, psiMI.MI, psiMI.AR, 1-psiMI.DB-psiMI.JB-psiMI.MI-psiMI.AR, psiAR.DB, psiAR.JB, psiAR.MI, psiAR.AR, 1-psiAR.DB-psiAR.JB-psiAR.MI-psiAR.AR,
                            psiUN.DB, psiUN.JB, psiUN.MI, psiUN.AR, 1-psiUN.DB-psiUN.JB-psiUN.MI-psiUN.AR))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4, mean.phi.5) %>% 
  pivot_longer(cols = 1:5, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phi.4", "AR")) %>% 
  mutate(site = str_replace(site, "mean.phi.5", "UN"))

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
  select(mean.p.1, mean.p.2, mean.p.3, mean.p.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.p.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.p.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.p.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.p.4", "AR"))

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
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiDB.5, psiJB.1, psiJB.2, psiJB.3, psiJB.4, psiJB.5,
         psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiMI.5, psiAR.1, psiAR.2, psiAR.3, psiAR.4, psiAR.5,
         psiUN.1, psiUN.2, psiUN.3, psiUN.4, psiUN.5) %>% 
  pivot_longer(cols = 1:25, names_to = "site1", values_to = "estimate") %>%
  mutate(site2 = site1) %>% 
  mutate(site2 = str_replace(site2, "psiDB.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiDB.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiDB.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiDB.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiDB.5", "UN")) %>%
  mutate(site2 = str_replace(site2, "psiJB.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiJB.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiJB.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiJB.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiJB.5", "UN")) %>%
  mutate(site2 = str_replace(site2, "psiMI.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiMI.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiMI.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiMI.4", "AR")) %>%
  mutate(site2 = str_replace(site2, "psiMI.5", "UN")) %>%
  mutate(site2 = str_replace(site2, "psiAR.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiAR.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiAR.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiAR.4", "AR")) %>% 
  mutate(site2 = str_replace(site2, "psiAR.5", "UN")) %>%
  mutate(site2 = str_replace(site2, "psiUN.1", "DB")) %>%
  mutate(site2 = str_replace(site2, "psiUN.2", "JB")) %>%
  mutate(site2 = str_replace(site2, "psiUN.3", "MI")) %>%
  mutate(site2 = str_replace(site2, "psiUN.4", "AR")) %>% 
  mutate(site2 = str_replace(site2, "psiUN.5", "UN")) %>%
  mutate(site1 = str_replace(site1, "psiDB.1", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.2", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.3", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.4", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiDB.5", "DB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.1", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.2", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.3", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.4", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiJB.5", "JB")) %>%
  mutate(site1 = str_replace(site1, "psiMI.1", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.2", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.3", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.4", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiMI.5", "MI")) %>%
  mutate(site1 = str_replace(site1, "psiAR.1", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.2", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.3", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiAR.4", "AR")) %>% 
  mutate(site1 = str_replace(site1, "psiAR.5", "AR")) %>%
  mutate(site1 = str_replace(site1, "psiUN.1", "UN")) %>%
  mutate(site1 = str_replace(site1, "psiUN.2", "UN")) %>%
  mutate(site1 = str_replace(site1, "psiUN.3", "UN")) %>%
  mutate(site1 = str_replace(site1, "psiUN.4", "UN")) %>% 
  mutate(site1 = str_replace(site1, "psiUN.5", "UN"))

# plot transitions
psi.plot <- ggplot() +
  geom_violin(psi.mod,
              mapping = aes(x = site2, y = estimate, group = site2,
                            fill = site2), alpha = 0.6) +
  geom_point(psi.sim,
             mapping = aes(x = site2, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(. ~ site1, ncol = 2) +
  #coord_flip() +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/scenario-3.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/scenario-3.marray.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(p.plot)

dev.off()

# transition
png(filename = "figures/scenario-3.marray.psi.png", width = 8, height = 8,
    units = "in", res = 600)

print(psi.plot)

dev.off()

################################################################################
#### SCENARIO 4 ####

# States: DB, JB, MI, CC, SE, BR, AR
# phi: all constant
# p: all constant
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 30
n.states <- 8
n.obs <- 8

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiCC <- 0.85
phiSE <- 0.90
phiBR <- 0.87
phiAR <- 0.93

pDB <- 0.6
pJB <- 0.4
pMI <- 0.5
pCC <- 0.6
pSE <- 0.6
pBR <- 0.4
pAR <- 0.5

psiDB.DB <- 0.20
psiDB.JB <- 0.10
psiDB.MI <- 0.10
psiDB.CC <- 0.15
psiDB.SE <- 0.05
psiDB.BR <- 0.25

psiJB.DB <- 0.10
psiJB.JB <- 0.30
psiJB.MI <- 0.05
psiJB.CC <- 0.15
psiJB.SE <- 0.20
psiJB.BR <- 0.10

psiMI.DB <- 0.30
psiMI.JB <- 0.05
psiMI.MI <- 0.20
psiMI.CC <- 0.05
psiMI.SE <- 0.15
psiMI.BR <- 0.20

psiCC.DB <- 0.10
psiCC.JB <- 0.15
psiCC.MI <- 0.20
psiCC.CC <- 0.25
psiCC.SE <- 0.05
psiCC.BR <- 0.05

psiSE.DB <- 0.05
psiSE.JB <- 0.25
psiSE.MI <- 0.15
psiSE.CC <- 0.10
psiSE.SE <- 0.20
psiSE.BR <- 0.10

psiBR.DB <- 0.05
psiBR.JB <- 0.10
psiBR.MI <- 0.10
psiBR.CC <- 0.15
psiBR.SE <- 0.20
psiBR.BR <- 0.10

psiAR.DB <- 0.20
psiAR.JB <- 0.10
psiAR.MI <- 0.15
psiAR.CC <- 0.05
psiAR.SE <- 0.10
psiAR.BR <- 0.10

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(80, n.occasions)
marked[,2] <- rep(20, n.occasions)
marked[,3] <- rep(30, n.occasions)
marked[,4] <- rep(20, n.occasions)
marked[,5] <- rep(30, n.occasions)
marked[,6] <- rep(20, n.occasions)
marked[,7] <- rep(50, n.occasions)
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
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.CC, phiDB*psiDB.SE, phiDB*psiDB.BR, phiDB*(1-psiDB.DB-psiDB.JB-psiDB.MI-psiDB.CC-psiDB.SE-psiDB.BR), 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.CC, phiJB*psiJB.SE, phiJB*psiJB.BR, phiJB*(1-psiJB.DB-psiJB.JB-psiJB.MI-psiJB.CC-psiJB.SE-psiJB.BR), 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.CC, phiMI*psiMI.SE, phiMI*psiMI.BR, phiMI*(1-psiMI.DB-psiMI.JB-psiMI.MI-psiMI.CC-psiMI.SE-psiMI.BR), 1-phiMI,
      phiCC*psiCC.DB, phiCC*psiCC.JB, phiCC*psiCC.MI, phiCC*psiCC.CC, phiCC*psiCC.SE, phiCC*psiCC.BR, phiCC*(1-psiCC.DB-psiCC.JB-psiCC.MI-psiCC.CC-psiCC.SE-psiCC.BR), 1-phiCC,
      phiSE*psiSE.DB, phiSE*psiSE.JB, phiSE*psiSE.MI, phiSE*psiSE.CC, phiSE*psiSE.SE, phiSE*psiSE.BR, phiSE*(1-psiSE.DB-psiSE.JB-psiSE.MI-psiSE.CC-psiSE.SE-psiSE.BR), 1-phiSE,
      phiBR*psiBR.DB, phiBR*psiBR.JB, phiBR*psiBR.MI, phiBR*psiBR.CC, phiBR*psiBR.SE, phiBR*psiBR.BR, phiBR*(1-psiBR.DB-psiBR.JB-psiBR.MI-psiBR.CC-psiBR.SE-psiBR.BR), 1-phiBR,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.CC, phiAR*psiAR.SE, phiAR*psiAR.BR, phiAR*(1-psiAR.DB-psiAR.JB-psiAR.MI-psiAR.CC-psiAR.SE-psiAR.BR), 1-phiAR,
      0,              0,              0,              0,              0,              0,              0,                                                               1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pDB, 0,   0,   0,   0,   0,   0,   1-pDB,
      0,   pJB, 0,   0,   0,   0,   0,   1-pJB,
      0,   0,   pMI, 0,   0,   0,   0,   1-pMI,
      0,   0,   0,   pCC, 0,   0,   0,   1-pCC,
      0,   0,   0,   0,   pSE, 0,   0,   1-pSE,
      0,   0,   0,   0,   0,   pBR, 0,   1-pBR,
      0,   0,   0,   0,   0,   0,   pAR, 1-pAR,
      0,   0,   0,   0,   0,   0,   0,   1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = NA)
CH <- sim$CH

marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phi[1]
   phiJB[t] <- mean.phi[2]
   phiMI[t] <- mean.phi[3]
   phiCC[t] <- mean.phi[4]
   phiSE[t] <- mean.phi[5]
   phiBR[t] <- mean.phi[6]
   phiAR[t] <- mean.phi[7]
   }

   for (u in 1:7){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   for (t in 1:(n.occasions-1)){
   pDB[t] <- mean.p[1]
   pJB[t] <- mean.p[2]
   pMI[t] <- mean.p[3]
   pCC[t] <- mean.p[4]
   pSE[t] <- mean.p[5]
   pBR[t] <- mean.p[6]
   pAR[t] <- mean.p[7]
   }
   
   for (u in 1:7){
   mean.p[u] ~ dunif(0, 1)
   }

# Transitions: gamma priors
   
   for (i in 1:7){
      db[i] ~ dgamma(1, 1)
      psiDB[i] <- db[i]/sum(db[])
      
      jb[i] ~ dgamma(1, 1)
      psiJB[i] <- jb[i]/sum(jb[])
      
      mi[i] ~ dgamma(1, 1)
      psiMI[i] <- mi[i]/sum(mi[])
      
      cc[i] ~ dgamma(1, 1)
      psiCC[i] <- cc[i]/sum(cc[])
      
      se[i] ~ dgamma(1, 1)
      psiSE[i] <- se[i]/sum(se[])
      
      br[i] ~ dgamma(1, 1)
      psiBR[i] <- br[i]/sum(br[])
      
      ar[i] ~ dgamma(1, 1)
      psiAR[i] <- ar[i]/sum(ar[])
   }

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
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(7, 0, 1), mean.p = runif(7, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
                "mean.p", "fit", "fit.new")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 70000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/scenario-4-marray-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/scenario-4-marray-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/marray-scenario-4/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
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

png(filename = "figures/scenario-4.marray.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                  value = c(phiDB, phiJB, phiMI, phiCC, phiSE, phiBR, phiAR))

# resighting
p.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                value = c(pDB, pJB, pMI, pCC, pSE, pBR, pAR))

# transition
psiDB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI, psiDB.CC, psiDB.SE, psiDB.BR,
                              1-psiDB.DB-psiDB.JB-psiDB.MI-psiDB.CC-psiDB.SE-psiDB.BR))
  
psiJB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI, psiJB.CC, psiJB.SE, psiJB.BR,
                              1-psiJB.DB-psiJB.JB-psiJB.MI-psiJB.CC-psiJB.SE-psiJB.BR))

psiMI.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI, psiMI.CC, psiMI.SE, psiMI.BR,
                              1-psiMI.DB-psiMI.JB-psiMI.MI-psiMI.CC-psiMI.SE-psiMI.BR))

psiCC.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiCC.DB, psiCC.JB, psiCC.MI, psiCC.CC, psiCC.SE, psiCC.BR,
                              1-psiCC.DB-psiCC.JB-psiCC.MI-psiCC.CC-psiCC.SE-psiCC.BR))

psiSE.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiSE.DB, psiSE.JB, psiSE.MI, psiSE.CC, psiSE.SE, psiSE.BR,
                              1-psiSE.DB-psiSE.JB-psiSE.MI-psiSE.CC-psiSE.SE-psiSE.BR))

psiBR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiBR.DB, psiBR.JB, psiBR.MI, psiBR.CC, psiBR.SE, psiBR.BR,
                              1-psiBR.DB-psiBR.JB-psiBR.MI-psiBR.CC-psiBR.SE-psiBR.BR))

psiAR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiAR.DB, psiAR.JB, psiAR.MI, psiAR.CC, psiAR.SE, psiAR.BR,
                              1-psiAR.DB-psiAR.JB-psiAR.MI-psiAR.CC-psiAR.SE-psiAR.BR))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4, mean.phi.5, mean.phi.6, mean.phi.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phi.4", "CC")) %>% 
  mutate(site = str_replace(site, "mean.phi.5", "SE")) %>%
  mutate(site = str_replace(site, "mean.phi.6", "BR")) %>%
  mutate(site = str_replace(site, "mean.phi.7", "AR"))

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
  select(mean.p.1, mean.p.2, mean.p.3, mean.p.4, mean.p.5, mean.p.6, mean.p.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.p.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.p.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.p.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.p.4", "CC")) %>% 
  mutate(site = str_replace(site, "mean.p.5", "SE")) %>%
  mutate(site = str_replace(site, "mean.p.6", "BR")) %>%
  mutate(site = str_replace(site, "mean.p.7", "AR"))

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
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiDB.5, psiDB.6, psiDB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiDB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiDB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiDB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiDB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiDB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiDB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiDB.7", "AR"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiJB.1, psiJB.2, psiJB.3, psiJB.4, psiJB.5, psiJB.6, psiJB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiJB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiJB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiJB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiJB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiJB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiJB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiJB.7", "AR"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiMI.5, psiMI.6, psiMI.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiMI.1", "DB")) %>%
  mutate(site = str_replace(site, "psiMI.2", "JB")) %>%
  mutate(site = str_replace(site, "psiMI.3", "MI")) %>%
  mutate(site = str_replace(site, "psiMI.4", "CC")) %>%
  mutate(site = str_replace(site, "psiMI.5", "SE")) %>%
  mutate(site = str_replace(site, "psiMI.6", "BR")) %>%
  mutate(site = str_replace(site, "psiMI.7", "AR"))

psiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiCC.1, psiCC.2, psiCC.3, psiCC.4, psiCC.5, psiCC.6, psiCC.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiCC.1", "DB")) %>%
  mutate(site = str_replace(site, "psiCC.2", "JB")) %>%
  mutate(site = str_replace(site, "psiCC.3", "MI")) %>%
  mutate(site = str_replace(site, "psiCC.4", "CC")) %>%
  mutate(site = str_replace(site, "psiCC.5", "SE")) %>%
  mutate(site = str_replace(site, "psiCC.6", "BR")) %>%
  mutate(site = str_replace(site, "psiCC.7", "AR"))

psiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiSE.1, psiSE.2, psiSE.3, psiSE.4, psiSE.5, psiSE.6, psiSE.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiSE.1", "DB")) %>%
  mutate(site = str_replace(site, "psiSE.2", "JB")) %>%
  mutate(site = str_replace(site, "psiSE.3", "MI")) %>%
  mutate(site = str_replace(site, "psiSE.4", "CC")) %>%
  mutate(site = str_replace(site, "psiSE.5", "SE")) %>%
  mutate(site = str_replace(site, "psiSE.6", "BR")) %>%
  mutate(site = str_replace(site, "psiSE.7", "AR"))

psiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiBR.1, psiBR.2, psiBR.3, psiBR.4, psiBR.5, psiBR.6, psiBR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiBR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiBR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiBR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiBR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiBR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiBR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiBR.7", "AR"))

psiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiAR.1, psiAR.2, psiAR.3, psiAR.4, psiAR.5, psiAR.6, psiAR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiAR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiAR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiAR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiAR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiAR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiAR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiAR.7", "AR"))

# plot transitions
psiDB.plot <- ggplot() +
  geom_violin(psiDB.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
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
  geom_point(psiAR.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/scenario-4.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/scenario-4.marray.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(p.plot)

dev.off()

# transition
png(filename = "figures/scenario-4.marray.psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/scenario-4.marray.psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/scenario-4.marray.psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

png(filename = "figures/scenario-4.marray.psiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiCC.plot)

dev.off()

png(filename = "figures/scenario-4.marray.psiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiSE.plot)

dev.off()

png(filename = "figures/scenario-4.marray.psiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiBR.plot)

dev.off()

png(filename = "figures/scenario-4.marray.psiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiAR.plot)

dev.off()

################################################################################
#### SCENARIO 5 ####

# States: DB, JB, MI, CC, SE, BR, AR
# phi: all constant
# p: all constant
# psi: biologically realistic ones, but possible at all times

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 27
n.states <- 8
n.obs <- 8

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiCC <- 0.85
phiSE <- 0.90
phiBR <- 0.87
phiAR <- 0.93

pDB <- 0.6
pJB <- 0.4
pMI <- 0.5
pCC <- 0.6
pSE <- 0.6
pBR <- 0.4
pAR <- 0.5

psiDB.DB <- 0
psiDB.JB <- 0.35
psiDB.MI <- 0.30
psiDB.CC <- 0.15
psiDB.SE <- 0.10
psiDB.BR <- 0
psiDB.AR <- 0

psiJB.DB <- 0
psiJB.JB <- 0
psiJB.MI <- 0
psiJB.CC <- 0.35
psiJB.SE <- 0.45
psiJB.BR <- 0.10
psiJB.AR <- 0.10

psiMI.DB <- 0
psiMI.JB <- 0
psiMI.MI <- 0
psiMI.CC <- 0.45
psiMI.SE <- 0.35
psiMI.BR <- 0.10
psiMI.AR <- 0.10

psiCC.DB <- 0
psiCC.JB <- 0
psiCC.MI <- 0
psiCC.CC <- 0
psiCC.SE <- 0.30
psiCC.BR <- 0.30
psiCC.AR <- 0.40

psiSE.DB <- 0
psiSE.JB <- 0
psiSE.MI <- 0
psiSE.CC <- 0
psiSE.SE <- 0.60
psiSE.BR <- 0.20
psiSE.AR <- 0.20

psiBR.DB <- 0.70
psiBR.JB <- 0
psiBR.MI <- 0
psiBR.CC <- 0
psiBR.SE <- 0
psiBR.BR <- 0
psiBR.AR <- 0.30

psiAR.DB <- 0.90
psiAR.JB <- 0
psiAR.MI <- 0
psiAR.CC <- 0
psiAR.SE <- 0
psiAR.BR <- 0.10
psiAR.AR <- 0

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(80, n.occasions)
marked[,2] <- rep(20, n.occasions)
marked[,3] <- rep(30, n.occasions)
marked[,4] <- rep(20, n.occasions)
marked[,5] <- rep(30, n.occasions)
marked[,6] <- rep(20, n.occasions)
marked[,7] <- rep(50, n.occasions)
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
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.CC, phiDB*psiDB.SE, phiDB*psiDB.BR, phiDB*psiDB.AR, 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.CC, phiJB*psiJB.SE, phiJB*psiJB.BR, phiJB*psiJB.AR, 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.CC, phiMI*psiMI.SE, phiMI*psiMI.BR, phiMI*psiMI.AR, 1-phiMI,
      phiCC*psiCC.DB, phiCC*psiCC.JB, phiCC*psiCC.MI, phiCC*psiCC.CC, phiCC*psiCC.SE, phiCC*psiCC.BR, phiCC*psiCC.AR, 1-phiCC,
      phiSE*psiSE.DB, phiSE*psiSE.JB, phiSE*psiSE.MI, phiSE*psiSE.CC, phiSE*psiSE.SE, phiSE*psiSE.BR, phiSE*psiSE.AR, 1-phiSE,
      phiBR*psiBR.DB, phiBR*psiBR.JB, phiBR*psiBR.MI, phiBR*psiBR.CC, phiBR*psiBR.SE, phiBR*psiBR.BR, phiBR*psiBR.AR, 1-phiBR,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.CC, phiAR*psiAR.SE, phiAR*psiAR.BR, phiAR*psiAR.AR, 1-phiAR,
      0,              0,              0,              0,              0,              0,              0,              1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pDB, 0,   0,   0,   0,   0,   0,   1-pDB,
      0,   pJB, 0,   0,   0,   0,   0,   1-pJB,
      0,   0,   pMI, 0,   0,   0,   0,   1-pMI,
      0,   0,   0,   pCC, 0,   0,   0,   1-pCC,
      0,   0,   0,   0,   pSE, 0,   0,   1-pSE,
      0,   0,   0,   0,   0,   pBR, 0,   1-pBR,
      0,   0,   0,   0,   0,   0,   pAR, 1-pAR,
      0,   0,   0,   0,   0,   0,   0,   1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked, unobservable = NA)
CH <- sim$CH

marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phi[1]
   phiJB[t] <- mean.phi[2]
   phiMI[t] <- mean.phi[3]
   phiCC[t] <- mean.phi[4]
   phiSE[t] <- mean.phi[5]
   phiBR[t] <- mean.phi[6]
   phiAR[t] <- mean.phi[7]
   }

   for (u in 1:7){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   for (t in 1:(n.occasions-1)){
   pDB[t] <- mean.p[1]
   pJB[t] <- mean.p[2]
   pMI[t] <- mean.p[3]
   pCC[t] <- mean.p[4]
   pSE[t] <- mean.p[5]
   pBR[t] <- mean.p[6]
   pAR[t] <- mean.p[7]
   }
   
   for (u in 1:7){
   mean.p[u] ~ dunif(0, 1)
   }

# Transitions: gamma priors
   for (i in 1:4){
      db[i] ~ dgamma(1, 1)
   }
      psiDB[1] <- 0
      psiDB[2] <- db[1]/sum(db[])
      psiDB[3] <- db[2]/sum(db[])
      psiDB[4] <- db[3]/sum(db[])
      psiDB[5] <- db[4]/sum(db[])
      psiDB[6] <- 0
      psiDB[7] <- 0
      
    for (i in 1:4){
      jb[i] ~ dgamma(1, 1)
    }
      psiJB[1] <- 0
      psiJB[2] <- 0
      psiJB[3] <- 0
      psiJB[4] <- jb[1]/sum(jb[])
      psiJB[5] <- jb[2]/sum(jb[])
      psiJB[6] <- jb[3]/sum(jb[])
      psiJB[7] <- jb[4]/sum(jb[])
      
    for (i in 1:4){
      mi[i] ~ dgamma(1, 1)
    }
      psiMI[1] <- 0
      psiMI[2] <- 0
      psiMI[3] <- 0
      psiMI[4] <- mi[1]/sum(mi[])
      psiMI[5] <- mi[2]/sum(mi[])
      psiMI[6] <- mi[3]/sum(mi[])
      psiMI[7] <- mi[4]/sum(mi[])
      
    for (i in 1:3){
      cc[i] ~ dgamma(1, 1)
    }
      psiCC[1] <- 0
      psiCC[2] <- 0
      psiCC[3] <- 0
      psiCC[4] <- 0
      psiCC[5] <- cc[1]/sum(cc[])
      psiCC[6] <- cc[2]/sum(cc[])
      psiCC[7] <- cc[3]/sum(cc[])
      
    for (i in 1:3){
      se[i] ~ dgamma(1, 1)
    }
      psiSE[1] <- 0
      psiSE[2] <- 0
      psiSE[3] <- 0
      psiSE[4] <- 0
      psiSE[5] <- se[1]/sum(se[])
      psiSE[6] <- se[2]/sum(se[])
      psiSE[7] <- se[3]/sum(se[])
      
    for (i in 1:2){
      br[i] ~ dgamma(1, 1)
    }
      psiBR[1] <- br[1]/sum(br[])
      psiBR[2] <- 0
      psiBR[3] <- 0
      psiBR[4] <- 0
      psiBR[5] <- 0
      psiBR[6] <- 0
      psiBR[7] <- br[2]/sum(br[])
      
    for (i in 1:2){
      ar[i] ~ dgamma(1, 1)
    }
      psiAR[1] <- ar[1]/sum(ar[])
      psiAR[2] <- 0
      psiAR[3] <- 0
      psiAR[4] <- 0
      psiAR[5] <- 0
      psiAR[6] <- ar[2]/sum(ar[])
      psiAR[7] <- 0

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- 0
      psi[1,t,2] <- phiDB[t] * psiDB[2]
      psi[1,t,3] <- phiDB[t] * psiDB[3]
      psi[1,t,4] <- phiDB[t] * psiDB[4]
      psi[1,t,5] <- phiDB[t] * psiDB[5]
      psi[1,t,6] <- 0
      psi[1,t,7] <- 0
      
      psi[2,t,1] <- 0
      psi[2,t,2] <- 0
      psi[2,t,3] <- 0
      psi[2,t,4] <- phiJB[t] * psiJB[4]
      psi[2,t,5] <- phiJB[t] * psiJB[5]
      psi[2,t,6] <- phiJB[t] * psiJB[6]
      psi[2,t,7] <- phiJB[t] * psiJB[7]

      psi[3,t,1] <- 0
      psi[3,t,2] <- 0
      psi[3,t,3] <- 0
      psi[3,t,4] <- phiMI[t] * psiMI[4]
      psi[3,t,5] <- phiMI[t] * psiMI[5]
      psi[3,t,6] <- phiMI[t] * psiMI[6]
      psi[3,t,7] <- phiMI[t] * psiMI[7]
      
      psi[4,t,1] <- 0
      psi[4,t,2] <- 0
      psi[4,t,3] <- 0
      psi[4,t,4] <- 0
      psi[4,t,5] <- phiCC[t] * psiCC[5]
      psi[4,t,6] <- phiCC[t] * psiCC[6]
      psi[4,t,7] <- phiCC[t] * psiCC[7]
      
      psi[5,t,1] <- 0
      psi[5,t,2] <- 0
      psi[5,t,3] <- 0
      psi[5,t,4] <- 0
      psi[5,t,5] <- phiSE[t] * psiSE[5]
      psi[5,t,6] <- phiSE[t] * psiSE[6]
      psi[5,t,7] <- phiSE[t] * psiSE[7]
      
      psi[6,t,1] <- phiBR[t] * psiBR[1]
      psi[6,t,2] <- 0
      psi[6,t,3] <- 0
      psi[6,t,4] <- 0
      psi[6,t,5] <- 0
      psi[6,t,6] <- 0
      psi[6,t,7] <- phiBR[t] * psiBR[7]
      
      psi[7,t,1] <- phiAR[t] * psiAR[1]
      psi[7,t,2] <- 0
      psi[7,t,3] <- 0
      psi[7,t,4] <- 0
      psi[7,t,5] <- 0
      psi[7,t,6] <- phiAR[t] * psiAR[6]
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
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(7, 0, 1), mean.p = runif(7, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
                "mean.p", "fit", "fit.new")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 70000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/scenario-5-marray-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/scenario-5-marray-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/marray-scenario-5/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
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

png(filename = "figures/scenario-5.marray.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                  value = c(phiDB, phiJB, phiMI, phiCC, phiSE, phiBR, phiAR))

# resighting
p.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                value = c(pDB, pJB, pMI, pCC, pSE, pBR, pAR))

# transition
psiDB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI, psiDB.CC, psiDB.SE, psiDB.BR, psiDB.AR))

psiJB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI, psiJB.CC, psiJB.SE, psiJB.BR, psiJB.AR))

psiMI.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI, psiMI.CC, psiMI.SE, psiMI.BR, psiMI.AR))

psiCC.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiCC.DB, psiCC.JB, psiCC.MI, psiCC.CC, psiCC.SE, psiCC.BR, psiCC.AR))

psiSE.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiSE.DB, psiSE.JB, psiSE.MI, psiSE.CC, psiSE.SE, psiSE.BR, psiSE.AR))

psiBR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiBR.DB, psiBR.JB, psiBR.MI, psiBR.CC, psiBR.SE, psiBR.BR,psiBR.AR))

psiAR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiAR.DB, psiAR.JB, psiAR.MI, psiAR.CC, psiAR.SE, psiAR.BR, psiAR.AR))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4, mean.phi.5, mean.phi.6, mean.phi.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phi.4", "CC")) %>% 
  mutate(site = str_replace(site, "mean.phi.5", "SE")) %>%
  mutate(site = str_replace(site, "mean.phi.6", "BR")) %>%
  mutate(site = str_replace(site, "mean.phi.7", "AR"))

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
  select(mean.p.1, mean.p.2, mean.p.3, mean.p.4, mean.p.5, mean.p.6, mean.p.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.p.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.p.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.p.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.p.4", "CC")) %>% 
  mutate(site = str_replace(site, "mean.p.5", "SE")) %>%
  mutate(site = str_replace(site, "mean.p.6", "BR")) %>%
  mutate(site = str_replace(site, "mean.p.7", "AR"))

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
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiDB.5, psiDB.6, psiDB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiDB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiDB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiDB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiDB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiDB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiDB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiDB.7", "AR"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiJB.1, psiJB.2, psiJB.3, psiJB.4, psiJB.5, psiJB.6, psiJB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiJB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiJB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiJB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiJB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiJB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiJB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiJB.7", "AR"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiMI.5, psiMI.6, psiMI.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiMI.1", "DB")) %>%
  mutate(site = str_replace(site, "psiMI.2", "JB")) %>%
  mutate(site = str_replace(site, "psiMI.3", "MI")) %>%
  mutate(site = str_replace(site, "psiMI.4", "CC")) %>%
  mutate(site = str_replace(site, "psiMI.5", "SE")) %>%
  mutate(site = str_replace(site, "psiMI.6", "BR")) %>%
  mutate(site = str_replace(site, "psiMI.7", "AR"))

psiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiCC.1, psiCC.2, psiCC.3, psiCC.4, psiCC.5, psiCC.6, psiCC.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiCC.1", "DB")) %>%
  mutate(site = str_replace(site, "psiCC.2", "JB")) %>%
  mutate(site = str_replace(site, "psiCC.3", "MI")) %>%
  mutate(site = str_replace(site, "psiCC.4", "CC")) %>%
  mutate(site = str_replace(site, "psiCC.5", "SE")) %>%
  mutate(site = str_replace(site, "psiCC.6", "BR")) %>%
  mutate(site = str_replace(site, "psiCC.7", "AR"))

psiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiSE.1, psiSE.2, psiSE.3, psiSE.4, psiSE.5, psiSE.6, psiSE.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiSE.1", "DB")) %>%
  mutate(site = str_replace(site, "psiSE.2", "JB")) %>%
  mutate(site = str_replace(site, "psiSE.3", "MI")) %>%
  mutate(site = str_replace(site, "psiSE.4", "CC")) %>%
  mutate(site = str_replace(site, "psiSE.5", "SE")) %>%
  mutate(site = str_replace(site, "psiSE.6", "BR")) %>%
  mutate(site = str_replace(site, "psiSE.7", "AR"))

psiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiBR.1, psiBR.2, psiBR.3, psiBR.4, psiBR.5, psiBR.6, psiBR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiBR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiBR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiBR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiBR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiBR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiBR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiBR.7", "AR"))

psiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiAR.1, psiAR.2, psiAR.3, psiAR.4, psiAR.5, psiAR.6, psiAR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiAR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiAR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiAR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiAR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiAR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiAR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiAR.7", "AR"))

# plot transitions
psiDB.plot <- ggplot() +
  geom_violin(psiDB.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
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
  geom_point(psiAR.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/scenario-5.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/scenario-5.marray.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(p.plot)

dev.off()

# transition
png(filename = "figures/scenario-5.marray.psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/scenario-5.marray.psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/scenario-5.marray.psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

png(filename = "figures/scenario-5.marray.psiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiCC.plot)

dev.off()

png(filename = "figures/scenario-5.marray.psiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiSE.plot)

dev.off()

png(filename = "figures/scenario-5.marray.psiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiBR.plot)

dev.off()

png(filename = "figures/scenario-5.marray.psiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiAR.plot)

dev.off()

################################################################################
#### SCENARIO 6 ####

# States: DB, JB, MI, CC, SE, BR, AR
# phi: all constant
# p: temporal random effect
# psi: biologically realistic ones, but possible at all times (each site can be transitioned to from at least two sites)

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 27
n.states <- 8
n.obs <- 8

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiCC <- 0.85
phiSE <- 0.90
phiBR <- 0.87
phiAR <- 0.93

mean.pDB <- 0.6
mean.pJB <- 0.4
mean.pMI <- 0.5
mean.pCC <- 0.6
mean.pSE <- 0.6
mean.pBR <- 0.4
mean.pAR <- 0.5

var.p <- 0.5

logit.pDB <- rnorm(n.occasions-1, qlogis(mean.pDB), var.p^0.5)
pDB <- plogis(logit.pDB)
logit.pJB <- rnorm(n.occasions-1, qlogis(mean.pJB), var.p^0.5)
pJB <- plogis(logit.pJB)
logit.pMI <- rnorm(n.occasions-1, qlogis(mean.pMI), var.p^0.5)
pMI <- plogis(logit.pMI)
logit.pCC <- rnorm(n.occasions-1, qlogis(mean.pCC), var.p^0.5)
pCC <- plogis(logit.pCC)
logit.pSE <- rnorm(n.occasions-1, qlogis(mean.pSE), var.p^0.5)
pSE <- plogis(logit.pSE)
logit.pBR <- rnorm(n.occasions-1, qlogis(mean.pBR), var.p^0.5)
pBR <- plogis(logit.pBR)
logit.pAR <- rnorm(n.occasions-1, qlogis(mean.pAR), var.p^0.5)
pAR <- plogis(logit.pAR)

psiDB.DB <- 0
psiDB.JB <- 0.35
psiDB.MI <- 0.30
psiDB.CC <- 0.15
psiDB.SE <- 0.10
psiDB.BR <- 0
psiDB.AR <- 0

psiJB.DB <- 0
psiJB.JB <- 0
psiJB.MI <- 0.10
psiJB.CC <- 0.35
psiJB.SE <- 0.35
psiJB.BR <- 0.10
psiJB.AR <- 0.10

psiMI.DB <- 0
psiMI.JB <- 0.10
psiMI.MI <- 0
psiMI.CC <- 0.35
psiMI.SE <- 0.35
psiMI.BR <- 0.10
psiMI.AR <- 0.10

psiCC.DB <- 0
psiCC.JB <- 0
psiCC.MI <- 0
psiCC.CC <- 0
psiCC.SE <- 0.30
psiCC.BR <- 0.30
psiCC.AR <- 0.40

psiSE.DB <- 0
psiSE.JB <- 0
psiSE.MI <- 0
psiSE.CC <- 0
psiSE.SE <- 0.60
psiSE.BR <- 0.20
psiSE.AR <- 0.20

psiBR.DB <- 0.70
psiBR.JB <- 0
psiBR.MI <- 0
psiBR.CC <- 0
psiBR.SE <- 0
psiBR.BR <- 0
psiBR.AR <- 0.30

psiAR.DB <- 0.90
psiAR.JB <- 0
psiAR.MI <- 0
psiAR.CC <- 0
psiAR.SE <- 0
psiAR.BR <- 0.10
psiAR.AR <- 0

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(80, n.occasions)
marked[,2] <- rep(20, n.occasions)
marked[,3] <- rep(30, n.occasions)
marked[,4] <- rep(20, n.occasions)
marked[,5] <- rep(30, n.occasions)
marked[,6] <- rep(20, n.occasions)
marked[,7] <- rep(50, n.occasions)
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
      phiDB*psiDB.DB, phiDB*psiDB.JB, phiDB*psiDB.MI, phiDB*psiDB.CC, phiDB*psiDB.SE, phiDB*psiDB.BR, phiDB*psiDB.AR, 1-phiDB,
      phiJB*psiJB.DB, phiJB*psiJB.JB, phiJB*psiJB.MI, phiJB*psiJB.CC, phiJB*psiJB.SE, phiJB*psiJB.BR, phiJB*psiJB.AR, 1-phiJB,
      phiMI*psiMI.DB, phiMI*psiMI.JB, phiMI*psiMI.MI, phiMI*psiMI.CC, phiMI*psiMI.SE, phiMI*psiMI.BR, phiMI*psiMI.AR, 1-phiMI,
      phiCC*psiCC.DB, phiCC*psiCC.JB, phiCC*psiCC.MI, phiCC*psiCC.CC, phiCC*psiCC.SE, phiCC*psiCC.BR, phiCC*psiCC.AR, 1-phiCC,
      phiSE*psiSE.DB, phiSE*psiSE.JB, phiSE*psiSE.MI, phiSE*psiSE.CC, phiSE*psiSE.SE, phiSE*psiSE.BR, phiSE*psiSE.AR, 1-phiSE,
      phiBR*psiBR.DB, phiBR*psiBR.JB, phiBR*psiBR.MI, phiBR*psiBR.CC, phiBR*psiBR.SE, phiBR*psiBR.BR, phiBR*psiBR.AR, 1-phiBR,
      phiAR*psiAR.DB, phiAR*psiAR.JB, phiAR*psiAR.MI, phiAR*psiAR.CC, phiAR*psiAR.SE, phiAR*psiAR.BR, phiAR*psiAR.AR, 1-phiAR,
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

marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phi[1]
   phiJB[t] <- mean.phi[2]
   phiMI[t] <- mean.phi[3]
   phiCC[t] <- mean.phi[4]
   phiSE[t] <- mean.phi[5]
   phiBR[t] <- mean.phi[6]
   phiAR[t] <- mean.phi[7]
   }

   for (u in 1:7){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   for (t in 1:(n.occasions-1)){
   logit(pDB[t]) <- muDB + epsilonDB[t]
   epsilonDB[t] ~ dnorm(0, tauDB)T(-15, 15)
   logit(pJB[t]) <- muJB + epsilonJB[t]
   epsilonJB[t] ~ dnorm(0, tauJB)T(-15, 15)
   logit(pMI[t]) <- muMI + epsilonMI[t]
   epsilonMI[t] ~ dnorm(0, tauMI)T(-15, 15)
   logit(pCC[t]) <- muCC + epsilonCC[t]
   epsilonCC[t] ~ dnorm(0, tauCC)T(-15, 15)
   logit(pSE[t]) <- muSE + epsilonSE[t]
   epsilonSE[t] ~ dnorm(0, tauSE)T(-15, 15)
   logit(pBR[t]) <- muBR + epsilonBR[t]
   epsilonBR[t] ~ dnorm(0, tauBR)T(-15, 15)
   logit(pAR[t]) <- muAR + epsilonAR[t]
   epsilonAR[t] ~ dnorm(0, tauAR)T(-15, 15)
   }
   
   muDB <- log(mean.pDB / (1-mean.pDB))     # Logit transformation
   mean.pDB ~ dunif(0, 1)                   # Prior for mean survival
   tauDB <- pow(sigmaDB, -2)
   sigmaDB ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2DB <- pow(sigmaDB, 2)              # Temporal variance
   sigma2DB.real <- sigma2DB * pow(mean.pDB, 2) * pow((1-mean.pDB), 2)
   
   muJB <- log(mean.pJB / (1-mean.pJB))     # Logit transformation
   mean.pJB ~ dunif(0, 1)                   # Prior for mean survival
   tauJB <- pow(sigmaJB, -2)
   sigmaJB ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2JB <- pow(sigmaJB, 2)              # Temporal variance
   sigma2JB.real <- sigma2JB * pow(mean.pJB, 2) * pow((1-mean.pJB), 2)
   
   muMI <- log(mean.pMI / (1-mean.pMI))     # Logit transformation
   mean.pMI ~ dunif(0, 1)                   # Prior for mean survival
   tauMI <- pow(sigmaMI, -2)
   sigmaMI ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2MI <- pow(sigmaMI, 2)              # Temporal variance
   sigma2MI.real <- sigma2MI * pow(mean.pMI, 2) * pow((1-mean.pMI), 2)
   
   muCC <- log(mean.pCC / (1-mean.pCC))     # Logit transformation
   mean.pCC ~ dunif(0, 1)                   # Prior for mean survival
   tauCC <- pow(sigmaCC, -2)
   sigmaCC ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2CC <- pow(sigmaCC, 2)              # Temporal variance
   sigma2CC.real <- sigma2CC * pow(mean.pCC, 2) * pow((1-mean.pCC), 2)
   
   muSE <- log(mean.pSE / (1-mean.pSE))     # Logit transformation
   mean.pSE ~ dunif(0, 1)                   # Prior for mean survival
   tauSE <- pow(sigmaSE, -2)
   sigmaSE ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2SE <- pow(sigmaSE, 2)              # Temporal variance
   sigma2SE.real <- sigma2SE * pow(mean.pSE, 2) * pow((1-mean.pSE), 2)
   
   muBR <- log(mean.pBR / (1-mean.pBR))     # Logit transformation
   mean.pBR ~ dunif(0, 1)                   # Prior for mean survival
   tauBR <- pow(sigmaBR, -2)
   sigmaBR ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2BR <- pow(sigmaBR, 2)              # Temporal variance
   sigma2BR.real <- sigma2BR * pow(mean.pBR, 2) * pow((1-mean.pBR), 2)
   
   muAR <- log(mean.pAR / (1-mean.pAR))     # Logit transformation
   mean.pAR ~ dunif(0, 1)                   # Prior for mean survival
   tauAR <- pow(sigmaAR, -2)
   sigmaAR ~ dunif(0, 5)                    # Prior on standard deviation
   sigma2AR <- pow(sigmaAR, 2)              # Temporal variance
   sigma2AR.real <- sigma2AR * pow(mean.pAR, 2) * pow((1-mean.pAR), 2)

# Transitions: gamma priors
   for (i in 1:4){
      db[i] ~ dgamma(1, 1)
   }
      psiDB[1] <- 0
      psiDB[2] <- db[1]/sum(db[])
      psiDB[3] <- db[2]/sum(db[])
      psiDB[4] <- db[3]/sum(db[])
      psiDB[5] <- db[4]/sum(db[])
      psiDB[6] <- 0
      psiDB[7] <- 0
      
    for (i in 1:5){
      jb[i] ~ dgamma(1, 1)
    }
      psiJB[1] <- 0
      psiJB[2] <- 0
      psiJB[3] <- jb[1]/sum(jb[])
      psiJB[4] <- jb[2]/sum(jb[])
      psiJB[5] <- jb[3]/sum(jb[])
      psiJB[6] <- jb[4]/sum(jb[])
      psiJB[7] <- jb[5]/sum(jb[])
      
    for (i in 1:5){
      mi[i] ~ dgamma(1, 1)
    }
      psiMI[1] <- 0
      psiMI[2] <- mi[1]/sum(mi[])
      psiMI[3] <- 0
      psiMI[4] <- mi[2]/sum(mi[])
      psiMI[5] <- mi[3]/sum(mi[])
      psiMI[6] <- mi[4]/sum(mi[])
      psiMI[7] <- mi[5]/sum(mi[])
      
    for (i in 1:3){
      cc[i] ~ dgamma(1, 1)
    }
      psiCC[1] <- 0
      psiCC[2] <- 0
      psiCC[3] <- 0
      psiCC[4] <- 0
      psiCC[5] <- cc[1]/sum(cc[])
      psiCC[6] <- cc[2]/sum(cc[])
      psiCC[7] <- cc[3]/sum(cc[])
      
    for (i in 1:3){
      se[i] ~ dgamma(1, 1)
    }
      psiSE[1] <- 0
      psiSE[2] <- 0
      psiSE[3] <- 0
      psiSE[4] <- 0
      psiSE[5] <- se[1]/sum(se[])
      psiSE[6] <- se[2]/sum(se[])
      psiSE[7] <- se[3]/sum(se[])
      
    for (i in 1:2){
      br[i] ~ dgamma(1, 1)
    }
      psiBR[1] <- br[1]/sum(br[])
      psiBR[2] <- 0
      psiBR[3] <- 0
      psiBR[4] <- 0
      psiBR[5] <- 0
      psiBR[6] <- 0
      psiBR[7] <- br[2]/sum(br[])
      
    for (i in 1:2){
      ar[i] ~ dgamma(1, 1)
    }
      psiAR[1] <- ar[1]/sum(ar[])
      psiAR[2] <- 0
      psiAR[3] <- 0
      psiAR[4] <- 0
      psiAR[5] <- 0
      psiAR[6] <- ar[2]/sum(ar[])
      psiAR[7] <- 0

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- 0
      psi[1,t,2] <- phiDB[t] * psiDB[2]
      psi[1,t,3] <- phiDB[t] * psiDB[3]
      psi[1,t,4] <- phiDB[t] * psiDB[4]
      psi[1,t,5] <- phiDB[t] * psiDB[5]
      psi[1,t,6] <- 0
      psi[1,t,7] <- 0
      
      psi[2,t,1] <- 0
      psi[2,t,2] <- 0
      psi[2,t,3] <- phiJB[t] * psiJB[3]
      psi[2,t,4] <- phiJB[t] * psiJB[4]
      psi[2,t,5] <- phiJB[t] * psiJB[5]
      psi[2,t,6] <- phiJB[t] * psiJB[6]
      psi[2,t,7] <- phiJB[t] * psiJB[7]

      psi[3,t,1] <- 0
      psi[3,t,2] <- phiMI[t] * psiMI[2]
      psi[3,t,3] <- 0
      psi[3,t,4] <- phiMI[t] * psiMI[4]
      psi[3,t,5] <- phiMI[t] * psiMI[5]
      psi[3,t,6] <- phiMI[t] * psiMI[6]
      psi[3,t,7] <- phiMI[t] * psiMI[7]
      
      psi[4,t,1] <- 0
      psi[4,t,2] <- 0
      psi[4,t,3] <- 0
      psi[4,t,4] <- 0
      psi[4,t,5] <- phiCC[t] * psiCC[5]
      psi[4,t,6] <- phiCC[t] * psiCC[6]
      psi[4,t,7] <- phiCC[t] * psiCC[7]
      
      psi[5,t,1] <- 0
      psi[5,t,2] <- 0
      psi[5,t,3] <- 0
      psi[5,t,4] <- 0
      psi[5,t,5] <- phiSE[t] * psiSE[5]
      psi[5,t,6] <- phiSE[t] * psiSE[6]
      psi[5,t,7] <- phiSE[t] * psiSE[7]
      
      psi[6,t,1] <- phiBR[t] * psiBR[1]
      psi[6,t,2] <- 0
      psi[6,t,3] <- 0
      psi[6,t,4] <- 0
      psi[6,t,5] <- 0
      psi[6,t,6] <- 0
      psi[6,t,7] <- phiBR[t] * psiBR[7]
      
      psi[7,t,1] <- phiAR[t] * psiAR[1]
      psi[7,t,2] <- 0
      psi[7,t,3] <- 0
      psi[7,t,4] <- 0
      psi[7,t,5] <- 0
      psi[7,t,6] <- phiAR[t] * psiAR[6]
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
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns))

# Initial values 
inits <- function(){list(mean.phi = runif(7, 0, 1),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pCC = runif(1, 0, 1),
                         mean.pSE = runif(1, 0, 1), mean.pBR = runif(1, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigmaDB = runif(1, 0, 5), sigmaJB = runif(1, 0, 5), sigmaMI = runif(1, 0, 5), sigmaCC = runif(1, 0, 5),
                         sigmaSE = runif(1, 0, 5), sigmaBR = runif(1, 0, 5), sigmaAR = runif(1, 0, 5))}  


# Parameters monitored
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
                "mean.pDB", "mean.pJB", "mean.pMI", "mean.pCC", "mean.pSE", "mean.pBR", "mean.pAR",
                "sigma2DB.real", "sigma2JB.real", "sigma2MI.real", "sigma2CC.real", "sigma2SE.real", "sigma2BR.real", "sigma2AR.real",
                "fit", "fit.new")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 70000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/scenario-6-marray-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/scenario-6-marray-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/marray-scenario-6/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
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

png(filename = "figures/scenario-6.marray.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                  value = c(phiDB, phiJB, phiMI, phiCC, phiSE, phiBR, phiAR))

# resighting
p.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                value = c(mean.pDB, mean.pJB, mean.pMI, mean.pCC, mean.pSE, mean.pBR, mean.pAR))

# transition
psiDB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI, psiDB.CC, psiDB.SE, psiDB.BR, psiDB.AR))

psiJB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI, psiJB.CC, psiJB.SE, psiJB.BR, psiJB.AR))

psiMI.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI, psiMI.CC, psiMI.SE, psiMI.BR, psiMI.AR))

psiCC.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiCC.DB, psiCC.JB, psiCC.MI, psiCC.CC, psiCC.SE, psiCC.BR, psiCC.AR))

psiSE.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiSE.DB, psiSE.JB, psiSE.MI, psiSE.CC, psiSE.SE, psiSE.BR, psiSE.AR))

psiBR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiBR.DB, psiBR.JB, psiBR.MI, psiBR.CC, psiBR.SE, psiBR.BR,psiBR.AR))

psiAR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiAR.DB, psiAR.JB, psiAR.MI, psiAR.CC, psiAR.SE, psiAR.BR, psiAR.AR))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4, mean.phi.5, mean.phi.6, mean.phi.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phi.4", "CC")) %>% 
  mutate(site = str_replace(site, "mean.phi.5", "SE")) %>%
  mutate(site = str_replace(site, "mean.phi.6", "BR")) %>%
  mutate(site = str_replace(site, "mean.phi.7", "AR"))

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
  select(mean.pDB, mean.pJB, mean.pMI, mean.pCC, mean.pSE, mean.pBR, mean.pAR) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pDB", "DB")) %>%
  mutate(site = str_replace(site, "mean.pJB", "JB")) %>%
  mutate(site = str_replace(site, "mean.pMI", "MI")) %>% 
  mutate(site = str_replace(site, "mean.pCC", "CC")) %>% 
  mutate(site = str_replace(site, "mean.pSE", "SE")) %>%
  mutate(site = str_replace(site, "mean.pBR", "BR")) %>%
  mutate(site = str_replace(site, "mean.pAR", "AR"))

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
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiDB.5, psiDB.6, psiDB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiDB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiDB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiDB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiDB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiDB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiDB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiDB.7", "AR"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiJB.1, psiJB.2, psiJB.3, psiJB.4, psiJB.5, psiJB.6, psiJB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiJB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiJB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiJB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiJB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiJB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiJB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiJB.7", "AR"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiMI.5, psiMI.6, psiMI.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiMI.1", "DB")) %>%
  mutate(site = str_replace(site, "psiMI.2", "JB")) %>%
  mutate(site = str_replace(site, "psiMI.3", "MI")) %>%
  mutate(site = str_replace(site, "psiMI.4", "CC")) %>%
  mutate(site = str_replace(site, "psiMI.5", "SE")) %>%
  mutate(site = str_replace(site, "psiMI.6", "BR")) %>%
  mutate(site = str_replace(site, "psiMI.7", "AR"))

psiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiCC.1, psiCC.2, psiCC.3, psiCC.4, psiCC.5, psiCC.6, psiCC.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiCC.1", "DB")) %>%
  mutate(site = str_replace(site, "psiCC.2", "JB")) %>%
  mutate(site = str_replace(site, "psiCC.3", "MI")) %>%
  mutate(site = str_replace(site, "psiCC.4", "CC")) %>%
  mutate(site = str_replace(site, "psiCC.5", "SE")) %>%
  mutate(site = str_replace(site, "psiCC.6", "BR")) %>%
  mutate(site = str_replace(site, "psiCC.7", "AR"))

psiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiSE.1, psiSE.2, psiSE.3, psiSE.4, psiSE.5, psiSE.6, psiSE.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiSE.1", "DB")) %>%
  mutate(site = str_replace(site, "psiSE.2", "JB")) %>%
  mutate(site = str_replace(site, "psiSE.3", "MI")) %>%
  mutate(site = str_replace(site, "psiSE.4", "CC")) %>%
  mutate(site = str_replace(site, "psiSE.5", "SE")) %>%
  mutate(site = str_replace(site, "psiSE.6", "BR")) %>%
  mutate(site = str_replace(site, "psiSE.7", "AR"))

psiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiBR.1, psiBR.2, psiBR.3, psiBR.4, psiBR.5, psiBR.6, psiBR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiBR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiBR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiBR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiBR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiBR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiBR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiBR.7", "AR"))

psiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiAR.1, psiAR.2, psiAR.3, psiAR.4, psiAR.5, psiAR.6, psiAR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiAR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiAR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiAR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiAR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiAR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiAR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiAR.7", "AR"))

# plot transitions
psiDB.plot <- ggplot() +
  geom_violin(psiDB.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
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
  geom_point(psiAR.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/scenario-6.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/scenario-6.marray.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(p.plot)

dev.off()

# transition
png(filename = "figures/scenario-6.marray.psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/scenario-6.marray.psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/scenario-6.marray.psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

png(filename = "figures/scenario-6.marray.psiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiCC.plot)

dev.off()

png(filename = "figures/scenario-6.marray.psiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiSE.plot)

dev.off()

png(filename = "figures/scenario-6.marray.psiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiBR.plot)

dev.off()

png(filename = "figures/scenario-6.marray.psiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiAR.plot)

dev.off()

################################################################################
#### SCENARIO 7 ####

# States: DB, JB, MI, CC, SE, BR, AR
# phi: all constant
# p: constant, but periods set to zero when no observation
# psi: biologically realistic ones with seasonality

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 38
n.states <- 8
n.obs <- 8

phiDB <- 0.97
phiJB <- 0.89
phiMI <- 0.89
phiCC <- 0.85
phiSE <- 0.90
phiBR <- 0.87
phiAR <- 0.93

pDB <- c(rep(c(0, 0, 0, 0.7), 9), 0, 0)
pJB <- c(rep(c(0.6, 0, 0, 0), 9), 0.6, 0)
pMI <- c(rep(c(0.5, 0, 0, 0), 9), 0.5, 0)
pCC <- c(rep(c(0.6, 0.6, 0, 0), 9), 0.6, 0.6)
pSE <- rep(0.6, n.occasions-1)
pBR <- c(rep(c(0, 0.4, 0.4, 0), 9), 0, 0.4)
pAR <- c(rep(c(0, 0, 0.5, 0), 9), 0, 0)

psiDB.DB <- rep(0, n.occasions)
psiDB.JB <- c(rep(c(0.35, 0, 0, 0), 9), 0.35, 0)
psiDB.MI <- c(rep(c(0.30, 0, 0, 0), 9), 0.30, 0)
psiDB.CC <- c(rep(c(0.15, 0, 0, 0), 9), 0.15, 0)
psiDB.SE <- c(rep(c(0.10, 0, 0, 0), 9), 0.10, 0)
psiDB.BR <- rep(0, n.occasions)
psiDB.AR <- rep(0, n.occasions)

psiJB.DB <- rep(0, n.occasions)
psiJB.JB <- rep(0, n.occasions)
psiJB.MI <- rep(0, n.occasions)
psiJB.CC <- c(rep(c(0, 0.45, 0, 0), 9), 0, 0.45)
psiJB.SE <- c(rep(c(0, 0.35, 0, 0), 9), 0, 0.35)
psiJB.BR <- c(rep(c(0, 0.20, 0, 0), 9), 0, 0.20)
psiJB.AR <- rep(0, n.occasions)

psiMI.DB <- rep(0, n.occasions)
psiMI.JB <- rep(0, n.occasions)
psiMI.MI <- rep(0, n.occasions)
psiMI.CC <- c(rep(c(0, 0.35, 0, 0), 9), 0, 0.35)
psiMI.SE <- c(rep(c(0, 0.45, 0, 0), 9), 0, 0.45)
psiMI.BR <- c(rep(c(0, 0.20, 0, 0), 9), 0, 0.20)
psiMI.AR <- rep(0, n.occasions)

psiCC.DB <- rep(0, n.occasions)
psiCC.JB <- rep(0, n.occasions)
psiCC.MI <- rep(0, n.occasions)
psiCC.CC <- c(rep(c(0, 0.20, 0, 0), 9), 0, 0.20)
psiCC.SE <- c(rep(c(0, 0.20, 0.20, 0), 9), 0, 0.20)
psiCC.BR <- c(rep(c(0, 0.35, 0.35, 0), 9), 0, 0.35)
psiCC.AR <- c(rep(c(0, 0, 0.25, 0), 9), 0, 0)

psiSE.DB <- c(rep(c(0, 0, 0, 0.20), 9), 0, 0)
psiSE.JB <- c(rep(c(0.05, 0, 0, 0), 9), 0.5, 0)
psiSE.MI <- c(rep(c(0.05, 0, 0, 0), 9), 0.5, 0)
psiSE.CC <- c(rep(c(0.20, 0.20, 0, 0), 9), 0.20, 0.20)
psiSE.SE <- rep(0.20, n.occasions)
psiSE.BR <- c(rep(c(0, 0, 0.25, 0.25), 9), 0, 0)
psiSE.AR <- c(rep(c(0, 0, 0, 0.05), 9), 0, 0)

psiBR.DB <- c(rep(c(0, 0, 0, 0.3), 9), 0, 0)
psiBR.JB <- rep(0, n.occasions)
psiBR.MI <- rep(0, n.occasions)
psiBR.CC <- rep(0, n.occasions)
psiBR.SE <- c(rep(c(0, 0, 0.5, 0.5), 9), 0, 0)
psiBR.BR <- c(rep(c(0, 0, 0.1, 0), 9), 0, 0)
psiBR.AR <- c(rep(c(0, 0, 0.1, 0), 9), 0, 0)

psiAR.DB <- c(rep(c(0, 0, 0, 0.90), 9), 0, 0)
psiAR.JB <- rep(0, n.occasions)
psiAR.MI <- rep(0, n.occasions)
psiAR.CC <- rep(0, n.occasions)
psiAR.SE <- c(rep(c(0, 0, 0, 0.10), 9), 0, 0)
psiAR.BR <- rep(0, n.occasions)
psiAR.AR <- rep(0, n.occasions)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- c(rep(c(150, 0, 0, 0), 9), 150, 0) # DB
marked[,2] <- c(rep(c(0, 40, 0, 0), 9), 0, 40) # JB
marked[,3] <- c(rep(c(0, 30, 0, 0), 9), 0, 30) # MI
marked[,4] <- c(rep(c(0, 30, 30, 0), 9), 0, 30) # CC
marked[,5] <- c(rep(c(5, 20, 10, 5), 9), 5, 20) # SE
marked[,6] <- c(rep(c(0, 0, 10, 15), 9), 0, 0) # BR
marked[,7] <- c(rep(c(0, 0, 0, 40), 9), 0, 0) # AR
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
      phiDB*psiDB.DB[t], phiDB*psiDB.JB[t], phiDB*psiDB.MI[t], phiDB*psiDB.CC[t], phiDB*psiDB.SE[t], phiDB*psiDB.BR[t], phiDB*psiDB.AR[t], 1-phiDB,
      phiJB*psiJB.DB[t], phiJB*psiJB.JB[t], phiJB*psiJB.MI[t], phiJB*psiJB.CC[t], phiJB*psiJB.SE[t], phiJB*psiJB.BR[t], phiJB*psiJB.AR[t], 1-phiJB,
      phiMI*psiMI.DB[t], phiMI*psiMI.JB[t], phiMI*psiMI.MI[t], phiMI*psiMI.CC[t], phiMI*psiMI.SE[t], phiMI*psiMI.BR[t], phiMI*psiMI.AR[t], 1-phiMI,
      phiCC*psiCC.DB[t], phiCC*psiCC.JB[t], phiCC*psiCC.MI[t], phiCC*psiCC.CC[t], phiCC*psiCC.SE[t], phiCC*psiCC.BR[t], phiCC*psiCC.AR[t], 1-phiCC,
      phiSE*psiSE.DB[t], phiSE*psiSE.JB[t], phiSE*psiSE.MI[t], phiSE*psiSE.CC[t], phiSE*psiSE.SE[t], phiSE*psiSE.BR[t], phiSE*psiSE.AR[t], 1-phiSE,
      phiBR*psiBR.DB[t], phiBR*psiBR.JB[t], phiBR*psiBR.MI[t], phiBR*psiBR.CC[t], phiBR*psiBR.SE[t], phiBR*psiBR.BR[t], phiBR*psiBR.AR[t], 1-phiBR,
      phiAR*psiAR.DB[t], phiAR*psiAR.JB[t], phiAR*psiAR.MI[t], phiAR*psiAR.CC[t], phiAR*psiAR.SE[t], phiAR*psiAR.BR[t], phiAR*psiAR.AR[t], 1-phiAR,
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

sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phi[1]
   phiJB[t] <- mean.phi[2]
   phiMI[t] <- mean.phi[3]
   phiCC[t] <- mean.phi[4]
   phiSE[t] <- mean.phi[5]
   phiBR[t] <- mean.phi[6]
   phiAR[t] <- mean.phi[7]
   }

   for (u in 1:7){
   mean.phi[u] ~ dunif(0, 1)
   }
   
pDB[1] <- 0
pDB[2] <- 0
pDB[3] <- 0
pDB[4] <- mean.pDB
pDB[5] <- 0
pDB[6] <- 0
pDB[7] <- 0
pDB[8] <- mean.pDB
pDB[9] <- 0
pDB[10] <- 0
pDB[11] <- 0
pDB[12] <- mean.pDB
pDB[13] <- 0
pDB[14] <- 0
pDB[15] <- 0
pDB[16] <- mean.pDB
pDB[17] <- 0
pDB[18] <- 0
pDB[19] <- 0
pDB[20] <- mean.pDB
pDB[21] <- 0
pDB[22] <- 0
pDB[23] <- 0
pDB[24] <- mean.pDB
pDB[25] <- 0
pDB[26] <- 0
pDB[27] <- 0
pDB[28] <- mean.pDB
pDB[29] <- 0
pDB[30] <- 0
pDB[31] <- 0
pDB[32] <- mean.pDB
pDB[33] <- 0
pDB[34] <- 0
pDB[35] <- 0
pDB[36] <- mean.pDB
pDB[37] <- 0
pDB[38] <- 0

pJB[1] <- mean.pJB
pJB[2] <- 0
pJB[3] <- 0
pJB[4] <- 0
pJB[5] <- mean.pJB
pJB[6] <- 0
pJB[7] <- 0
pJB[8] <- 0
pJB[9] <- mean.pJB
pJB[10] <- 0
pJB[11] <- 0
pJB[12] <- 0
pJB[13] <- mean.pJB
pJB[14] <- 0
pJB[15] <- 0
pJB[16] <- 0
pJB[17] <- mean.pJB
pJB[18] <- 0
pJB[19] <- 0
pJB[20] <- 0
pJB[21] <- mean.pJB
pJB[22] <- 0
pJB[23] <- 0
pJB[24] <- 0
pJB[25] <- mean.pJB
pJB[26] <- 0
pJB[27] <- 0
pJB[28] <- 0
pJB[29] <- mean.pJB
pJB[30] <- 0
pJB[31] <- 0
pJB[32] <- 0
pJB[33] <- mean.pJB
pJB[34] <- 0
pJB[35] <- 0
pJB[36] <- 0
pJB[37] <- mean.pJB
pJB[38] <- 0

pMI[1] <- mean.pMI
pMI[2] <- 0
pMI[3] <- 0
pMI[4] <- 0
pMI[5] <- mean.pMI
pMI[6] <- 0
pMI[7] <- 0
pMI[8] <- 0
pMI[9] <- mean.pMI
pMI[10] <- 0
pMI[11] <- 0
pMI[12] <- 0
pMI[13] <- mean.pMI
pMI[14] <- 0
pMI[15] <- 0
pMI[16] <- 0
pMI[17] <- mean.pMI
pMI[18] <- 0
pMI[19] <- 0
pMI[20] <- 0
pMI[21] <- mean.pMI
pMI[22] <- 0
pMI[23] <- 0
pMI[24] <- 0
pMI[25] <- mean.pMI
pMI[26] <- 0
pMI[27] <- 0
pMI[28] <- 0
pMI[29] <- mean.pMI
pMI[30] <- 0
pMI[31] <- 0
pMI[32] <- 0
pMI[33] <- mean.pMI
pMI[34] <- 0
pMI[35] <- 0
pMI[36] <- 0
pMI[37] <- mean.pMI
pMI[38] <- 0

pCC[1] <- mean.pCC
pCC[2] <- mean.pCC
pCC[3] <- 0
pCC[4] <- 0
pCC[5] <- mean.pCC
pCC[6] <- mean.pCC
pCC[7] <- 0
pCC[8] <- 0
pCC[9] <- mean.pCC
pCC[10] <- mean.pCC
pCC[11] <- 0
pCC[12] <- 0
pCC[13] <- mean.pCC
pCC[14] <- mean.pCC
pCC[15] <- 0
pCC[16] <- 0
pCC[17] <- mean.pCC
pCC[18] <- mean.pCC
pCC[19] <- 0
pCC[20] <- 0
pCC[21] <- mean.pCC
pCC[22] <- mean.pCC
pCC[23] <- 0
pCC[24] <- 0
pCC[25] <- mean.pCC
pCC[26] <- mean.pCC
pCC[27] <- 0
pCC[28] <- 0
pCC[29] <- mean.pCC
pCC[30] <- mean.pCC
pCC[31] <- 0
pCC[32] <- 0
pCC[33] <- mean.pCC
pCC[34] <- mean.pCC
pCC[35] <- 0
pCC[36] <- 0
pCC[37] <- mean.pCC
pCC[38] <- mean.pCC

pSE[1] <- mean.pSE
pSE[2] <- mean.pSE
pSE[3] <- mean.pSE
pSE[4] <- mean.pSE
pSE[5] <- mean.pSE
pSE[6] <- mean.pSE
pSE[7] <- mean.pSE
pSE[8] <- mean.pSE
pSE[9] <- mean.pSE
pSE[10] <- mean.pSE
pSE[11] <- mean.pSE
pSE[12] <- mean.pSE
pSE[13] <- mean.pSE
pSE[14] <- mean.pSE
pSE[15] <- mean.pSE
pSE[16] <- mean.pSE
pSE[17] <- mean.pSE
pSE[18] <- mean.pSE
pSE[19] <- mean.pSE
pSE[20] <- mean.pSE
pSE[21] <- mean.pSE
pSE[22] <- mean.pSE
pSE[23] <- mean.pSE
pSE[24] <- mean.pSE
pSE[25] <- mean.pSE
pSE[26] <- mean.pSE
pSE[27] <- mean.pSE
pSE[28] <- mean.pSE
pSE[29] <- mean.pSE
pSE[30] <- mean.pSE
pSE[31] <- mean.pSE
pSE[32] <- mean.pSE
pSE[33] <- mean.pSE
pSE[34] <- mean.pSE
pSE[35] <- mean.pSE
pSE[36] <- mean.pSE
pSE[37] <- mean.pSE
pSE[38] <- mean.pSE

pBR[1] <- 0
pBR[2] <- mean.pBR
pBR[3] <- mean.pBR
pBR[4] <- 0
pBR[5] <- 0
pBR[6] <- mean.pBR
pBR[7] <- mean.pBR
pBR[8] <- 0
pBR[9] <- 0
pBR[10] <- mean.pBR
pBR[11] <- mean.pBR
pBR[12] <- 0
pBR[13] <- 0
pBR[14] <- mean.pBR
pBR[15] <- mean.pBR
pBR[16] <- 0
pBR[17] <- 0
pBR[18] <- mean.pBR
pBR[19] <- mean.pBR
pBR[20] <- 0
pBR[21] <- 0
pBR[22] <- mean.pBR
pBR[23] <- mean.pBR
pBR[24] <- 0
pBR[25] <- 0
pBR[26] <- mean.pBR
pBR[27] <- mean.pBR
pBR[28] <- 0
pBR[29] <- 0
pBR[30] <- mean.pBR
pBR[31] <- mean.pBR
pBR[32] <- 0
pBR[33] <- 0
pBR[34] <- mean.pBR
pBR[35] <- mean.pBR
pBR[36] <- 0
pBR[37] <- 0
pBR[38] <- mean.pBR

pAR[1] <- 0
pAR[2] <- 0
pAR[3] <- mean.pAR
pAR[4] <- 0
pAR[5] <- 0
pAR[6] <- 0
pAR[7] <- mean.pAR
pAR[8] <- 0
pAR[9] <- 0
pAR[10] <- 0
pAR[11] <- mean.pAR
pAR[12] <- 0
pAR[13] <- 0
pAR[14] <- 0
pAR[15] <- mean.pAR
pAR[16] <- 0
pAR[17] <- 0
pAR[18] <- 0
pAR[19] <- mean.pAR
pAR[20] <- 0
pAR[21] <- 0
pAR[22] <- 0
pAR[23] <- mean.pAR
pAR[24] <- 0
pAR[25] <- 0
pAR[26] <- 0
pAR[27] <- mean.pAR
pAR[28] <- 0
pAR[29] <- 0
pAR[30] <- 0
pAR[31] <- mean.pAR
pAR[32] <- 0
pAR[33] <- 0
pAR[34] <- 0
pAR[35] <- mean.pAR
pAR[36] <- 0
pAR[37] <- 0
pAR[38] <- 0

mean.pDB ~ dunif(0, 1)
mean.pJB ~ dunif(0, 1)
mean.pMI ~ dunif(0, 1)
mean.pCC ~ dunif(0, 1)
mean.pSE ~ dunif(0, 1)
mean.pBR ~ dunif(0, 1)
mean.pAR ~ dunif(0, 1)

# Transitions: gamma priors
   
   for (i in 1:4){
     db[i] ~ dgamma(1, 1)
   } # i
   
   for (t in 1:(n.occasions-1)){
      psiDB[1,t] <- 0
      psiDB[2,t] <- (db[1]/sum(db[]))*transDB[t]
      psiDB[3,t] <- (db[2]/sum(db[]))*transDB[t]
      psiDB[4,t] <- (db[3]/sum(db[]))*transDB[t]
      psiDB[5,t] <- (db[4]/sum(db[]))*transDB[t]
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
      psiJB[4,t] <- (jb[1]/sum(jb[]))*transJB[t]
      psiJB[5,t] <- (jb[2]/sum(jb[]))*transJB[t]
      psiJB[6,t] <- (jb[3]/sum(jb[]))*transJB[t]
      psiJB[7,t] <- 0
    }
    
    for (i in 1:3){
      mi[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){
      psiMI[1,t] <- 0
      psiMI[2,t] <- 0
      psiMI[3,t] <- 0
      psiMI[4,t] <- (mi[1]/sum(mi[]))*transMI[t]
      psiMI[5,t] <- (mi[2]/sum(mi[]))*transMI[t]
      psiMI[6,t] <- (mi[3]/sum(mi[]))*transMI[t]
      psiMI[7,t] <- 0
    }  
    
    for (i in 1:4){
      cc[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){
      psiCC[1,t] <- 0
      psiCC[2,t] <- 0
      psiCC[3,t] <- 0
      psiCC[4,t] <- (cc[1]/sum(cc[]))*transCCCC[t]
      psiCC[5,t] <- (cc[2]/sum(cc[]))*transCCSE[t]
      psiCC[6,t] <- (cc[3]/sum(cc[]))*transCCBR[t]
      psiCC[7,t] <- (cc[4]/sum(cc[]))*transCCAR[t]
    }
      
    for (i in 1:7){
      se[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){  
      psiSE[1,t] <- (se[1]/sum(se[]))*transSEDB[t]
      psiSE[2,t] <- (se[2]/sum(se[]))*transSEJB[t]
      psiSE[3,t] <- (se[3]/sum(se[]))*transSEMI[t]
      psiSE[4,t] <- (se[4]/sum(se[]))*transSECC[t]
      psiSE[5,t] <- (se[5]/sum(se[]))*transSESE[t]
      psiSE[6,t] <- (se[6]/sum(se[]))*transSEBR[t]
      psiSE[7,t] <- (se[7]/sum(se[]))*transSEAR[t]
    }  
    
    for (i in 1:4){
      br[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){
      psiBR[1,t] <- (br[1]/sum(br[]))*transBRDB[t]
      psiBR[2,t] <- 0
      psiBR[3,t] <- 0
      psiBR[4,t] <- 0
      psiBR[5,t] <- (br[2]/sum(br[]))*transBRSE[t]
      psiBR[6,t] <- (br[3]/sum(br[]))*transBRBR[t]
      psiBR[7,t] <- (br[4]/sum(br[]))*transBRAR[t]
    }  
    
    for (i in 1:2){
      ar[i] ~ dgamma(1, 1)
    }
      
    for (t in 1:(n.occasions-1)){  
      psiAR[1,t] <- (ar[1]/sum(ar[]))*transAR[t]
      psiAR[2,t] <- 0
      psiAR[3,t] <- 0
      psiAR[4,t] <- 0
      psiAR[5,t] <- (ar[2]/sum(ar[]))*transAR[t]
      psiAR[6,t] <- 0
      psiAR[7,t] <- 0
  } # t
   
  # for (i in 1:4){ 
  #   for (s in 1:4){
  #     mean.psiDB[i,s] <- (db[i]/sum(db[]))*transDB[s]
  #   } 
  # }
   
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
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns),
                  ones = diag(ns),
                  transDB = rep(c(1, 0, 0, 0), length.out = n.occasions-1),
                  transJB = rep(c(0, 1, 0, 0), length.out = n.occasions-1),
                  transMI = rep(c(0, 1, 0, 0), length.out = n.occasions-1),
                  transCCCC = rep(c(0, 1, 0, 0), length.out = n.occasions-1),
                  transCCSE = rep(c(0, 1, 1, 0), length.out = n.occasions-1),
                  transCCBR = rep(c(0, 1, 1, 0), length.out = n.occasions-1),
                  transCCAR = rep(c(0, 0, 1, 0), length.out = n.occasions-1),
                  transSEDB = rep(c(0, 0, 0, 1), length.out = n.occasions-1),
                  transSEJB = rep(c(1, 0, 0, 0), length.out = n.occasions-1),
                  transSEMI = rep(c(1, 0, 0, 0), length.out = n.occasions-1),
                  transSECC = rep(c(1, 1, 0, 0), length.out = n.occasions-1),
                  transSESE = rep(c(1, 1, 1, 1), length.out = n.occasions-1),
                  transSEBR = rep(c(0, 1, 1, 0), length.out = n.occasions-1),
                  transSEAR = rep(c(0, 0, 1, 0), length.out = n.occasions-1),
                  transBRDB = rep(c(0, 0, 0, 1), length.out = n.occasions-1),
                  transBRSE = rep(c(0, 0, 1, 1), length.out = n.occasions-1),
                  transBRBR = rep(c(0, 0, 1, 0), length.out = n.occasions-1),
                  transBRAR = rep(c(0, 0, 1, 0), length.out = n.occasions-1),
                  transAR = rep(c(0, 0, 0, 1), length.out = n.occasions-1))

# Initial values 
inits <- function(){list(mean.phi = runif(7, 0, 1),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1), mean.pCC = runif(1, 0, 1),
                         mean.pSE = runif(1, 0, 1), mean.pBR = runif(1, 0, 1), mean.pAR = runif(1, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
                "mean.pDB", "mean.pJB", "mean.pMI", "mean.pCC", "mean.pSE", "mean.pBR", "mean.pAR",
                "fit", "fit.new")

# MCMC settings
ni <- 100000
nt <- 5
nb <- 70000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/scenario-7-marray-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/scenario-7-marray-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/marray-scenario-7/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
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

png(filename = "figures/scenario-7.marray.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                  value = c(phiDB, phiJB, phiMI, phiCC, phiSE, phiBR, phiAR))

# resighting
p.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                value = c(mean.pDB, mean.pJB, mean.pMI, mean.pCC, mean.pSE, mean.pBR, mean.pAR))

# transition
psiDB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI, psiDB.CC, psiDB.SE, psiDB.BR, psiDB.AR))

psiJB.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI, psiJB.CC, psiJB.SE, psiJB.BR, psiJB.AR))

psiMI.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI, psiMI.CC, psiMI.SE, psiMI.BR, psiMI.AR))

psiCC.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiCC.DB, psiCC.JB, psiCC.MI, psiCC.CC, psiCC.SE, psiCC.BR, psiCC.AR))

psiSE.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiSE.DB, psiSE.JB, psiSE.MI, psiSE.CC, psiSE.SE, psiSE.BR, psiSE.AR))

psiBR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiBR.DB, psiBR.JB, psiBR.MI, psiBR.CC, psiBR.SE, psiBR.BR,psiBR.AR))

psiAR.sim <- tibble(site = c("DB", "JB", "MI", "CC", "SE", "BR", "AR"),
                    value = c(psiAR.DB, psiAR.JB, psiAR.MI, psiAR.CC, psiAR.SE, psiAR.BR, psiAR.AR))

# format survival
phi.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3, mean.phi.4, mean.phi.5, mean.phi.6, mean.phi.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phi.1", "DB")) %>%
  mutate(site = str_replace(site, "mean.phi.2", "JB")) %>%
  mutate(site = str_replace(site, "mean.phi.3", "MI")) %>% 
  mutate(site = str_replace(site, "mean.phi.4", "CC")) %>% 
  mutate(site = str_replace(site, "mean.phi.5", "SE")) %>%
  mutate(site = str_replace(site, "mean.phi.6", "BR")) %>%
  mutate(site = str_replace(site, "mean.phi.7", "AR"))

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
  select(mean.pDB, mean.pJB, mean.pMI, mean.pCC, mean.pSE, mean.pBR, mean.pAR) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pDB", "DB")) %>%
  mutate(site = str_replace(site, "mean.pJB", "JB")) %>%
  mutate(site = str_replace(site, "mean.pMI", "MI")) %>% 
  mutate(site = str_replace(site, "mean.pCC", "CC")) %>% 
  mutate(site = str_replace(site, "mean.pSE", "SE")) %>%
  mutate(site = str_replace(site, "mean.pBR", "BR")) %>%
  mutate(site = str_replace(site, "mean.pAR", "AR"))

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
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiDB.1, psiDB.2, psiDB.3, psiDB.4, psiDB.5, psiDB.6, psiDB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiDB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiDB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiDB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiDB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiDB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiDB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiDB.7", "AR"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiJB.1, psiJB.2, psiJB.3, psiJB.4, psiJB.5, psiJB.6, psiJB.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiJB.1", "DB")) %>%
  mutate(site = str_replace(site, "psiJB.2", "JB")) %>%
  mutate(site = str_replace(site, "psiJB.3", "MI")) %>%
  mutate(site = str_replace(site, "psiJB.4", "CC")) %>%
  mutate(site = str_replace(site, "psiJB.5", "SE")) %>%
  mutate(site = str_replace(site, "psiJB.6", "BR")) %>%
  mutate(site = str_replace(site, "psiJB.7", "AR"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiMI.1, psiMI.2, psiMI.3, psiMI.4, psiMI.5, psiMI.6, psiMI.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiMI.1", "DB")) %>%
  mutate(site = str_replace(site, "psiMI.2", "JB")) %>%
  mutate(site = str_replace(site, "psiMI.3", "MI")) %>%
  mutate(site = str_replace(site, "psiMI.4", "CC")) %>%
  mutate(site = str_replace(site, "psiMI.5", "SE")) %>%
  mutate(site = str_replace(site, "psiMI.6", "BR")) %>%
  mutate(site = str_replace(site, "psiMI.7", "AR"))

psiCC.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiCC.1, psiCC.2, psiCC.3, psiCC.4, psiCC.5, psiCC.6, psiCC.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiCC.1", "DB")) %>%
  mutate(site = str_replace(site, "psiCC.2", "JB")) %>%
  mutate(site = str_replace(site, "psiCC.3", "MI")) %>%
  mutate(site = str_replace(site, "psiCC.4", "CC")) %>%
  mutate(site = str_replace(site, "psiCC.5", "SE")) %>%
  mutate(site = str_replace(site, "psiCC.6", "BR")) %>%
  mutate(site = str_replace(site, "psiCC.7", "AR"))

psiSE.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiSE.1, psiSE.2, psiSE.3, psiSE.4, psiSE.5, psiSE.6, psiSE.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiSE.1", "DB")) %>%
  mutate(site = str_replace(site, "psiSE.2", "JB")) %>%
  mutate(site = str_replace(site, "psiSE.3", "MI")) %>%
  mutate(site = str_replace(site, "psiSE.4", "CC")) %>%
  mutate(site = str_replace(site, "psiSE.5", "SE")) %>%
  mutate(site = str_replace(site, "psiSE.6", "BR")) %>%
  mutate(site = str_replace(site, "psiSE.7", "AR"))

psiBR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiBR.1, psiBR.2, psiBR.3, psiBR.4, psiBR.5, psiBR.6, psiBR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiBR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiBR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiBR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiBR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiBR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiBR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiBR.7", "AR"))

psiAR.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(psiAR.1, psiAR.2, psiAR.3, psiAR.4, psiAR.5, psiAR.6, psiAR.7) %>% 
  pivot_longer(cols = 1:7, names_to = "site", values_to = "estimate") %>%
  mutate(site = str_replace(site, "psiAR.1", "DB")) %>%
  mutate(site = str_replace(site, "psiAR.2", "JB")) %>%
  mutate(site = str_replace(site, "psiAR.3", "MI")) %>%
  mutate(site = str_replace(site, "psiAR.4", "CC")) %>%
  mutate(site = str_replace(site, "psiAR.5", "SE")) %>%
  mutate(site = str_replace(site, "psiAR.6", "BR")) %>%
  mutate(site = str_replace(site, "psiAR.7", "AR"))

# plot transitions
psiDB.plot <- ggplot() +
  geom_violin(psiDB.mod,
              mapping = aes(x = site, y = estimate, group = site,
                            fill = site), alpha = 0.6) +
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
  geom_point(psiAR.sim,
             mapping = aes(x = site, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/scenario-7.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/scenario-7.marray.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(p.plot)

dev.off()

# transition
png(filename = "figures/scenario-7.marray.psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/scenario-7.marray.psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/scenario-7.marray.psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

png(filename = "figures/scenario-7.marray.psiCC.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiCC.plot)

dev.off()

png(filename = "figures/scenario-7.marray.psiSE.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiSE.plot)

dev.off()

png(filename = "figures/scenario-7.marray.psiBR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiBR.plot)

dev.off()

png(filename = "figures/scenario-7.marray.psiAR.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiAR.plot)

dev.off()
