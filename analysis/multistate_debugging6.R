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
  for (i in 1:(dim(CH)[1])){
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
# Transitions: constant with JB-JB set to 0

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 40
n.years <- 10
n.states <- 4
n.obs <- 4

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

# transition probabilities (sum to 1)
psiDB.DB <- rep(c(0.3, 0.5, 0.6, 0.3), length.out = n.occasions-1)
psiDB.JB <- rep(c(0.2, 0.1, 0.2, 0.6), length.out = n.occasions-1)
psiDB.MI <- rep(c(0.5, 0.4, 0.2, 0.1), length.out = n.occasions-1)

psiJB.DB <- rep(c(0.4, 0.7, 0.3, 0.2), length.out = n.occasions-1)
psiJB.JB <- rep(c(0.5, 0.1, 0.4, 0.1), length.out = n.occasions-1)
psiJB.MI <- rep(c(0.1, 0.2, 0.3, 0.7), length.out = n.occasions-1)

psiMI.DB <- rep(c(0.1, 0.2, 0.5, 0.4), length.out = n.occasions-1)
psiMI.JB <- rep(c(0.4, 0.7, 0.2, 0.1), length.out = n.occasions-1)
psiMI.MI <- rep(c(0.5, 0.1, 0.3, 0.5), length.out = n.occasions-1)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(30, n.occasions) # DB
marked[,2] <- rep(30, n.occasions) # JB
marked[,3] <- rep(30, n.occasions) # MI
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
      phiDB[t]*psiDB.DB[t], phiDB[t]*psiDB.JB[t], phiDB[t]*psiDB.MI[t], 1-phiDB[t],
      phiJB[t]*psiJB.DB[t], phiJB[t]*psiJB.JB[t], phiJB[t]*psiJB.MI[t], 1-phiJB[t],
      phiMI[t]*psiMI.DB[t], phiMI[t]*psiMI.JB[t], phiMI[t]*psiMI.MI[t], 1-phiMI[t],
      0,                    0,                    0,                    1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pDB[t], 0,      0,      1-pDB[t],
      0,      pJB[t], 0,      1-pJB[t],
      0,      0,      pMI[t], 1-pMI[t],
      0,      0,      0,      1),
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
  
  # Transitions: gamma priors
  for (i in 1:3){
  for (s in 1:4){
    db[i,s] ~ dgamma(1, 1)
    jb[i,s] ~ dgamma(1, 1)
    mi[i,s] ~ dgamma(1, 1)
    }
  }
  
  for (s in 1:4){
    intDB1[s] <- (db[1,s]/sum(db[,s])) 
    intDB2[s] <- (db[2,s]/sum(db[,s]))
    intDB3[s] <- (db[3,s]/sum(db[,s]))
    
    intJB1[s] <- (jb[1,s]/sum(jb[,s]))
    intJB2[s] <- (jb[2,s]/sum(jb[,s]))
    intJB3[s] <- (jb[3,s]/sum(jb[,s]))
    
    intMI1[s] <- (mi[1,s]/sum(mi[,s]))
    intMI2[s] <- (mi[2,s]/sum(mi[,s]))
    intMI3[s] <- (mi[3,s]/sum(mi[,s]))
  }
  
  for (t in 1:(n.occasions-1)){  
    psiDB[1,t] <- intDB1[season[t]]
    psiDB[2,t] <- intDB2[season[t]]
    psiDB[3,t] <- intDB3[season[t]]
    
    psiJB[1,t] <- intJB1[season[t]]
    psiJB[2,t] <- intJB2[season[t]]
    psiJB[3,t] <- intJB3[season[t]]
    
    psiMI[1,t] <- intMI1[season[t]]
    psiMI[2,t] <- intMI2[season[t]]
    psiMI[3,t] <- intMI3[season[t]]
  } 
  
  # Define state-transition and reencounter probabilities - note no i index - is no longer individual 
  for (t in 1:(n.occasions-1)){
    psi[1,t,1] <- phiDB[t] * psiDB[1,t]
    psi[1,t,2] <- phiDB[t] * psiDB[2,t]
    psi[1,t,3] <- phiDB[t] * psiDB[3,t]

    psi[2,t,1] <- phiJB[t] * psiJB[1,t]
    psi[2,t,2] <- phiJB[t] * psiJB[2,t]
    psi[2,t,3] <- phiJB[t] * psiJB[3,t]

    psi[3,t,1] <- phiMI[t] * psiMI[1,t]
    psi[3,t,2] <- phiMI[t] * psiMI[2,t]
    psi[3,t,3] <- phiMI[t] * psiMI[3,t]

    po[1,t] <- pDB[t]
    po[2,t] <- pJB[t]
    po[3,t] <- pMI[t]
    
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
                         mean.pDB = runif(4, 0, 1), mean.pJB = runif(4, 0, 1), mean.pMI = runif(4, 0, 1))}  


# Parameters monitored
parameters <- c("mean.phiDB", "mean.phiJB", "mean.phiMI",
                "mean.pDB", "mean.pJB", "mean.pMI",
                "intDB1", "intDB2", "intDB3", "intJB1", "intJB2", "intJB3", "intMI1", "intMI2", "intMI3",
                "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 5
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 2 days)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/debugging/multistate-debugging6-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/debugging/multistate-debugging6-simslist", Sys.Date(), ".rds"))

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/debugging6/",
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

png(filename = "figures/debugging/debugging6/ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phiDB.sim <- tibble(site = c("DB-1", "DB-2", "DB-3", "DB-4"),
                    value = c(mean.phiDB1, mean.phiDB2, mean.phiDB3, mean.phiDB4))

phiJB.sim <- tibble(site = c("JB-1", "JB-2", "JB-3", "JB-4"),
                    value = c(mean.phiJB1, mean.phiJB2, mean.phiJB3, mean.phiJB4))

phiMI.sim <- tibble(site = c("MI-1", "MI-2", "MI-3", "MI-4"),
                    value = c(mean.phiMI1, mean.phiMI2, mean.phiMI3, mean.phiMI4))

# resighting
pDB.sim <- tibble(site = c("DB-1", "DB-2", "DB-3", "DB-4"),
                  value = c(mean.pDB1, mean.pDB2, mean.pDB3, mean.pDB4))

pJB.sim <- tibble(site = c("JB-1", "JB-2", "JB-3", "JB-4"),
                  value = c(mean.pJB1, mean.pJB2, mean.pJB3, mean.pJB4))

pMI.sim <- tibble(site = c("MI-1", "MI-2", "MI-3", "MI-4"),
                  value = c(mean.pMI1, mean.pMI2, mean.pMI3, mean.pMI4))

# transitions
psiDB.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 12),
                    transition = c(rep("DB-DB", 4), rep("DB-JB", 4), rep("DB-MI", 4)),
                    value = c(psiDB.DB[1:4], psiDB.JB[1:4], psiDB.MI[1:4]))

psiJB.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 12),
                    transition = c(rep("JB-DB", 4), rep("JB-JB", 4), rep("JB-MI", 4)),
                    value = c(psiJB.DB[1:4], psiJB.JB[1:4], psiJB.MI[1:4]))

psiMI.sim <- tibble(season = rep(c(1, 2, 3, 4), length.out = 12),
                    transition = c(rep("MI-DB", 4), rep("MI-JB", 4), rep("MI-MI", 4)),
                    value = c(psiMI.DB[1:4], psiMI.JB[1:4], psiMI.MI[1:4]))

# format survival
phiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiDB.1, mean.phiDB.2, mean.phiDB.3, mean.phiDB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiDB.1", "DB-1")) %>%
  mutate(site = str_replace(site, "mean.phiDB.2", "DB-2")) %>%
  mutate(site = str_replace(site, "mean.phiDB.3", "DB-3")) %>% 
  mutate(site = str_replace(site, "mean.phiDB.4", "DB-4"))

phiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiJB.1, mean.phiJB.2, mean.phiJB.3, mean.phiJB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiJB.1", "JB-1")) %>%
  mutate(site = str_replace(site, "mean.phiJB.2", "JB-2")) %>%
  mutate(site = str_replace(site, "mean.phiJB.3", "JB-3")) %>% 
  mutate(site = str_replace(site, "mean.phiJB.4", "JB-4"))

phiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.phiMI.1, mean.phiMI.2, mean.phiMI.3, mean.phiMI.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.phiMI.1", "MI-1")) %>%
  mutate(site = str_replace(site, "mean.phiMI.2", "MI-2")) %>%
  mutate(site = str_replace(site, "mean.phiMI.3", "MI-3")) %>% 
  mutate(site = str_replace(site, "mean.phiMI.4", "MI-4"))

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

# format resighting
pDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pDB.1, mean.pDB.2, mean.pDB.3, mean.pDB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pDB.1", "DB-1")) %>%
  mutate(site = str_replace(site, "mean.pDB.2", "DB-2")) %>%
  mutate(site = str_replace(site, "mean.pDB.3", "DB-3")) %>% 
  mutate(site = str_replace(site, "mean.pDB.4", "DB-4"))

pJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pJB.1, mean.pJB.2, mean.pJB.3, mean.pJB.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pJB.1", "JB-1")) %>%
  mutate(site = str_replace(site, "mean.pJB.2", "JB-2")) %>%
  mutate(site = str_replace(site, "mean.pJB.3", "JB-3")) %>% 
  mutate(site = str_replace(site, "mean.pJB.4", "JB-4"))

pMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(mean.pMI.1, mean.pMI.2, mean.pMI.3, mean.pMI.4) %>% 
  pivot_longer(cols = 1:4, names_to = "site", values_to = "estimate") %>% 
  mutate(site = str_replace(site, "mean.pMI.1", "MI-1")) %>%
  mutate(site = str_replace(site, "mean.pMI.2", "MI-2")) %>%
  mutate(site = str_replace(site, "mean.pMI.3", "MI-3")) %>% 
  mutate(site = str_replace(site, "mean.pMI.4", "MI-4"))

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
  ylab("Resighting probability") +
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
  ylab("Resighting probability") +
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
  ylab("Resighting probability") +
  theme(legend.position = "none")

# format transitions
psiDB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(intDB1.1, intDB1.2, intDB1.3, intDB1.4, intDB2.1, intDB2.2, intDB2.3,
         intDB2.4, intDB3.1, intDB3.2, intDB3.3, intDB3.4) %>% 
  pivot_longer(cols = 1:12, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "intDB1.1" ~ "1",
                            transition == "intDB1.2" ~ "2",
                            transition == "intDB1.3" ~ "3",
                            transition == "intDB1.4" ~ "4",
                            transition == "intDB2.1" ~ "1",
                            transition == "intDB2.2" ~ "2",
                            transition == "intDB2.3" ~ "3",
                            transition == "intDB2.4" ~ "4",
                            transition == "intDB3.1" ~ "1",
                            transition == "intDB3.2" ~ "2",
                            transition == "intDB3.3" ~ "3",
                            transition == "intDB3.4" ~ "4")) %>% 
  mutate(transition = str_replace(transition, "intDB1.1", "DB-DB")) %>%
  mutate(transition = str_replace(transition, "intDB1.2", "DB-DB")) %>%
  mutate(transition = str_replace(transition, "intDB1.3", "DB-DB")) %>% 
  mutate(transition = str_replace(transition, "intDB1.4", "DB-DB")) %>% 
  mutate(transition = str_replace(transition, "intDB2.1", "DB-JB")) %>%
  mutate(transition = str_replace(transition, "intDB2.2", "DB-JB")) %>%
  mutate(transition = str_replace(transition, "intDB2.3", "DB-JB")) %>% 
  mutate(transition = str_replace(transition, "intDB2.4", "DB-JB")) %>% 
  mutate(transition = str_replace(transition, "intDB3.1", "DB-MI")) %>%
  mutate(transition = str_replace(transition, "intDB3.2", "DB-MI")) %>%
  mutate(transition = str_replace(transition, "intDB3.3", "DB-MI")) %>% 
  mutate(transition = str_replace(transition, "intDB3.4", "DB-MI"))

psiJB.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(intJB1.1, intJB1.2, intJB1.3, intJB1.4, intJB2.1, intJB2.2, intJB2.3,
         intJB2.4, intJB3.1, intJB3.2, intJB3.3, intJB3.4) %>% 
  pivot_longer(cols = 1:12, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "intJB1.1" ~ "1",
                            transition == "intJB1.2" ~ "2",
                            transition == "intJB1.3" ~ "3",
                            transition == "intJB1.4" ~ "4",
                            transition == "intJB2.1" ~ "1",
                            transition == "intJB2.2" ~ "2",
                            transition == "intJB2.3" ~ "3",
                            transition == "intJB2.4" ~ "4",
                            transition == "intJB3.1" ~ "1",
                            transition == "intJB3.2" ~ "2",
                            transition == "intJB3.3" ~ "3",
                            transition == "intJB3.4" ~ "4")) %>% 
  mutate(transition = str_replace(transition, "intJB1.1", "JB-DB")) %>%
  mutate(transition = str_replace(transition, "intJB1.2", "JB-DB")) %>%
  mutate(transition = str_replace(transition, "intJB1.3", "JB-DB")) %>% 
  mutate(transition = str_replace(transition, "intJB1.4", "JB-DB")) %>% 
  mutate(transition = str_replace(transition, "intJB2.1", "JB-JB")) %>%
  mutate(transition = str_replace(transition, "intJB2.2", "JB-JB")) %>%
  mutate(transition = str_replace(transition, "intJB2.3", "JB-JB")) %>% 
  mutate(transition = str_replace(transition, "intJB2.4", "JB-JB")) %>% 
  mutate(transition = str_replace(transition, "intJB3.1", "JB-MI")) %>%
  mutate(transition = str_replace(transition, "intJB3.2", "JB-MI")) %>%
  mutate(transition = str_replace(transition, "intJB3.3", "JB-MI")) %>% 
  mutate(transition = str_replace(transition, "intJB3.4", "JB-MI"))

psiMI.mod <- sims.list %>% 
  as.data.frame() %>% 
  select(intMI1.1, intMI1.2, intMI1.3, intMI1.4, intMI2.1, intMI2.2, intMI2.3,
         intMI2.4, intMI3.1, intMI3.2, intMI3.3, intMI3.4) %>% 
  pivot_longer(cols = 1:12, names_to = "transition", values_to = "estimate") %>%
  mutate(season = case_when(transition == "intMI1.1" ~ "1",
                            transition == "intMI1.2" ~ "2",
                            transition == "intMI1.3" ~ "3",
                            transition == "intMI1.4" ~ "4",
                            transition == "intMI2.1" ~ "1",
                            transition == "intMI2.2" ~ "2",
                            transition == "intMI2.3" ~ "3",
                            transition == "intMI2.4" ~ "4",
                            transition == "intMI3.1" ~ "1",
                            transition == "intMI3.2" ~ "2",
                            transition == "intMI3.3" ~ "3",
                            transition == "intMI3.4" ~ "4")) %>% 
  mutate(transition = str_replace(transition, "intMI1.1", "MI-DB")) %>%
  mutate(transition = str_replace(transition, "intMI1.2", "MI-DB")) %>%
  mutate(transition = str_replace(transition, "intMI1.3", "MI-DB")) %>% 
  mutate(transition = str_replace(transition, "intMI1.4", "MI-DB")) %>% 
  mutate(transition = str_replace(transition, "intMI2.1", "MI-JB")) %>%
  mutate(transition = str_replace(transition, "intMI2.2", "MI-JB")) %>%
  mutate(transition = str_replace(transition, "intMI2.3", "MI-JB")) %>% 
  mutate(transition = str_replace(transition, "intMI2.4", "MI-JB")) %>% 
  mutate(transition = str_replace(transition, "intMI3.1", "MI-MI")) %>%
  mutate(transition = str_replace(transition, "intMI3.2", "MI-MI")) %>%
  mutate(transition = str_replace(transition, "intMI3.3", "MI-MI")) %>% 
  mutate(transition = str_replace(transition, "intMI3.4", "MI-MI"))


# plot transitions
psiDB.plot <- ggplot() +
  geom_violin(psiDB.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiDB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiDB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiDB.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiJB.plot <- ggplot() +
  geom_violin(psiJB.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiJB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiJB.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiJB.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

psiMI.plot <- ggplot() +
  geom_violin(psiMI.mod,
              mapping = aes(x = transition, y = estimate, group = transition,
                            fill = transition), alpha = 0.6) +
  geom_boxplot(psiMI.mod, mapping = aes(x = transition, y = estimate, group = transition),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(psiMI.mod, mapping = aes(x = transition, y = estimate, group = transition),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(psiMI.sim,
             mapping = aes(x = transition, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(. ~ season, ncol = 2) +
  xlab("Transition") +
  ylab("Transition probability") +
  theme(legend.position = "none")

# plot all
# survival
png(filename = "figures/debugging/debugging6/phiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiDB.plot)

dev.off()

png(filename = "figures/debugging/debugging6/phiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiJB.plot)

dev.off()

png(filename = "figures/debugging/debugging6/phiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiMI.plot)

dev.off()

# resighting
png(filename = "figures/debugging/debugging6/pDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pDB.plot)

dev.off()

png(filename = "figures/debugging/debugging6/pJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pJB.plot)

dev.off()

png(filename = "figures/debugging/debugging6/pMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(pMI.plot)

dev.off()

# transition
png(filename = "figures/debugging/debugging6/psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/debugging/debugging6/psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/debugging/debugging6/psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()
