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

#### SCENARIO 5 ####

# States: DB, JB, MI, AR
# phi: all constant
# p: constant, following annual cycle
# psi: following annual cycle

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

pDB <- c(0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0)
pJB <- c(0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0)
pMI <- c(0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0)
pAR <- c(0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5)

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
marked[,1] <- rep(c(700, 0, 0), n.occasions/3)
marked[,2] <- rep(c(0, 300, 0), n.occasions/3)
marked[,3] <- rep(c(0, 300, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 300), n.occasions/3)
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
#CH.TRUE <- sim$CH.TRUE

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

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

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
   }

   for (u in 1:4){
   mean.phi[u] ~ dunif(0, 1)
   }
   
   pDB[1] <- 0
   pDB[2] <- 0
   pDB[3] <- mean.p[1]
   pDB[4] <- 0
   pDB[5] <- 0
   pDB[6] <- mean.p[1]
   pDB[7] <- 0
   pDB[8] <- 0
   pDB[9] <- mean.p[1]
   pDB[10] <- 0
   pDB[11] <- 0
   pDB[12] <- mean.p[1]
   pDB[13] <- 0
   pDB[14] <- 0
   pDB[15] <- mean.p[1]
   pDB[16] <- 0
   pDB[17] <- 0
   pDB[18] <- mean.p[1]
   pDB[19] <- 0
   pDB[20] <- 0
   pDB[21] <- mean.p[1]
   pDB[22] <- 0
   pDB[23] <- 0
   pDB[24] <- mean.p[1]
   pDB[25] <- 0
   pDB[26] <- 0
   pDB[27] <- mean.p[1]
   pDB[28] <- 0
   pDB[29] <- 0
   
   pJB[1] <- mean.p[2]
   pJB[2] <- 0
   pJB[3] <- 0
   pJB[4] <- mean.p[2]
   pJB[5] <- 0
   pJB[6] <- 0
   pJB[7] <- mean.p[2]
   pJB[8] <- 0
   pJB[9] <- 0
   pJB[10] <- mean.p[2]
   pJB[11] <- 0
   pJB[12] <- 0
   pJB[13] <- mean.p[2]
   pJB[14] <- 0
   pJB[15] <- 0
   pJB[16] <- mean.p[2]
   pJB[17] <- 0
   pJB[18] <- 0
   pJB[19] <- mean.p[2]
   pJB[20] <- 0
   pJB[21] <- 0
   pJB[22] <- mean.p[2]
   pJB[23] <- 0
   pJB[24] <- 0
   pJB[25] <- mean.p[2]
   pJB[26] <- 0
   pJB[27] <- 0
   pJB[28] <- mean.p[2]
   pJB[29] <- 0
   
   pMI[1] <- mean.p[3]
   pMI[2] <- 0
   pMI[3] <- 0
   pMI[4] <- mean.p[3]
   pMI[5] <- 0
   pMI[6] <- 0
   pMI[7] <- mean.p[3]
   pMI[8] <- 0
   pMI[9] <- 0
   pMI[10] <- mean.p[3]
   pMI[11] <- 0
   pMI[12] <- 0
   pMI[13] <- mean.p[3]
   pMI[14] <- 0
   pMI[15] <- 0
   pMI[16] <- mean.p[3]
   pMI[17] <- 0
   pMI[18] <- 0
   pMI[19] <- mean.p[3]
   pMI[20] <- 0
   pMI[21] <- 0
   pMI[22] <- mean.p[3]
   pMI[23] <- 0
   pMI[24] <- 0
   pMI[25] <- mean.p[3]
   pMI[26] <- 0
   pMI[27] <- 0
   pMI[28] <- mean.p[3]
   pMI[29] <- 0
   
   pAR[1] <- 0
   pAR[2] <- mean.p[4]
   pAR[3] <- 0
   pAR[4] <- 0
   pAR[5] <- mean.p[4]
   pAR[6] <- 0
   pAR[7] <- 0
   pAR[8] <- mean.p[4]
   pAR[9] <- 0
   pAR[10] <- 0
   pAR[11] <- mean.p[4]
   pAR[12] <- 0
   pAR[13] <- 0
   pAR[14] <- mean.p[4]
   pAR[15] <- 0
   pAR[16] <- 0
   pAR[17] <- mean.p[4]
   pAR[18] <- 0
   pAR[19] <- 0
   pAR[20] <- mean.p[4]
   pAR[21] <- 0
   pAR[22] <- 0
   pAR[23] <- mean.p[4]
   pAR[24] <- 0
   pAR[25] <- 0
   pAR[26] <- mean.p[4]
   pAR[27] <- 0
   pAR[28] <- 0
   pAR[29] <- mean.p[4]
   
   for (u in 1:4){
   mean.p[u] ~ dunif(0, 1)
   }

# Transitions: gamma priors
   
   for (i in 1:4){
      db[i] ~ dgamma(1, 1)
   }
   
   psiDB[1] <- 0
   psiDB[2] <- db[2]/sum(db[])
   psiDB[3] <- db[3]/sum(db[])
   psiDB[4] <- 0
   
   psiJB[1] <- 0
   psiJB[2] <- 0
   psiJB[3] <- 0
      #jb ~ dgamma(1, 1)
   psiJB[4] <- 1
      
   psiMI[1] <- 0
   psiMI[2] <- 0
   psiMI[3] <- 0
      #mi ~ dgamma(1, 1)
   psiMI[4] <- 1
   
      #ar ~ dgamma(1, 1)
   psiAR[1] <- 1
   psiAR[2] <- 0
   psiAR[3] <- 0
   psiAR[4] <- 0

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- 0
      psi[1,t,2] <- phiDB[t] * psiDB[2]
      psi[1,t,3] <- phiDB[t] * psiDB[3]
      psi[1,t,4] <- 0
      
      psi[2,t,1] <- 0
      psi[2,t,2] <- 0
      psi[2,t,3] <- 0
      psi[2,t,4] <- phiJB[t] * psiJB[4]

      psi[3,t,1] <- 0
      psi[3,t,2] <- 0
      psi[3,t,3] <- 0
      psi[3,t,4] <- phiMI[t] * psiMI[4]
      
      psi[4,t,1] <- phiAR[t] * psiAR[1]
      psi[4,t,2] <- 0
      psi[4,t,3] <- 0
      psi[4,t,4] <- 0
      
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
parameters <- c("mean.phi", "psiDB", "psiJB", "psiMI", "psiAR", "mean.p")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 50000
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/ms-simulation-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/ms-simulation-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################
library(tidyverse)
library(cowplot)
library(viridis)

# parameter identifiability checks
sims.list <- readRDS("analysis-output/ms-simulation-simslist2021-02-17.rds")
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
   
   png(filename = paste0("analysis-output/parameter-identifiability/2021-02-17-marray/",
                         i,"-","check.png"),
       width=4, height=3, units="in", res=600)
   
   print(ggplot(sims.list, aes(sims.list[,i])) +
            geom_density() +
            geom_hline(yintercept = 1, linetype = "dashed") +
            xlab(i))
   
   dev.off()
   
}

# create dataframes of set values
# survival
phi.sim <- tibble(site = c("DB", "JB", "MI", "AR"),
                  value = c(phiDB, phiJB, phiMI, phiAR))

# resighting
p.sim <- tibble(site = c("DB", "JB", "MI", "AR"),
                value = c(0.6, 0.4, 0.5, 0.5))

# transition
psi.sim <- tibble(transition = c("DB-JB", "DB-MI"),
                  value = c(0.5, 0.5))

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
   select(psiDB.2, psiDB.3) %>% 
   pivot_longer(cols = 1:2, names_to = "transition", values_to = "estimate") %>% 
   mutate(transition = str_replace(transition, "psiDB.2", "DB-JB")) %>%
   mutate(transition = str_replace(transition, "psiDB.3", "DB-MI"))

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

# plot all
# survival
png(filename = "figures/sim.reduced.marray.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(phi.plot)

dev.off()

# resighting
png(filename = "figures/sim.reduced.marray.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(p.plot)

dev.off()

# transition
png(filename = "figures/sim.reduced.marray.psi.png", width = 8, height = 8,
    units = "in", res = 600)

print(psi.plot)

dev.off()
