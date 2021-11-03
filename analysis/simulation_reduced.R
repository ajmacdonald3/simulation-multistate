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

#### SCENARIO 5 ####

# States: DB, JB, MI, AR
# phi: all constant
# p: constant, following annual cycle
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

pDB <- c(0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0, 0.6, 0, 0)
pJB <- c(0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0, 0, 0.4, 0)
pMI <- c(0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0)
pAR <- c(0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5)

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
marked[,2] <- rep(c(0, 200, 0), n.occasions/3)
marked[,3] <- rep(c(0, 200, 0), n.occasions/3)
marked[,4] <- rep(c(0, 0, 200), n.occasions/3)
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
   
   for (u in 1:4){
   mean.p[u] ~ dunif(0, 1)
   }

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probs
      for (i in 1:3){
         lpsiDB[i] ~ dnorm(0, 0.001)
         lpsiJB[i] ~ dnorm(0, 0.001)
         lpsiMI[i] ~ dnorm(0, 0.001)
         lpsiAR[i] ~ dnorm(0, 0.001)
      }
         
      # Constrain the transitions such that their sum is < 1
         psiDB[1] <- 0
         psiDB[2] <- exp(lpsiDB[2]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]))
         psiDB[3] <- exp(lpsiDB[3]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]))
         
         psiJB[1] <- 0
         psiJB[2] <- 0
         psiJB[3] <- 0
         
         psiMI[1] <- 0
         psiMI[2] <- 0
         psiMI[3] <- 0
         
         psiAR[1] <- exp(lpsiAR[1]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]))
         psiAR[2] <- 0
         psiAR[3] <- 0
         
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
                         mean.p = runif(4, 0, 1),
                         z = mod.init.z2(rCH, 5))} 


# Parameters monitored
parameters <- c("mean.phi",
                "psiDB", "psiJB", "psiMI", "psiAR",
                "mean.p")

# MCMC settings
ni <- 100000
nt <- 10
nb <- 50000
nc <- 3

# Call JAGS from R (BRT 56 min)
rekn.ms <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(rekn.ms, digits = 3)

rekn.ms.update <- autojags(jags.data, inits, parameters, "simulation_scenarios.jags",
                           n.chains = nc, iter.increment = 50000, n.burnin = nb,
                           n.thin = nt, save.all.iter = FALSE, Rhat.limit = 1.2,
                           max.iter = 200000, parallel = TRUE)

print(rekn.ms.update, digits = 3)

saveRDS(rekn.ms.update$summary, file = paste0("./analysis-output/ms-simulation-summary", Sys.Date(), ".rds"))

saveRDS(rekn.ms.update$sims.list, file = paste0("./analysis-output/ms-simulation-simslist", Sys.Date(), ".rds"))

################################################################################
# PLOT SIMULATION RESULTS
#
################################################################################
library(tidyverse)
library(cowplot)
library(viridis)

# parameter identifiability checks
sims.list <- readRDS("analysis-output/ms-simulation-simslist2021-02-16.rds")
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
   
   png(filename = paste0("analysis-output/parameter-identifiability/2021-02-16/",
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
psi.sim <- tibble(transition = c("DB-JB", "DB-MI", "JB-AR", "MI-AR", "AR-DB"),
                  value = c(0.5, 0.5, 1, 1, 1))

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


# p.plot <- ggplot() +
#    geom_violin(p.mod, mapping = aes(x = as.numeric(occasion), y = estimate,
#                                     group = occasion, fill = occasion), alpha = 0.6) +
#    geom_point(p.sim, mapping = aes(x = as.numeric(occasion), y = value),
#               size = 2, colour = "black") +
#    scale_fill_viridis(discrete = TRUE) +
#    scale_x_continuous(breaks = seq(1, 14, 1)) +
#    #coord_flip() +
#    facet_wrap(site ~ ., nrow = 4) +
#    xlab("Occasion") +
#    ylab("Resighting probability") +
#    theme(legend.position = "none")

# format transitions
psi.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiDB.2, psiDB.3, psiJB.4, psiMI.4, psiAR.1) %>% 
   pivot_longer(cols = 1:5, names_to = "transition", values_to = "estimate") %>% 
   mutate(transition = str_replace(transition, "psiDB.2", "DB-JB")) %>%
   mutate(transition = str_replace(transition, "psiDB.3", "DB-MI")) %>%
   mutate(transition = str_replace(transition, "psiJB.4", "JB-AR")) %>%
   mutate(transition = str_replace(transition, "psiMI.4", "MI-AR")) %>%
   mutate(transition = str_replace(transition, "psiAR.1", "AR-DB"))

# plot transitions
psiDB.plot <- ggplot() +
   geom_violin(psi.mod %>% filter(transition %in% c("DB-JB", "DB-MI")),
               mapping = aes(x = transition, y = estimate, group = transition,
                             fill = transition), alpha = 0.6) +
   geom_point(psi.sim %>% filter(transition %in% c("DB-JB", "DB-MI")),
              mapping = aes(x = transition, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Transition") +
   ylab("Transition probability") +
   theme(legend.position = "none")

psi.plot <- ggplot() +
   geom_violin(psi.mod %>% filter(!transition %in% c("DB-JB", "DB-MI")),
               mapping = aes(x = transition, y = estimate, group = transition,
                             fill = transition), alpha = 0.6) +
   geom_point(psi.sim %>% filter(!transition %in% c("DB-JB", "DB-MI")),
              mapping = aes(x = transition, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Transition") +
   ylab("Transition probability") +
   theme(legend.position = "none")

# plot all in grid
sim.plots <- plot_grid(phi.plot, p.plot, psiDB.plot, psi.plot, labels = c("A", "B", "C", "D"), ncol = 2)

png(filename = "figures/plot.sim.reduced.results.png", width = 8, height = 8,
    units = "in", res = 600)

print(sim.plots)

dev.off()

# plot_grid(bottom.row, p.plot, nrow = 2, labels = c("", "D"), rel_heights = c(2/5, 3/5))


