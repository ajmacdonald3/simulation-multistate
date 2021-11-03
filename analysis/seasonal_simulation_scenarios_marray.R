################################################################################
# MULTISTATE SEASONAL SURVIVAL SIMULATION SCENARIOS
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
# States: DB, JB, MI
# phi: seasonal with temporal random effect
# p: seasonal with temporal random effect
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 40
n.years <- 10
n.states <- 4
n.obs <- 4

mean.phiDB1 <- 0.9
mean.phiDB2 <- 0.7
mean.phiDB3 <- 0.8
mean.phiDB4 <- 0.6
var.phiDB <- 0.3                       # Temporal variance of survival

mean.pDB1 <- 0.7
mean.pDB2 <- 0.7
mean.pDB3 <- 0.7
mean.pDB4 <- 0.7
var.pDB <- 0.5

# Determine annual survival probabilities
logit.phiDB1 <- rnorm(n.years, qlogis(mean.phiDB1), var.phiDB^0.5)
phiDB1 <- plogis(logit.phiDB1)

logit.phiDB2 <- rnorm(n.years, qlogis(mean.phiDB2), var.phiDB^0.5)
phiDB2 <- plogis(logit.phiDB2)

logit.phiDB3 <- rnorm(n.years, qlogis(mean.phiDB3), var.phiDB^0.5)
phiDB3 <- plogis(logit.phiDB3)

logit.phiDB4 <- rnorm(n.years, qlogis(mean.phiDB4), var.phiDB^0.5)
phiDB4 <- plogis(logit.phiDB4)

phiDB <- c(rbind(phiDB1, phiDB2, phiDB3, phiDB4))
phiDB <- phiDB[1:39]

# Determine annual resighting probabilities
logit.pDB1 <- rnorm(n.years, qlogis(mean.pDB1), var.pDB^0.5)
pDB1 <- plogis(logit.pDB1)

logit.pDB2 <- rnorm(n.years, qlogis(mean.pDB2), var.pDB^0.5)
pDB2 <- plogis(logit.pDB2)

logit.pDB3 <- rnorm(n.years, qlogis(mean.pDB3), var.pDB^0.5)
pDB3 <- plogis(logit.pDB3)

logit.pDB4 <- rnorm(n.years, qlogis(mean.pDB4), var.pDB^0.5)
pDB4 <- plogis(logit.pDB4)

pDB <- c(rbind(pDB1, pDB2, pDB3, pDB4))
pDB <- pDB[1:39]


mean.phiJB1 <- 0.5
mean.phiJB2 <- 0.8
mean.phiJB3 <- 0.6
mean.phiJB4 <- 0.7
var.phiJB <- 0.3                       # Temporal variance of survival

mean.pJB1 <- 0.4
mean.pJB2 <- 0.4
mean.pJB3 <- 0.4
mean.pJB4 <- 0.4
var.pJB <- 0.5

# Determine annual survival probabilities
logit.phiJB1 <- rnorm(n.years, qlogis(mean.phiJB1), var.phiJB^0.5)
phiJB1 <- plogis(logit.phiJB1)

logit.phiJB2 <- rnorm(n.years, qlogis(mean.phiJB2), var.phiJB^0.5)
phiJB2 <- plogis(logit.phiJB2)

logit.phiJB3 <- rnorm(n.years, qlogis(mean.phiJB3), var.phiJB^0.5)
phiJB3 <- plogis(logit.phiJB3)

logit.phiJB4 <- rnorm(n.years, qlogis(mean.phiJB4), var.phiJB^0.5)
phiJB4 <- plogis(logit.phiJB4)

phiJB <- c(rbind(phiJB1, phiJB2, phiJB3, phiJB4))
phiJB <- phiJB[1:39]

# Determine annual resighting probabilities
logit.pJB1 <- rnorm(n.years, qlogis(mean.pJB1), var.pJB^0.5)
pJB1 <- plogis(logit.pJB1)

logit.pJB2 <- rnorm(n.years, qlogis(mean.pJB2), var.pJB^0.5)
pJB2 <- plogis(logit.pJB2)

logit.pJB3 <- rnorm(n.years, qlogis(mean.pJB3), var.pJB^0.5)
pJB3 <- plogis(logit.pJB3)

logit.pJB4 <- rnorm(n.years, qlogis(mean.pJB4), var.pJB^0.5)
pJB4 <- plogis(logit.pJB4)

pJB <- c(rbind(pJB1, pJB2, pJB3, pJB4))
pJB <- pJB[1:39]

mean.phiMI1 <- 0.6
mean.phiMI2 <- 0.8
mean.phiMI3 <- 0.7
mean.phiMI4 <- 0.5
var.phiMI <- 0.3                       # Temporal variance of survival

mean.pMI1 <- 0.5
mean.pMI2 <- 0.5
mean.pMI3 <- 0.5
mean.pMI4 <- 0.5
var.pMI <- 0.5

# Determine annual survival probabilities
logit.phiMI1 <- rnorm(n.years, qlogis(mean.phiMI1), var.phiMI^0.5)
phiMI1 <- plogis(logit.phiMI1)

logit.phiMI2 <- rnorm(n.years, qlogis(mean.phiMI2), var.phiMI^0.5)
phiMI2 <- plogis(logit.phiMI2)

logit.phiMI3 <- rnorm(n.years, qlogis(mean.phiMI3), var.phiMI^0.5)
phiMI3 <- plogis(logit.phiMI3)

logit.phiMI4 <- rnorm(n.years, qlogis(mean.phiMI4), var.phiMI^0.5)
phiMI4 <- plogis(logit.phiMI4)

phiMI <- c(rbind(phiMI1, phiMI2, phiMI3, phiMI4))
phiMI <- phiMI[1:39]

# Determine annual resighting probabilities
logit.pMI1 <- rnorm(n.years, qlogis(mean.pMI1), var.pMI^0.5)
pMI1 <- plogis(logit.pMI1)

logit.pMI2 <- rnorm(n.years, qlogis(mean.pMI2), var.pMI^0.5)
pMI2 <- plogis(logit.pMI2)

logit.pMI3 <- rnorm(n.years, qlogis(mean.pMI3), var.pMI^0.5)
pMI3 <- plogis(logit.pMI3)

logit.pMI4 <- rnorm(n.years, qlogis(mean.pMI4), var.pMI^0.5)
pMI4 <- plogis(logit.pMI4)

pMI <- c(rbind(pMI1, pMI2, pMI3, pMI4))
pMI <- pMI[1:39]

psiDB.DB <- 0.3
psiDB.JB <- 0.3
psiDB.MI <- 0.4

psiJB.DB <- 0.5
psiJB.JB <- 0.3
psiJB.MI <- 0.2

psiMI.DB <- 0.1
psiMI.JB <- 0.4
psiMI.MI <- 0.5

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(20, n.occasions)
marked[,2] <- rep(20, n.occasions)
marked[,3] <- rep(20, n.occasions)
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
         phiDB[t]*psiDB.DB, phiDB[t]*psiDB.JB, phiDB[t]*psiDB.MI, 1-phiDB[t],
         phiJB[t]*psiJB.DB, phiJB[t]*psiJB.JB, phiJB[t]*psiJB.MI, 1-phiJB[t],
         phiMI[t]*psiMI.DB, phiMI[t]*psiMI.JB, phiMI[t]*psiMI.MI, 1-phiMI[t],
         0,              0,              0,              1),
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

marr <- marray(CH)

# Calculate the number of states
ns <- length(unique(as.numeric(CH))) - 1

# JAGS model code
sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   for (t in 1:(n.occasions-1)){
      for (u in 1:3){
         logit(phi[u,t]) <- mu.phi[u, season[t]] + eps.phi[u,t]
         eps.phi[u,t] ~ dnorm(0, tau.phi[u, season[t]])T(-10,10)
         
         logit(p[u,t]) <- mu.p[u, season[t]] + eps.p[u,t]
         eps.p[u,t] ~ dnorm(0, tau.p[u, season[t]])T(-10,10)
      }
   }

   for (u in 1:3){
      for (s in 1:4){
         mean.phi[u,s] ~ dunif(0, 1)
         mu.phi[u,s] <- logit(mean.phi[u,s])
         sigma.phi[u,s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phi[u,s] <- pow(sigma.phi[u,s], -2)
         sigma2.phi[u,s] <- pow(sigma.phi[u,s], 2)
         
         mean.p[u,s] ~ dunif(0, 1)
         mu.p[u,s] <- logit(mean.p[u,s])
         sigma.p[u,s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.p[u,s] <- pow(sigma.p[u,s], -2)
         sigma2.p[u,s] <- pow(sigma.p[u,s], 2)
      }
   }

# Transitions: gamma priors
   for (i in 1:3){
      db[i] ~ dgamma(1, 1)
      psiDB[i] <- db[i]/sum(db[])
      
      jb[i] ~ dgamma(1, 1)
      psiJB[i] <- jb[i]/sum(jb[])
      
      mi[i] ~ dgamma(1, 1)
      psiMI[i] <- mi[i]/sum(mi[])
   }

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- phi[1,t] * psiDB[1]
      psi[1,t,2] <- phi[1,t] * psiDB[2]
      psi[1,t,3] <- phi[1,t] * psiDB[3]

      
      psi[2,t,1] <- phi[2,t] * psiJB[1]
      psi[2,t,2] <- phi[2,t] * psiJB[2]
      psi[2,t,3] <- phi[2,t] * psiJB[3]


      psi[3,t,1] <- phi[3,t] * psiMI[1]
      psi[3,t,2] <- phi[3,t] * psiMI[2]
      psi[3,t,3] <- phi[3,t] * psiMI[3]
      
      po[1,t] <- p[1,t]
      po[2,t] <- p[2,t]
      po[3,t] <- p[3,t]


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
                  season = rep(c(1,2,3,4), length.out = n.occasions-1))

# Initial values 
inits <- function(){list(mean.phi = matrix(runif(1, 0, 1), nrow = 3, ncol = 4),
                         mean.p = matrix(runif(1, 0, 1), nrow = 3, ncol = 4),
                         sigma.phi = matrix(runif(1, 0, 10), nrow = 3, ncol = 4),
                         sigma.p = matrix(runif(1, 0, 10), nrow = 3, ncol = 4))}  


# Parameters monitored
parameters <- c("mean.phi", "mean.p", "psiDB", "psiJB", "psiMI", "sigma2.phi", "sigma2.p",
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

saveRDS(sim.marr$summary, file = paste0("./analysis-output/seasonal-scenario-1-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/seasonal-scenario-1-simslist", Sys.Date(), ".rds"))

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
   
   png(filename = paste0("analysis-output/parameter-identifiability/seasonal-scenario-1/",
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

png(filename = "figures/seasonal-scenario-1.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phiPreB.sim <- tibble(site = c("DB-PreB", "JB-PreB", "MI-PreB"),
                      value = c(mean.phiDB1, mean.phiJB1, mean.phiMI1))

phiPostB1.sim <- tibble(site = c("DB-PostB1", "JB-PostB1", "MI-PostB1"),
                        value = c(mean.phiDB2, mean.phiJB2, mean.phiMI2))

phiPostB2.sim <- tibble(site = c("DB-PostB2", "JB-PostB2", "MI-PostB2"),
                        value = c(mean.phiDB3, mean.phiJB3, mean.phiMI3))

phiNonB.sim <- tibble(site = c("DB-NonB", "JB-NonB", "MI-NonB"),
                      value = c(mean.phiDB4, mean.phiJB4, mean.phiMI4))

# resighting
pPreB.sim <- tibble(site = c("DB-PreB", "JB-PreB", "MI-PreB"),
                    value = c(mean.pDB1, mean.pJB1, mean.pMI1))

pPostB1.sim <- tibble(site = c("DB-PostB1", "JB-PostB1", "MI-PostB1"),
                      value = c(mean.pDB2, mean.pJB2, mean.pMI2))

pPostB2.sim <- tibble(site = c("DB-PostB2", "JB-PostB2", "MI-PostB2"),
                      value = c(mean.pDB3, mean.pJB3, mean.pMI3))

pNonB.sim <- tibble(site = c("DB-NonB", "JB-NonB", "MI-NonB"),
                    value = c(mean.pDB4, mean.pJB4, mean.pMI4))

# transition
psiDB.sim <- tibble(site = c("DB-DB", "DB-JB", "DB-MI"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI))

psiJB.sim <- tibble(site = c("JB-DB", "JB-JB", "JB-MI"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI))

psiMI.sim <- tibble(site = c("MI-DB", "MI-JB", "MI-MI"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI))

# format survival
phiPreB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.1, mean.phi.2, mean.phi.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.1", "DB-PreB")) %>%
   mutate(site = str_replace(site, "mean.phi.2", "JB-PreB")) %>%
   mutate(site = str_replace(site, "mean.phi.3", "MI-PreB"))

phiPostB1.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.4, mean.phi.5, mean.phi.6) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.4", "DB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.phi.5", "JB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.phi.6", "MI-PostB1"))

phiPostB2.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.7, mean.phi.8, mean.phi.9) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.7", "DB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.phi.8", "JB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.phi.9", "MI-PostB2"))

phiNonB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.10, mean.phi.11, mean.phi.12) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.10", "DB-NonB")) %>%
   mutate(site = str_replace(site, "mean.phi.11", "JB-NonB")) %>%
   mutate(site = str_replace(site, "mean.phi.12", "MI-NonB"))

# plot survival
phiPreB.plot <- ggplot() +
   geom_violin(phiPreB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPreB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiPostB1.plot <- ggplot() +
   geom_violin(phiPostB1.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPostB1.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiPostB2.plot <- ggplot() +
   geom_violin(phiPostB2.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPostB2.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiNonB.plot <- ggplot() +
   geom_violin(phiNonB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiNonB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

# format resighting
pPreB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.1, mean.p.2, mean.p.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.1", "DB-PreB")) %>%
   mutate(site = str_replace(site, "mean.p.2", "JB-PreB")) %>%
   mutate(site = str_replace(site, "mean.p.3", "MI-PreB"))

pPostB1.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.4, mean.p.5, mean.p.6) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.4", "DB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.p.5", "JB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.p.6", "MI-PostB1"))

pPostB2.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.7, mean.p.8, mean.p.9) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.7", "DB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.p.8", "JB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.p.9", "MI-PostB2"))

pNonB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.10, mean.p.11, mean.p.12) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.10", "DB-NonB")) %>%
   mutate(site = str_replace(site, "mean.p.11", "JB-NonB")) %>%
   mutate(site = str_replace(site, "mean.p.12", "MI-NonB"))

# plot resighting
pPreB.plot <- ggplot() +
   geom_violin(pPreB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPreB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pPostB1.plot <- ggplot() +
   geom_violin(pPostB1.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPostB1.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pPostB2.plot <- ggplot() +
   geom_violin(pPostB2.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPostB2.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pNonB.plot <- ggplot() +
   geom_violin(pNonB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pNonB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

# format transitions
psiDB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiDB.1, psiDB.2, psiDB.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiDB.1", "DB-DB")) %>%
   mutate(site = str_replace(site, "psiDB.2", "DB-JB")) %>%
   mutate(site = str_replace(site, "psiDB.3", "DB-MI"))

psiJB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiJB.1, psiJB.2, psiJB.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiJB.1", "JB-DB")) %>%
   mutate(site = str_replace(site, "psiJB.2", "JB-JB")) %>%
   mutate(site = str_replace(site, "psiJB.3", "JB-MI"))

psiMI.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiMI.1, psiMI.2, psiMI.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiMI.1", "MI-DB")) %>%
   mutate(site = str_replace(site, "psiMI.2", "MI-JB")) %>%
   mutate(site = str_replace(site, "psiMI.3", "MI-MI"))

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

# plot all
# survival
png(filename = "figures/seasonal-scenario-1.phiPreB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPreB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiPostB1.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPostB1.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiPostB2.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPostB2.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiNonB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiNonB.plot)

dev.off()

# resighting
png(filename = "figures/seasonal-scenario-1.pPreB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPreB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pPostB1.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPostB1.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pPostB2.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPostB2.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pNonB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pNonB.plot)

dev.off()

# transition
png(filename = "figures/seasonal-scenario-1.psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

#### SCENARIO 2 ####
# States: DB, JB, MI
# phi: seasonal with temporal random effect
# p: seasonal with temporal random effect
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 40
n.years <- 10
n.states <- 8
n.obs <- 8

# DB
mean.phiDB1 <- 0.97
phiDB2 <- rep(0, n.years)
phiDB3 <- rep(0, n.years)
phiDB4 <- rep(0, n.years)

var.phiDB <- 0.3                       # Temporal variance of survival

pDB1 <- rep(0, n.years)
pDB2 <- rep(0, n.years)
pDB3 <- rep(0, n.years)
mean.pDB4 <- 0.7
var.pDB <- 0.5

# Determine annual survival probabilities
logit.phiDB1 <- rnorm(n.years, qlogis(mean.phiDB1), var.phiDB^0.5)
phiDB1 <- plogis(logit.phiDB1)

phiDB <- c(rbind(phiDB1, phiDB2, phiDB3, phiDB4))
phiDB <- phiDB[1:39]

# Determine annual resighting probabilities
logit.pDB4 <- rnorm(n.years, qlogis(mean.pDB4), var.pDB^0.5)
pDB4 <- plogis(logit.pDB4)

pDB <- c(rbind(pDB1, pDB2, pDB3, pDB4))
pDB <- pDB[1:39]

# JB
phiJB1 <- rep(0, n.years)
mean.phiJB2 <- 0.89
phiJB3 <- rep(0, n.years)
phiJB4 <- rep(0, n.years)
var.phiJB <- 0.3                       # Temporal variance of survival

mean.pJB1 <- 0.4
pJB2 <- rep(0, n.years)
pJB3 <- rep(0, n.years)
pJB4 <- rep(0, n.years)
var.pJB <- 0.5

# Determine annual survival probabilities
logit.phiJB2 <- rnorm(n.years, qlogis(mean.phiJB2), var.phiJB^0.5)
phiJB2 <- plogis(logit.phiJB2)

phiJB <- c(rbind(phiJB1, phiJB2, phiJB3, phiJB4))
phiJB <- phiJB[1:39]

# Determine annual resighting probabilities
logit.pJB1 <- rnorm(n.years, qlogis(mean.pJB1), var.pJB^0.5)
pJB1 <- plogis(logit.pJB1)

pJB <- c(rbind(pJB1, pJB2, pJB3, pJB4))
pJB <- pJB[1:39]

# MI
phiMI1 <- rep(0, n.years)
mean.phiMI2 <- 0.89
phiMI3 <- rep(0, n.years)
phiMI4 <- rep(0, n.years)
var.phiMI <- 0.3                       # Temporal variance of survival

mean.pMI1 <- 0.5
pMI2 <- rep(0, n.years)
pMI3 <- rep(0, n.years)
pMI4 <- rep(0, n.years)
var.pMI <- 0.5

# Determine annual survival probabilities
logit.phiMI2 <- rnorm(n.years, qlogis(mean.phiMI2), var.phiMI^0.5)
phiMI2 <- plogis(logit.phiMI2)

phiMI <- c(rbind(phiMI1, phiMI2, phiMI3, phiMI4))
phiMI <- phiMI[1:39]

# Determine annual resighting probabilities
logit.pMI1 <- rnorm(n.years, qlogis(mean.pMI1), var.pMI^0.5)
pMI1 <- plogis(logit.pMI1)

pMI <- c(rbind(pMI1, pMI2, pMI3, pMI4))
pMI <- pMI[1:39]

# CC
phiCC1 <- rep(0, n.years)
mean.phiCC2 <- 0.85
mean.phiCC3 <- 0.85
phiCC4 <- rep(0, n.years)
var.phiCC <- 0.3                       # Temporal variance of survival

mean.pCC1 <- 0.6
mean.pCC2 <- 0.6
pCC3 <- rep(0, n.years)
pCC4 <- rep(0, n.years)
var.pCC <- 0.5

# Determine annual survival probabilities
logit.phiCC2 <- rnorm(n.years, qlogis(mean.phiCC2), var.phiCC^0.5)
phiCC2 <- plogis(logit.phiCC2)

logit.phiCC3 <- rnorm(n.years, qlogis(mean.phiCC3), var.phiCC^0.5)
phiCC3 <- plogis(logit.phiCC3)

phiCC <- c(rbind(phiCC1, phiCC2, phiCC3, phiCC4))
phiCC <- phiCC[1:39]

# Determine annual resighting probabilities
logit.pCC1 <- rnorm(n.years, qlogis(mean.pCC1), var.pCC^0.5)
pCC1 <- plogis(logit.pCC1)

logit.pCC2 <- rnorm(n.years, qlogis(mean.pCC2), var.pCC^0.5)
pCC2 <- plogis(logit.pCC2)

pCC <- c(rbind(pCC1, pCC2, pCC3, pCC4))
pCC <- pCC[1:39]

# SE
mean.phiSE1 <- 0.85
mean.phiSE2 <- 0.9
mean.phiSE3 <- 0.9
mean.phiSE4 <- 0.95
var.phiSE <- 0.3                       # Temporal variance of survival

mean.pSE1 <- 0.6
mean.pSE2 <- 0.6
mean.pSE3 <- 0.4
mean.pSE4 <- 0.3
var.pSE <- 0.5

# Determine annual survival probabilities
logit.phiSE1 <- rnorm(n.years, qlogis(mean.phiSE1), var.phiSE^0.5)
phiSE1 <- plogis(logit.phiSE1)

logit.phiSE2 <- rnorm(n.years, qlogis(mean.phiSE2), var.phiSE^0.5)
phiSE2 <- plogis(logit.phiSE2)

logit.phiSE3 <- rnorm(n.years, qlogis(mean.phiSE3), var.phiSE^0.5)
phiSE3 <- plogis(logit.phiSE3)

logit.phiSE4 <- rnorm(n.years, qlogis(mean.phiSE4), var.phiSE^0.5)
phiSE4 <- plogis(logit.phiSE4)

phiSE <- c(rbind(phiSE1, phiSE2, phiSE3, phiSE4))
phiSE <- phiSE[1:39]

# Determine annual resighting probabilities
logit.pSE1 <- rnorm(n.years, qlogis(mean.pSE1), var.pSE^0.5)
pSE1 <- plogis(logit.pSE1)

logit.pSE2 <- rnorm(n.years, qlogis(mean.pSE2), var.pSE^0.5)
pSE2 <- plogis(logit.pSE2)

logit.pSE3 <- rnorm(n.years, qlogis(mean.pSE3), var.pSE^0.5)
pSE3 <- plogis(logit.pSE3)

logit.pSE4 <- rnorm(n.years, qlogis(mean.pSE4), var.pSE^0.5)
pSE4 <- plogis(logit.pSE4)

pSE <- c(rbind(pSE1, pSE2, pSE3, pSE4))
pSE <- pSE[1:39]

# BR
phiBR1 <- rep(0, n.years)
phiBR2 <- rep(0, n.years)
mean.phiBR3 <- 0.85
mean.phiBR4 <- 0.92
var.phiBR <- 0.3                       # Temporal variance of survival

pBR1 <- rep(0, n.years)
mean.pBR2 <- 0.5
mean.pBR3 <- 0.3
pBR4 <- rep(0, n.years)
var.pBR <- 0.5

# Determine annual survival probabilities
logit.phiBR3 <- rnorm(n.years, qlogis(mean.phiBR3), var.phiBR^0.5)
phiBR3 <- plogis(logit.phiBR3)

logit.phiBR4 <- rnorm(n.years, qlogis(mean.phiBR4), var.phiBR^0.5)
phiBR4 <- plogis(logit.phiBR4)

phiBR <- c(rbind(phiBR1, phiBR2, phiBR3, phiBR4))
phiBR <- phiBR[1:39]

# Determine annual resighting probabilities
logit.pBR2 <- rnorm(n.years, qlogis(mean.pBR2), var.pBR^0.5)
pBR2 <- plogis(logit.pBR2)

logit.pBR3 <- rnorm(n.years, qlogis(mean.pBR3), var.pBR^0.5)
pBR3 <- plogis(logit.pBR3)

pBR <- c(rbind(pBR1, pBR2, pBR3, pBR4))
pBR <- pBR[1:39]

# AR
phiAR1 <- rep(0, n.years)
phiAR2 <- rep(0, n.years)
phiAR3 <- rep(0, n.years)
mean.phiAR4 <- 0.93
var.phiAR <- 0.3                       # Temporal variance of survival

pAR1 <- rep(0, n.years)
pAR2 <- rep(0, n.years)
mean.pAR3 <- 0.6
pAR4 <- rep(0, n.years)
var.pAR <- 0.5

# Determine annual survival probabilities
logit.phiAR4 <- rnorm(n.years, qlogis(mean.phiAR4), var.phiAR^0.5)
phiAR4 <- plogis(logit.phiAR4)

phiAR <- c(rbind(phiAR1, phiAR2, phiAR3, phiAR4))
phiAR <- phiAR[1:39]

# Determine annual resighting probabilities
logit.pAR3 <- rnorm(n.years, qlogis(mean.pAR3), var.pAR^0.5)
pAR3 <- plogis(logit.pAR3)

pAR <- c(rbind(pAR1, pAR2, pAR3, pAR4))
pAR <- pAR[1:39]

# transition probabilities
psiDB.DB <- rep(0, n.occasions)
psiDB.JB <- rep(c(0.35, 0, 0, 0), length.out = n.occasions)
psiDB.MI <- rep(c(0.30, 0, 0, 0), length.out = n.occasions)
psiDB.CC <- rep(c(0.15, 0, 0, 0), length.out = n.occasions)
psiDB.SE <- rep(c(0.10, 0, 0, 0), length.out = n.occasions)
psiDB.BR <- rep(0, n.occasions)
psiDB.AR <- rep(0, n.occasions)

psiJB.DB <- rep(0, n.occasions)
psiJB.JB <- rep(0, n.occasions)
psiJB.MI <- rep(0, n.occasions)
psiJB.CC <- rep(c(0, 0.45, 0, 0), length.out = n.occasions)
psiJB.SE <- rep(c(0, 0.35, 0, 0), length.out = n.occasions)
psiJB.BR <- rep(c(0, 0.20, 0, 0), length.out = n.occasions)
psiJB.AR <- rep(0, n.occasions)

psiMI.DB <- rep(0, n.occasions)
psiMI.JB <- rep(0, n.occasions)
psiMI.MI <- rep(0, n.occasions)
psiMI.CC <- rep(c(0, 0.35, 0, 0), length.out = n.occasions)
psiMI.SE <- rep(c(0, 0.45, 0, 0), length.out = n.occasions)
psiMI.BR <- rep(c(0, 0.20, 0, 0), length.out = n.occasions)
psiMI.AR <- rep(0, n.occasions)

psiCC.DB <- rep(0, n.occasions)
psiCC.JB <- rep(0, n.occasions)
psiCC.MI <- rep(0, n.occasions)
psiCC.CC <- rep(c(0, 0.20, 0, 0), length.out = n.occasions)
psiCC.SE <- rep(c(0, 0.20, 0.20, 0), length.out = n.occasions)
psiCC.BR <- rep(c(0, 0.35, 0.35, 0), length.out = n.occasions)
psiCC.AR <- rep(c(0, 0, 0.25, 0), length.out = n.occasions)

psiSE.DB <- rep(c(0, 0, 0, 0.20), length.out = n.occasions)
psiSE.JB <- rep(c(0.05, 0, 0, 0), length.out = n.occasions)
psiSE.MI <- rep(c(0.05, 0, 0, 0), length.out = n.occasions)
psiSE.CC <- rep(c(0.20, 0.20, 0, 0), length.out = n.occasions)
psiSE.SE <- rep(0.20, n.occasions)
psiSE.BR <- rep(c(0, 0, 0.25, 0.25), length.out = n.occasions)
psiSE.AR <- rep(c(0, 0, 0, 0.05), length.out = n.occasions)

psiBR.DB <- rep(c(0, 0, 0, 0.3), length.out = n.occasions)
psiBR.JB <- rep(0, n.occasions)
psiBR.MI <- rep(0, n.occasions)
psiBR.CC <- rep(0, n.occasions)
psiBR.SE <- rep(c(0, 0, 0.5, 0.5), length.out = n.occasions)
psiBR.BR <- rep(c(0, 0, 0.1, 0), length.out = n.occasions)
psiBR.AR <- rep(c(0, 0, 0.1, 0), length.out = n.occasions)

psiAR.DB <- rep(c(0, 0, 0, 0.90), length.out = n.occasions)
psiAR.JB <- rep(0, n.occasions)
psiAR.MI <- rep(0, n.occasions)
psiAR.CC <- rep(0, n.occasions)
psiAR.SE <- rep(c(0, 0, 0, 0.10), length.out = n.occasions)
psiAR.BR <- rep(0, n.occasions)
psiAR.AR <- rep(0, n.occasions)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(100, 0, 0, 0), length.out = n.occasions) # DB
marked[,2] <- rep(c(0, 40, 0, 0), length.out = n.occasions) # JB
marked[,3] <- rep(c(0, 30, 0, 0), length.out = n.occasions) # MI
marked[,4] <- rep(c(0, 30, 30, 0), length.out = n.occasions) # CC
marked[,5] <- rep(c(5, 20, 10, 5), length.out = n.occasions) # SE
marked[,6] <- rep(c(0, 0, 10, 15), length.out = n.occasions) # BR
marked[,7] <- rep(c(0, 0, 0, 40), length.out = n.occasions) # AR
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
         phiDB[t]*psiDB.DB[t], phiDB[t]*psiDB.JB[t], phiDB[t]*psiDB.MI[t], phiDB[t]*psiDB.CC[t], phiDB[t]*psiDB.SE[t], phiDB[t]*psiDB.BR[t], phiDB[t]*psiDB.AR[t], 1-phiDB[t],
         phiJB[t]*psiJB.DB[t], phiJB[t]*psiJB.JB[t], phiJB[t]*psiJB.MI[t], phiJB[t]*psiJB.CC[t], phiJB[t]*psiJB.SE[t], phiJB[t]*psiJB.BR[t], phiJB[t]*psiJB.AR[t], 1-phiJB[t],
         phiMI[t]*psiMI.DB[t], phiMI[t]*psiMI.JB[t], phiMI[t]*psiMI.MI[t], phiMI[t]*psiMI.CC[t], phiMI[t]*psiMI.SE[t], phiMI[t]*psiMI.BR[t], phiMI[t]*psiMI.AR[t], 1-phiMI[t],
         phiCC[t]*psiCC.DB[t], phiCC[t]*psiCC.JB[t], phiCC[t]*psiCC.MI[t], phiCC[t]*psiCC.CC[t], phiCC[t]*psiCC.SE[t], phiCC[t]*psiCC.BR[t], phiCC[t]*psiCC.AR[t], 1-phiCC[t],
         phiSE[t]*psiSE.DB[t], phiSE[t]*psiSE.JB[t], phiSE[t]*psiSE.MI[t], phiSE[t]*psiSE.CC[t], phiSE[t]*psiSE.SE[t], phiSE[t]*psiSE.BR[t], phiSE[t]*psiSE.AR[t], 1-phiSE[t],
         phiBR[t]*psiBR.DB[t], phiBR[t]*psiBR.JB[t], phiBR[t]*psiBR.MI[t], phiBR[t]*psiBR.CC[t], phiBR[t]*psiBR.SE[t], phiBR[t]*psiBR.BR[t], phiBR[t]*psiBR.AR[t], 1-phiBR[t],
         phiAR[t]*psiAR.DB[t], phiAR[t]*psiAR.JB[t], phiAR[t]*psiAR.MI[t], phiAR[t]*psiAR.CC[t], phiAR[t]*psiAR.SE[t], phiAR[t]*psiAR.BR[t], phiAR[t]*psiAR.AR[t], 1-phiAR[t],
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

# JAGS model code
sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   
   # DB
   for (t in 1:(n.occasions-1)){
         logit(phiDB[t]) <- mu.phiDB[season[t]] + eps.phiDB[t]
         eps.phiDB[t] ~ dnorm(0, tau.phiDB[season[t]])T(-10,10)
         
         logit(pDB[t]) <- mu.pDB[season[t]] + eps.pDB[t]
         eps.pDB[t] ~ dnorm(0, tau.pDB[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiDB[s] ~ dunif(0, 1)
         mu.phiDB[s] <- logit(mean.phiDB[s])
         sigma.phiDB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiDB[s] <- pow(sigma.phiDB[s], -2)
         sigma2.phiDB[s] <- pow(sigma.phiDB[s], 2)
         
         mean.pDB[s] ~ dunif(0, 1)
         mu.pDB[s] <- logit(mean.pDB[s])
         sigma.pDB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pDB[s] <- pow(sigma.pDB[s], -2)
         sigma2.pDB[s] <- pow(sigma.pDB[s], 2)
      }
   
   # JB
   for (t in 1:(n.occasions-1)){
         logit(phiJB[t]) <- mu.phiJB[season[t]] + eps.phiJB[t]
         eps.phiJB[t] ~ dnorm(0, tau.phiJB[season[t]])T(-10,10)
         
         logit(pJB[t]) <- mu.pJB[season[t]] + eps.pJB[t]
         eps.pJB[t] ~ dnorm(0, tau.pJB[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiJB[s] ~ dunif(0, 1)
         mu.phiJB[s] <- logit(mean.phiJB[s])
         sigma.phiJB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiJB[s] <- pow(sigma.phiJB[s], -2)
         sigma2.phiJB[s] <- pow(sigma.phiJB[s], 2)
         
         mean.pJB[s] ~ dunif(0, 1)
         mu.pJB[s] <- logit(mean.pJB[s])
         sigma.pJB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pJB[s] <- pow(sigma.pJB[s], -2)
         sigma2.pJB[s] <- pow(sigma.pJB[s], 2)
      }   
   
   # MI
   for (t in 1:(n.occasions-1)){
         logit(phiMI[t]) <- mu.phiMI[season[t]] + eps.phiMI[t]
         eps.phiMI[t] ~ dnorm(0, tau.phiMI[season[t]])T(-10,10)
         
         logit(pMI[t]) <- mu.pMI[season[t]] + eps.pMI[t]
         eps.pMI[t] ~ dnorm(0, tau.pMI[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiMI[s] ~ dunif(0, 1)
         mu.phiMI[s] <- logit(mean.phiMI[s])
         sigma.phiMI[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiMI[s] <- pow(sigma.phiMI[s], -2)
         sigma2.phiMI[s] <- pow(sigma.phiMI[s], 2)
         
         mean.pMI[s] ~ dunif(0, 1)
         mu.pMI[s] <- logit(mean.pMI[s])
         sigma.pMI[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pMI[s] <- pow(sigma.pMI[s], -2)
         sigma2.pMI[s] <- pow(sigma.pMI[s], 2)
      }
   
   # CC
   for (t in 1:(n.occasions-1)){
         logit(phiCC[t]) <- mu.phiCC[season[t]] + eps.phiCC[t]
         eps.phiCC[t] ~ dnorm(0, tau.phiCC[season[t]])T(-10,10)
         
         logit(pCC[t]) <- mu.pCC[season[t]] + eps.pCC[t]
         eps.pCC[t] ~ dnorm(0, tau.pCC[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiCC[s] ~ dunif(0, 1)
         mu.phiCC[s] <- logit(mean.phiCC[s])
         sigma.phiCC[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiCC[s] <- pow(sigma.phiCC[s], -2)
         sigma2.phiCC[s] <- pow(sigma.phiCC[s], 2)
         
         mean.pCC[s] ~ dunif(0, 1)
         mu.pCC[s] <- logit(mean.pCC[s])
         sigma.pCC[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pCC[s] <- pow(sigma.pCC[s], -2)
         sigma2.pCC[s] <- pow(sigma.pCC[s], 2)
      }
   
   # SE
   for (t in 1:(n.occasions-1)){
         logit(phiSE[t]) <- mu.phiSE[season[t]] + eps.phiSE[t]
         eps.phiSE[t] ~ dnorm(0, tau.phiSE[season[t]])T(-10,10)
         
         logit(pSE[t]) <- mu.pSE[season[t]] + eps.pSE[t]
         eps.pSE[t] ~ dnorm(0, tau.pSE[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiSE[s] ~ dunif(0, 1)
         mu.phiSE[s] <- logit(mean.phiSE[s])
         sigma.phiSE[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiSE[s] <- pow(sigma.phiSE[s], -2)
         sigma2.phiSE[s] <- pow(sigma.phiSE[s], 2)
         
         mean.pSE[s] ~ dunif(0, 1)
         mu.pSE[s] <- logit(mean.pSE[s])
         sigma.pSE[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pSE[s] <- pow(sigma.pSE[s], -2)
         sigma2.pSE[s] <- pow(sigma.pSE[s], 2)
      }
   
   # BR
   for (t in 1:(n.occasions-1)){
         logit(phiBR[t]) <- mu.phiBR[season[t]] + eps.phiBR[t]
         eps.phiBR[t] ~ dnorm(0, tau.phiBR[season[t]])T(-10,10)
         
         logit(pBR[t]) <- mu.pBR[season[t]] + eps.pBR[t]
         eps.pBR[t] ~ dnorm(0, tau.pBR[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiBR[s] ~ dunif(0, 1)
         mu.phiBR[s] <- logit(mean.phiBR[s])
         sigma.phiBR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiBR[s] <- pow(sigma.phiBR[s], -2)
         sigma2.phiBR[s] <- pow(sigma.phiBR[s], 2)
         
         mean.pBR[s] ~ dunif(0, 1)
         mu.pBR[s] <- logit(mean.pBR[s])
         sigma.pBR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pBR[s] <- pow(sigma.pBR[s], -2)
         sigma2.pBR[s] <- pow(sigma.pBR[s], 2)
      }
      
      # AR
      for (t in 1:(n.occasions-1)){
         logit(phiAR[t]) <- mu.phiAR[season[t]] + eps.phiAR[t]
         eps.phiAR[t] ~ dnorm(0, tau.phiAR[season[t]])T(-10,10)
         
         logit(pAR[t]) <- mu.pAR[season[t]] + eps.pAR[t]
         eps.pAR[t] ~ dnorm(0, tau.pAR[season[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiAR[s] ~ dunif(0, 1)
         mu.phiAR[s] <- logit(mean.phiAR[s])
         sigma.phiAR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiAR[s] <- pow(sigma.phiAR[s], -2)
         sigma2.phiAR[s] <- pow(sigma.phiAR[s], 2)
         
         mean.pAR[s] ~ dunif(0, 1)
         mu.pAR[s] <- logit(mean.pAR[s])
         sigma.pAR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pAR[s] <- pow(sigma.pAR[s], -2)
         sigma2.pAR[s] <- pow(sigma.pAR[s], 2)
      }

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
                  ns = ns, zero = matrix(0, ncol = ns, nrow = ns), ones = diag(ns),
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
                  transAR = rep(c(0, 0, 0, 1), length.out = n.occasions-1),
                  season = rep(c(1, 2, 3, 4), length.out = n.occasions-1))

# seasonDBphi = rep(c(1,0,0,0), length.out = n.occasions-1),
# seasonJBphi = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonMIphi = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonCCphi = rep(c(0,2,3,0), length.out = n.occasions-1),
# seasonSEphi = rep(c(1,2,3,4), length.out = n.occasions-1),
# seasonBRphi = rep(c(0,0,3,4), length.out = n.occasions-1),
# seasonARphi = rep(c(0,0,0,4), length.out = n.occasions-1),
# seasonDBp = rep(c(1,0,0,0), length.out = n.occasions-1),
# seasonJBp = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonMIp = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonCCp = rep(c(0,2,3,0), length.out = n.occasions-1),
# seasonSEp = rep(c(1,2,3,4), length.out = n.occasions-1),
# seasonBRp = rep(c(0,0,3,4), length.out = n.occasions-1),
# seasonARp = rep(c(0,0,0,4), length.out = n.occasions-1)

# Initial values 
inits <- function(){list()}  


# Parameters monitored
parameters <- c("mean.phiDB", "mean.phiJB", "mean.phiMI", "mean.phiCC", "mean.phiSE", "mean.phiBR", "mean.phiAR",
                "mean.pDB", "mean.pJB", "mean.pMI", "mean.pCC", "mean.pSE", "mean.pBR", "mean.pAR",
                "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
                "fit", "fit.new")

# MCMC settings
ni <- 100
nt <- 1
nb <- 70
nc <- 3

# Call JAGS from R (BRT 56 min)
sim.marr <- jags(jags.data, inits, parameters, "simulation_scenarios.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(sim.marr, digits = 3)

saveRDS(sim.marr$summary, file = paste0("./analysis-output/seasonal-scenario-2-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/seasonal-scenario-2-simslist", Sys.Date(), ".rds"))

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
   
   png(filename = paste0("analysis-output/parameter-identifiability/seasonal-scenario-2/",
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

png(filename = "figures/seasonal-scenario-2.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phiPreB.sim <- tibble(site = c("DB-PreB", "JB-PreB", "MI-PreB"),
                      value = c(mean.phiDB1, mean.phiJB1, mean.phiMI1))

phiPostB1.sim <- tibble(site = c("DB-PostB1", "JB-PostB1", "MI-PostB1"),
                        value = c(mean.phiDB2, mean.phiJB2, mean.phiMI2))

phiPostB2.sim <- tibble(site = c("DB-PostB2", "JB-PostB2", "MI-PostB2"),
                        value = c(mean.phiDB3, mean.phiJB3, mean.phiMI3))

phiNonB.sim <- tibble(site = c("DB-NonB", "JB-NonB", "MI-NonB"),
                      value = c(mean.phiDB4, mean.phiJB4, mean.phiMI4))

# resighting
pPreB.sim <- tibble(site = c("DB-PreB", "JB-PreB", "MI-PreB"),
                    value = c(mean.pDB1, mean.pJB1, mean.pMI1))

pPostB1.sim <- tibble(site = c("DB-PostB1", "JB-PostB1", "MI-PostB1"),
                      value = c(mean.pDB2, mean.pJB2, mean.pMI2))

pPostB2.sim <- tibble(site = c("DB-PostB2", "JB-PostB2", "MI-PostB2"),
                      value = c(mean.pDB3, mean.pJB3, mean.pMI3))

pNonB.sim <- tibble(site = c("DB-NonB", "JB-NonB", "MI-NonB"),
                    value = c(mean.pDB4, mean.pJB4, mean.pMI4))

# transition
psiDB.sim <- tibble(site = c("DB-DB", "DB-JB", "DB-MI"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI))

psiJB.sim <- tibble(site = c("JB-DB", "JB-JB", "JB-MI"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI))

psiMI.sim <- tibble(site = c("MI-DB", "MI-JB", "MI-MI"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI))

# format survival
phiPreB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.1, mean.phi.2, mean.phi.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.1", "DB-PreB")) %>%
   mutate(site = str_replace(site, "mean.phi.2", "JB-PreB")) %>%
   mutate(site = str_replace(site, "mean.phi.3", "MI-PreB"))

phiPostB1.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.4, mean.phi.5, mean.phi.6) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.4", "DB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.phi.5", "JB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.phi.6", "MI-PostB1"))

phiPostB2.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.7, mean.phi.8, mean.phi.9) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.7", "DB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.phi.8", "JB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.phi.9", "MI-PostB2"))

phiNonB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.10, mean.phi.11, mean.phi.12) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.10", "DB-NonB")) %>%
   mutate(site = str_replace(site, "mean.phi.11", "JB-NonB")) %>%
   mutate(site = str_replace(site, "mean.phi.12", "MI-NonB"))

# plot survival
phiPreB.plot <- ggplot() +
   geom_violin(phiPreB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPreB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiPostB1.plot <- ggplot() +
   geom_violin(phiPostB1.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPostB1.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiPostB2.plot <- ggplot() +
   geom_violin(phiPostB2.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPostB2.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiNonB.plot <- ggplot() +
   geom_violin(phiNonB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiNonB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

# format resighting
pPreB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.1, mean.p.2, mean.p.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.1", "DB-PreB")) %>%
   mutate(site = str_replace(site, "mean.p.2", "JB-PreB")) %>%
   mutate(site = str_replace(site, "mean.p.3", "MI-PreB"))

pPostB1.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.4, mean.p.5, mean.p.6) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.4", "DB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.p.5", "JB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.p.6", "MI-PostB1"))

pPostB2.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.7, mean.p.8, mean.p.9) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.7", "DB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.p.8", "JB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.p.9", "MI-PostB2"))

pNonB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.10, mean.p.11, mean.p.12) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.10", "DB-NonB")) %>%
   mutate(site = str_replace(site, "mean.p.11", "JB-NonB")) %>%
   mutate(site = str_replace(site, "mean.p.12", "MI-NonB"))

# plot resighting
pPreB.plot <- ggplot() +
   geom_violin(pPreB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPreB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pPostB1.plot <- ggplot() +
   geom_violin(pPostB1.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPostB1.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pPostB2.plot <- ggplot() +
   geom_violin(pPostB2.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPostB2.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pNonB.plot <- ggplot() +
   geom_violin(pNonB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pNonB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

# format transitions
psiDB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiDB.1, psiDB.2, psiDB.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiDB.1", "DB-DB")) %>%
   mutate(site = str_replace(site, "psiDB.2", "DB-JB")) %>%
   mutate(site = str_replace(site, "psiDB.3", "DB-MI"))

psiJB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiJB.1, psiJB.2, psiJB.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiJB.1", "JB-DB")) %>%
   mutate(site = str_replace(site, "psiJB.2", "JB-JB")) %>%
   mutate(site = str_replace(site, "psiJB.3", "JB-MI"))

psiMI.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiMI.1, psiMI.2, psiMI.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiMI.1", "MI-DB")) %>%
   mutate(site = str_replace(site, "psiMI.2", "MI-JB")) %>%
   mutate(site = str_replace(site, "psiMI.3", "MI-MI"))

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

# plot all
# survival
png(filename = "figures/seasonal-scenario-1.phiPreB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPreB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiPostB1.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPostB1.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiPostB2.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPostB2.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiNonB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiNonB.plot)

dev.off()

# resighting
png(filename = "figures/seasonal-scenario-1.pPreB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPreB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pPostB1.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPostB1.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pPostB2.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPostB2.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pNonB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pNonB.plot)

dev.off()

# transition
png(filename = "figures/seasonal-scenario-1.psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()

#### SCENARIO 3 ####
# States: DB, JB, MI
# phi: seasonal with temporal random effect
# p: seasonal with temporal random effect
# psi: all possible

# 1. Model description
# 2. Generation of simulated data

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
n.occasions <- 40
n.years <- 10
n.states <- 8
n.obs <- 8

# DB
mean.phiDB1 <- 0.97
phiDB2 <- rep(0, n.years)
phiDB3 <- rep(0, n.years)
phiDB4 <- rep(0, n.years)

var.phiDB <- 0.3                       # Temporal variance of survival

pDB1 <- rep(0, n.years)
pDB2 <- rep(0, n.years)
pDB3 <- rep(0, n.years)
mean.pDB4 <- 0.7
var.pDB <- 0.5

# Determine annual survival probabilities
logit.phiDB1 <- rnorm(n.years, qlogis(mean.phiDB1), var.phiDB^0.5)
phiDB1 <- plogis(logit.phiDB1)

phiDB <- c(rbind(phiDB1, phiDB2, phiDB3, phiDB4))
phiDB <- phiDB[1:39]

# Determine annual resighting probabilities
logit.pDB4 <- rnorm(n.years, qlogis(mean.pDB4), var.pDB^0.5)
pDB4 <- plogis(logit.pDB4)

pDB <- c(rbind(pDB1, pDB2, pDB3, pDB4))
pDB <- pDB[1:39]

# JB
phiJB1 <- rep(0, n.years)
mean.phiJB2 <- 0.89
phiJB3 <- rep(0, n.years)
phiJB4 <- rep(0, n.years)
var.phiJB <- 0.3                       # Temporal variance of survival

mean.pJB1 <- 0.4
pJB2 <- rep(0, n.years)
pJB3 <- rep(0, n.years)
pJB4 <- rep(0, n.years)
var.pJB <- 0.5

# Determine annual survival probabilities
logit.phiJB2 <- rnorm(n.years, qlogis(mean.phiJB2), var.phiJB^0.5)
phiJB2 <- plogis(logit.phiJB2)

phiJB <- c(rbind(phiJB1, phiJB2, phiJB3, phiJB4))
phiJB <- phiJB[1:39]

# Determine annual resighting probabilities
logit.pJB1 <- rnorm(n.years, qlogis(mean.pJB1), var.pJB^0.5)
pJB1 <- plogis(logit.pJB1)

pJB <- c(rbind(pJB1, pJB2, pJB3, pJB4))
pJB <- pJB[1:39]

# MI
phiMI1 <- rep(0, n.years)
mean.phiMI2 <- 0.89
phiMI3 <- rep(0, n.years)
phiMI4 <- rep(0, n.years)
var.phiMI <- 0.3                       # Temporal variance of survival

mean.pMI1 <- 0.5
pMI2 <- rep(0, n.years)
pMI3 <- rep(0, n.years)
pMI4 <- rep(0, n.years)
var.pMI <- 0.5

# Determine annual survival probabilities
logit.phiMI2 <- rnorm(n.years, qlogis(mean.phiMI2), var.phiMI^0.5)
phiMI2 <- plogis(logit.phiMI2)

phiMI <- c(rbind(phiMI1, phiMI2, phiMI3, phiMI4))
phiMI <- phiMI[1:39]

# Determine annual resighting probabilities
logit.pMI1 <- rnorm(n.years, qlogis(mean.pMI1), var.pMI^0.5)
pMI1 <- plogis(logit.pMI1)

pMI <- c(rbind(pMI1, pMI2, pMI3, pMI4))
pMI <- pMI[1:39]

# CC
phiCC1 <- rep(0, n.years)
mean.phiCC2 <- 0.85
mean.phiCC3 <- 0.85
phiCC4 <- rep(0, n.years)
var.phiCC <- 0.3                       # Temporal variance of survival

mean.pCC1 <- 0.6
mean.pCC2 <- 0.6
pCC3 <- rep(0, n.years)
pCC4 <- rep(0, n.years)
var.pCC <- 0.5

# Determine annual survival probabilities
logit.phiCC2 <- rnorm(n.years, qlogis(mean.phiCC2), var.phiCC^0.5)
phiCC2 <- plogis(logit.phiCC2)

logit.phiCC3 <- rnorm(n.years, qlogis(mean.phiCC3), var.phiCC^0.5)
phiCC3 <- plogis(logit.phiCC3)

phiCC <- c(rbind(phiCC1, phiCC2, phiCC3, phiCC4))
phiCC <- phiCC[1:39]

# Determine annual resighting probabilities
logit.pCC1 <- rnorm(n.years, qlogis(mean.pCC1), var.pCC^0.5)
pCC1 <- plogis(logit.pCC1)

logit.pCC2 <- rnorm(n.years, qlogis(mean.pCC2), var.pCC^0.5)
pCC2 <- plogis(logit.pCC2)

pCC <- c(rbind(pCC1, pCC2, pCC3, pCC4))
pCC <- pCC[1:39]

# SE
mean.phiSE1 <- 0.85
mean.phiSE2 <- 0.9
mean.phiSE3 <- 0.9
mean.phiSE4 <- 0.95
var.phiSE <- 0.3                       # Temporal variance of survival

mean.pSE1 <- 0.6
mean.pSE2 <- 0.6
mean.pSE3 <- 0.4
mean.pSE4 <- 0.3
var.pSE <- 0.5

# Determine annual survival probabilities
logit.phiSE1 <- rnorm(n.years, qlogis(mean.phiSE1), var.phiSE^0.5)
phiSE1 <- plogis(logit.phiSE1)

logit.phiSE2 <- rnorm(n.years, qlogis(mean.phiSE2), var.phiSE^0.5)
phiSE2 <- plogis(logit.phiSE2)

logit.phiSE3 <- rnorm(n.years, qlogis(mean.phiSE3), var.phiSE^0.5)
phiSE3 <- plogis(logit.phiSE3)

logit.phiSE4 <- rnorm(n.years, qlogis(mean.phiSE4), var.phiSE^0.5)
phiSE4 <- plogis(logit.phiSE4)

phiSE <- c(rbind(phiSE1, phiSE2, phiSE3, phiSE4))
phiSE <- phiSE[1:39]

# Determine annual resighting probabilities
logit.pSE1 <- rnorm(n.years, qlogis(mean.pSE1), var.pSE^0.5)
pSE1 <- plogis(logit.pSE1)

logit.pSE2 <- rnorm(n.years, qlogis(mean.pSE2), var.pSE^0.5)
pSE2 <- plogis(logit.pSE2)

logit.pSE3 <- rnorm(n.years, qlogis(mean.pSE3), var.pSE^0.5)
pSE3 <- plogis(logit.pSE3)

logit.pSE4 <- rnorm(n.years, qlogis(mean.pSE4), var.pSE^0.5)
pSE4 <- plogis(logit.pSE4)

pSE <- c(rbind(pSE1, pSE2, pSE3, pSE4))
pSE <- pSE[1:39]

# BR
phiBR1 <- rep(0, n.years)
phiBR2 <- rep(0, n.years)
mean.phiBR3 <- 0.85
mean.phiBR4 <- 0.92
var.phiBR <- 0.3                       # Temporal variance of survival

pBR1 <- rep(0, n.years)
mean.pBR2 <- 0.5
mean.pBR3 <- 0.3
pBR4 <- rep(0, n.years)
var.pBR <- 0.5

# Determine annual survival probabilities
logit.phiBR3 <- rnorm(n.years, qlogis(mean.phiBR3), var.phiBR^0.5)
phiBR3 <- plogis(logit.phiBR3)

logit.phiBR4 <- rnorm(n.years, qlogis(mean.phiBR4), var.phiBR^0.5)
phiBR4 <- plogis(logit.phiBR4)

phiBR <- c(rbind(phiBR1, phiBR2, phiBR3, phiBR4))
phiBR <- phiBR[1:39]

# Determine annual resighting probabilities
logit.pBR2 <- rnorm(n.years, qlogis(mean.pBR2), var.pBR^0.5)
pBR2 <- plogis(logit.pBR2)

logit.pBR3 <- rnorm(n.years, qlogis(mean.pBR3), var.pBR^0.5)
pBR3 <- plogis(logit.pBR3)

pBR <- c(rbind(pBR1, pBR2, pBR3, pBR4))
pBR <- pBR[1:39]

# AR
phiAR1 <- rep(0, n.years)
phiAR2 <- rep(0, n.years)
phiAR3 <- rep(0, n.years)
mean.phiAR4 <- 0.93
var.phiAR <- 0.3                       # Temporal variance of survival

pAR1 <- rep(0, n.years)
pAR2 <- rep(0, n.years)
mean.pAR3 <- 0.6
pAR4 <- rep(0, n.years)
var.pAR <- 0.5

# Determine annual survival probabilities
logit.phiAR4 <- rnorm(n.years, qlogis(mean.phiAR4), var.phiAR^0.5)
phiAR4 <- plogis(logit.phiAR4)

phiAR <- c(rbind(phiAR1, phiAR2, phiAR3, phiAR4))
phiAR <- phiAR[1:39]

# Determine annual resighting probabilities
logit.pAR3 <- rnorm(n.years, qlogis(mean.pAR3), var.pAR^0.5)
pAR3 <- plogis(logit.pAR3)

pAR <- c(rbind(pAR1, pAR2, pAR3, pAR4))
pAR <- pAR[1:39]

# transition probabilities
psiDB.DB <- rep(0, n.occasions)
psiDB.JB <- rep(c(0.35, 0, 0, 0), length.out = n.occasions)
psiDB.MI <- rep(c(0.30, 0, 0, 0), length.out = n.occasions)
psiDB.CC <- rep(c(0.15, 0, 0, 0), length.out = n.occasions)
psiDB.SE <- rep(c(0.10, 0, 0, 0), length.out = n.occasions)
psiDB.BR <- rep(0, n.occasions)
psiDB.AR <- rep(0, n.occasions)

psiJB.DB <- rep(0, n.occasions)
psiJB.JB <- rep(0, n.occasions)
psiJB.MI <- rep(0, n.occasions)
psiJB.CC <- rep(c(0, 0.45, 0, 0), length.out = n.occasions)
psiJB.SE <- rep(c(0, 0.35, 0, 0), length.out = n.occasions)
psiJB.BR <- rep(c(0, 0.20, 0, 0), length.out = n.occasions)
psiJB.AR <- rep(0, n.occasions)

psiMI.DB <- rep(0, n.occasions)
psiMI.JB <- rep(0, n.occasions)
psiMI.MI <- rep(0, n.occasions)
psiMI.CC <- rep(c(0, 0.35, 0, 0), length.out = n.occasions)
psiMI.SE <- rep(c(0, 0.45, 0, 0), length.out = n.occasions)
psiMI.BR <- rep(c(0, 0.20, 0, 0), length.out = n.occasions)
psiMI.AR <- rep(0, n.occasions)

psiCC.DB <- rep(0, n.occasions)
psiCC.JB <- rep(0, n.occasions)
psiCC.MI <- rep(0, n.occasions)
psiCC.CC <- rep(c(0, 0.20, 0, 0), length.out = n.occasions)
psiCC.SE <- rep(c(0, 0.20, 0.20, 0), length.out = n.occasions)
psiCC.BR <- rep(c(0, 0.35, 0.35, 0), length.out = n.occasions)
psiCC.AR <- rep(c(0, 0, 0.25, 0), length.out = n.occasions)

psiSE.DB <- rep(c(0, 0, 0, 0.20), length.out = n.occasions)
psiSE.JB <- rep(c(0.05, 0, 0, 0), length.out = n.occasions)
psiSE.MI <- rep(c(0.05, 0, 0, 0), length.out = n.occasions)
psiSE.CC <- rep(c(0.20, 0.20, 0, 0), length.out = n.occasions)
psiSE.SE <- rep(0.20, n.occasions)
psiSE.BR <- rep(c(0, 0, 0.25, 0.25), length.out = n.occasions)
psiSE.AR <- rep(c(0, 0, 0, 0.05), length.out = n.occasions)

psiBR.DB <- rep(c(0, 0, 0, 0.3), length.out = n.occasions)
psiBR.JB <- rep(0, n.occasions)
psiBR.MI <- rep(0, n.occasions)
psiBR.CC <- rep(0, n.occasions)
psiBR.SE <- rep(c(0, 0, 0.5, 0.5), length.out = n.occasions)
psiBR.BR <- rep(c(0, 0, 0.1, 0), length.out = n.occasions)
psiBR.AR <- rep(c(0, 0, 0.1, 0), length.out = n.occasions)

psiAR.DB <- rep(c(0, 0, 0, 0.90), length.out = n.occasions)
psiAR.JB <- rep(0, n.occasions)
psiAR.MI <- rep(0, n.occasions)
psiAR.CC <- rep(0, n.occasions)
psiAR.SE <- rep(c(0, 0, 0, 0.10), length.out = n.occasions)
psiAR.BR <- rep(0, n.occasions)
psiAR.AR <- rep(0, n.occasions)

# released individuals
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(c(100, 0, 0, 0), length.out = n.occasions) # DB
marked[,2] <- rep(c(0, 40, 0, 0), length.out = n.occasions) # JB
marked[,3] <- rep(c(0, 30, 0, 0), length.out = n.occasions) # MI
marked[,4] <- rep(c(0, 30, 30, 0), length.out = n.occasions) # CC
marked[,5] <- rep(c(5, 20, 10, 5), length.out = n.occasions) # SE
marked[,6] <- rep(c(0, 0, 10, 15), length.out = n.occasions) # BR
marked[,7] <- rep(c(0, 0, 0, 40), length.out = n.occasions) # AR
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
         phiDB[t]*psiDB.DB[t], phiDB[t]*psiDB.JB[t], phiDB[t]*psiDB.MI[t], phiDB[t]*psiDB.CC[t], phiDB[t]*psiDB.SE[t], phiDB[t]*psiDB.BR[t], phiDB[t]*psiDB.AR[t], 1-phiDB[t],
         phiJB[t]*psiJB.DB[t], phiJB[t]*psiJB.JB[t], phiJB[t]*psiJB.MI[t], phiJB[t]*psiJB.CC[t], phiJB[t]*psiJB.SE[t], phiJB[t]*psiJB.BR[t], phiJB[t]*psiJB.AR[t], 1-phiJB[t],
         phiMI[t]*psiMI.DB[t], phiMI[t]*psiMI.JB[t], phiMI[t]*psiMI.MI[t], phiMI[t]*psiMI.CC[t], phiMI[t]*psiMI.SE[t], phiMI[t]*psiMI.BR[t], phiMI[t]*psiMI.AR[t], 1-phiMI[t],
         phiCC[t]*psiCC.DB[t], phiCC[t]*psiCC.JB[t], phiCC[t]*psiCC.MI[t], phiCC[t]*psiCC.CC[t], phiCC[t]*psiCC.SE[t], phiCC[t]*psiCC.BR[t], phiCC[t]*psiCC.AR[t], 1-phiCC[t],
         phiSE[t]*psiSE.DB[t], phiSE[t]*psiSE.JB[t], phiSE[t]*psiSE.MI[t], phiSE[t]*psiSE.CC[t], phiSE[t]*psiSE.SE[t], phiSE[t]*psiSE.BR[t], phiSE[t]*psiSE.AR[t], 1-phiSE[t],
         phiBR[t]*psiBR.DB[t], phiBR[t]*psiBR.JB[t], phiBR[t]*psiBR.MI[t], phiBR[t]*psiBR.CC[t], phiBR[t]*psiBR.SE[t], phiBR[t]*psiBR.BR[t], phiBR[t]*psiBR.AR[t], 1-phiBR[t],
         phiAR[t]*psiAR.DB[t], phiAR[t]*psiAR.JB[t], phiAR[t]*psiAR.MI[t], phiAR[t]*psiAR.CC[t], phiAR[t]*psiAR.SE[t], phiAR[t]*psiAR.BR[t], phiAR[t]*psiAR.AR[t], 1-phiAR[t],
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

# JAGS model code
sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   
   # DB
   for (t in 1:10){
         logit(phiDB[t]) <- mu.phiDB + eps.phiDB[t]
         eps.phiDB[t] ~ dnorm(0, tau.phiDB)T(-10,10)
   }
   
   for (t in 1:9){
         logit(pDB[t]) <- mu.pDB + eps.pDB[t]
         eps.pDB[t] ~ dnorm(0, tau.pDB)T(-10,10)
   }

         mean.phiDB ~ dunif(0, 1)
         mu.phiDB <- logit(mean.phiDB)
         sigma.phiDB ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiDB <- pow(sigma.phiDB, -2)
         sigma2.phiDB <- pow(sigma.phiDB, 2)
      
         mean.pDB ~ dunif(0, 1)
         mu.pDB <- logit(mean.pDB)
         sigma.pDB ~ dunif(0, 10)               # Prior for standard deviation
         tau.pDB <- pow(sigma.pDB, -2)
         sigma2.pDB <- pow(sigma.pDB, 2)
         
phitDB[1] <- phiDB[1]
phitDB[2] <- 0
phitDB[3] <- 0
phitDB[4] <- 0
phitDB[5] <- phiDB[2]
phitDB[6] <- 0
phitDB[7] <- 0
phitDB[8] <- 0
phitDB[9] <- phiDB[3]
phitDB[10] <- 0
phitDB[11] <- 0
phitDB[12] <- 0
phitDB[13] <- phiDB[4]
phitDB[14] <- 0
phitDB[15] <- 0
phitDB[16] <- 0
phitDB[17] <- phiDB[5]
phitDB[18] <- 0
phitDB[19] <- 0
phitDB[20] <- 0
phitDB[21] <- phiDB[6]
phitDB[22] <- 0
phitDB[23] <- 0
phitDB[24] <- 0
phitDB[25] <- phiDB[7]
phitDB[26] <- 0
phitDB[27] <- 0
phitDB[28] <- 0
phitDB[29] <- phiDB[8]
phitDB[30] <- 0
phitDB[31] <- 0
phitDB[32] <- 0
phitDB[33] <- phiDB[9]
phitDB[34] <- 0
phitDB[35] <- 0
phitDB[36] <- 0
phitDB[37] <- phiDB[10]
phitDB[38] <- 0
phitDB[39] <- 0

psightDB[1] <- 0
psightDB[2] <- 0
psightDB[3] <- 0
psightDB[4] <- pDB[1]
psightDB[5] <- 0
psightDB[6] <- 0
psightDB[7] <- 0
psightDB[8] <- pDB[2]
psightDB[9] <- 0
psightDB[10] <- 0
psightDB[11] <- 0
psightDB[12] <- pDB[3]
psightDB[13] <- 0
psightDB[14] <- 0
psightDB[15] <- 0
psightDB[16] <- pDB[4]
psightDB[17] <- 0
psightDB[18] <- 0
psightDB[19] <- 0
psightDB[20] <- pDB[5]
psightDB[21] <- 0
psightDB[22] <- 0
psightDB[23] <- 0
psightDB[24] <- pDB[6]
psightDB[25] <- 0
psightDB[26] <- 0
psightDB[27] <- 0
psightDB[28] <- pDB[7]
psightDB[29] <- 0
psightDB[30] <- 0
psightDB[31] <- 0
psightDB[32] <- pDB[8]
psightDB[33] <- 0
psightDB[34] <- 0
psightDB[35] <- 0
psightDB[36] <- pDB[9]
psightDB[37] <- 0
psightDB[38] <- 0
psightDB[39] <- 0
   
   # JB
   for (t in 1:10){
         logit(phiJB[t]) <- mu.phiJB + eps.phiJB[t]
         eps.phiJB[t] ~ dnorm(0, tau.phiJB)T(-10,10)
   }
   
   for(t in 1:10){
         logit(pJB[t]) <- mu.pJB + eps.pJB[t]
         eps.pJB[t] ~ dnorm(0, tau.pJB)T(-10,10)
   }

         mean.phiJB ~ dunif(0, 1)
         mu.phiJB <- logit(mean.phiJB)
         sigma.phiJB ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiJB <- pow(sigma.phiJB, -2)
         sigma2.phiJB <- pow(sigma.phiJB, 2)
      
         mean.pJB ~ dunif(0, 1)
         mu.pJB <- logit(mean.pJB)
         sigma.pJB ~ dunif(0, 10)               # Prior for standard deviation
         tau.pJB <- pow(sigma.pJB, -2)
         sigma2.pJB <- pow(sigma.pJB, 2)

phitJB[1] <- 0
phitJB[2] <- phiJB[1]
phitJB[3] <- 0
phitJB[4] <- 0
phitJB[5] <- 0
phitJB[6] <- phiJB[2]
phitJB[7] <- 0
phitJB[8] <- 0
phitJB[9] <- 0
phitJB[10] <- phiJB[3]
phitJB[11] <- 0
phitJB[12] <- 0
phitJB[13] <- 0
phitJB[14] <- phiJB[4]
phitJB[15] <- 0
phitJB[16] <- 0
phitJB[17] <- 0
phitJB[18] <- phiJB[5]
phitJB[19] <- 0
phitJB[20] <- 0
phitJB[21] <- 0
phitJB[22] <- phiJB[6]
phitJB[23] <- 0
phitJB[24] <- 0
phitJB[25] <- 0
phitJB[26] <- phiJB[7]
phitJB[27] <- 0
phitJB[28] <- 0
phitJB[29] <- 0
phitJB[30] <- phiJB[8]
phitJB[31] <- 0
phitJB[32] <- 0
phitJB[33] <- 0
phitJB[34] <- phiJB[9]
phitJB[35] <- 0
phitJB[36] <- 0
phitJB[37] <- 0
phitJB[38] <- phiJB[10]
phitJB[39] <- 0
         
psightJB[1] <- pJB[1]
psightJB[2] <- 0
psightJB[3] <- 0
psightJB[4] <- 0
psightJB[5] <- pJB[2]
psightJB[6] <- 0
psightJB[7] <- 0
psightJB[8] <- 0
psightJB[9] <- pJB[3]
psightJB[10] <- 0
psightJB[11] <- 0
psightJB[12] <- 0
psightJB[13] <- pJB[4]
psightJB[14] <- 0
psightJB[15] <- 0
psightJB[16] <- 0
psightJB[17] <- pJB[5]
psightJB[18] <- 0
psightJB[19] <- 0
psightJB[20] <- 0
psightJB[21] <- pJB[6]
psightJB[22] <- 0
psightJB[23] <- 0
psightJB[24] <- 0
psightJB[25] <- pJB[7]
psightJB[26] <- 0
psightJB[27] <- 0
psightJB[28] <- 0
psightJB[29] <- pJB[8]
psightJB[30] <- 0
psightJB[31] <- 0
psightJB[32] <- 0
psightJB[33] <- pJB[9]
psightJB[34] <- 0
psightJB[35] <- 0
psightJB[36] <- 0
psightJB[37] <- pJB[10]
psightJB[38] <- 0
psightJB[39] <- 0
   
   # MI
   for (t in 1:10){
         logit(phiMI[t]) <- mu.phiMI + eps.phiMI[t]
         eps.phiMI[t] ~ dnorm(0, tau.phiMI)T(-10,10)
   }
   
   for (t in 1:10){
         logit(pMI[t]) <- mu.pMI + eps.pMI[t]
         eps.pMI[t] ~ dnorm(0, tau.pMI)T(-10,10)
   }

         mean.phiMI ~ dunif(0, 1)
         mu.phiMI <- logit(mean.phiMI)
         sigma.phiMI ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiMI <- pow(sigma.phiMI, -2)
         sigma2.phiMI <- pow(sigma.phiMI, 2)

         mean.pMI ~ dunif(0, 1)
         mu.pMI <- logit(mean.pMI)
         sigma.pMI ~ dunif(0, 10)               # Prior for standard deviation
         tau.pMI <- pow(sigma.pMI, -2)
         sigma2.pMI <- pow(sigma.pMI, 2)

phitMI[1] <- 0
phitMI[2] <- phiMI[1]
phitMI[3] <- 0
phitMI[4] <- 0
phitMI[5] <- 0
phitMI[6] <- phiMI[2]
phitMI[7] <- 0
phitMI[8] <- 0
phitMI[9] <- 0
phitMI[10] <- phiMI[3]
phitMI[11] <- 0
phitMI[12] <- 0
phitMI[13] <- 0
phitMI[14] <- phiMI[4]
phitMI[15] <- 0
phitMI[16] <- 0
phitMI[17] <- 0
phitMI[18] <- phiMI[5]
phitMI[19] <- 0
phitMI[20] <- 0
phitMI[21] <- 0
phitMI[22] <- phiMI[6]
phitMI[23] <- 0
phitMI[24] <- 0
phitMI[25] <- 0
phitMI[26] <- phiMI[7]
phitMI[27] <- 0
phitMI[28] <- 0
phitMI[29] <- 0
phitMI[30] <- phiMI[8]
phitMI[31] <- 0
phitMI[32] <- 0
phitMI[33] <- 0
phitMI[34] <- phiMI[9]
phitMI[35] <- 0
phitMI[36] <- 0
phitMI[37] <- 0
phitMI[38] <- phiMI[10]
phitMI[39] <- 0

psightMI[1] <- pMI[1]
psightMI[2] <- 0
psightMI[3] <- 0
psightMI[4] <- 0
psightMI[5] <- pMI[2]
psightMI[6] <- 0
psightMI[7] <- 0
psightMI[8] <- 0
psightMI[9] <- pMI[3]
psightMI[10] <- 0
psightMI[11] <- 0
psightMI[12] <- 0
psightMI[13] <- pMI[4]
psightMI[14] <- 0
psightMI[15] <- 0
psightMI[16] <- 0
psightMI[17] <- pMI[5]
psightMI[18] <- 0
psightMI[19] <- 0
psightMI[20] <- 0
psightMI[21] <- pMI[6]
psightMI[22] <- 0
psightMI[23] <- 0
psightMI[24] <- 0
psightMI[25] <- pMI[7]
psightMI[26] <- 0
psightMI[27] <- 0
psightMI[28] <- 0
psightMI[29] <- pMI[8]
psightMI[30] <- 0
psightMI[31] <- 0
psightMI[32] <- 0
psightMI[33] <- pMI[9]
psightMI[34] <- 0
psightMI[35] <- 0
psightMI[36] <- 0
psightMI[37] <- pMI[10]
psightMI[38] <- 0
psightMI[39] <- 0
   
   # CC
   for (t in 1:20){
         logit(phiCC[t]) <- mu.phiCC[seasonphiCC[t]] + eps.phiCC[t]
         eps.phiCC[t] ~ dnorm(0, tau.phiCC[seasonphiCC[t]])T(-10,10)
   }
   
   for (t in 1:20){
         logit(pCC[t]) <- mu.pCC[seasonpCC[t]] + eps.pCC[t]
         eps.pCC[t] ~ dnorm(0, tau.pCC[seasonpCC[t]])T(-10,10)
   }

      for (s in 1:2){
         mean.phiCC[s] ~ dunif(0, 1)
         mu.phiCC[s] <- logit(mean.phiCC[s])
         sigma.phiCC[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiCC[s] <- pow(sigma.phiCC[s], -2)
         sigma2.phiCC[s] <- pow(sigma.phiCC[s], 2)
      }
      
      for (s in 1:2){
         mean.pCC[s] ~ dunif(0, 1)
         mu.pCC[s] <- logit(mean.pCC[s])
         sigma.pCC[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pCC[s] <- pow(sigma.pCC[s], -2)
         sigma2.pCC[s] <- pow(sigma.pCC[s], 2)
      }

phitCC[1] <- 0
phitCC[2] <- phiCC[1]
phitCC[3] <- phiCC[2]
phitCC[4] <- 0
phitCC[5] <- 0
phitCC[6] <- phiCC[3]
phitCC[7] <- phiCC[4]
phitCC[8] <- 0
phitCC[9] <- 0
phitCC[10] <- phiCC[5]
phitCC[11] <- phiCC[6]
phitCC[12] <- 0
phitCC[13] <- 0
phitCC[14] <- phiCC[7]
phitCC[15] <- phiCC[8]
phitCC[16] <- 0
phitCC[17] <- 0
phitCC[18] <- phiCC[9]
phitCC[19] <- phiCC[10]
phitCC[20] <- 0
phitCC[21] <- 0
phitCC[22] <- phiCC[11]
phitCC[23] <- phiCC[12]
phitCC[24] <- 0
phitCC[25] <- 0
phitCC[26] <- phiCC[13]
phitCC[27] <- phiCC[14]
phitCC[28] <- 0
phitCC[29] <- 0
phitCC[30] <- phiCC[15]
phitCC[31] <- phiCC[16]
phitCC[32] <- 0
phitCC[33] <- 0
phitCC[34] <- phiCC[17]
phitCC[35] <- phiCC[18]
phitCC[36] <- 0
phitCC[37] <- 0
phitCC[38] <- phiCC[19]
phitCC[39] <- phiCC[20]
 
psightCC[1] <- pCC[1]
psightCC[2] <- pCC[2]
psightCC[3] <- 0
psightCC[4] <- 0
psightCC[5] <- pCC[3]
psightCC[6] <- pCC[4]
psightCC[7] <- 0
psightCC[8] <- 0
psightCC[9] <- pCC[5]
psightCC[10] <- pCC[6]
psightCC[11] <- 0
psightCC[12] <- 0
psightCC[13] <- pCC[7]
psightCC[14] <- pCC[8]
psightCC[15] <- 0
psightCC[16] <- 0
psightCC[17] <- pCC[9]
psightCC[18] <- pCC[10]
psightCC[19] <- 0
psightCC[20] <- 0
psightCC[21] <- pCC[11]
psightCC[22] <- pCC[12]
psightCC[23] <- 0
psightCC[24] <- 0
psightCC[25] <- pCC[13]
psightCC[26] <- pCC[14]
psightCC[27] <- 0
psightCC[28] <- 0
psightCC[29] <- pCC[15]
psightCC[30] <- pCC[16]
psightCC[31] <- 0
psightCC[32] <- 0
psightCC[33] <- pCC[17]
psightCC[34] <- pCC[18]
psightCC[35] <- 0
psightCC[36] <- 0
psightCC[37] <- pCC[19]
psightCC[38] <- pCC[20]
psightCC[39] <- 0

   # SE
   for (t in 1:(n.occasions-1)){
         logit(phiSE[t]) <- mu.phiSE[seasonphiSE[t]] + eps.phiSE[t]
         eps.phiSE[t] ~ dnorm(0, tau.phiSE[seasonphiSE[t]])T(-10,10)
         
         logit(pSE[t]) <- mu.pSE[seasonpSE[t]] + eps.pSE[t]
         eps.pSE[t] ~ dnorm(0, tau.pSE[seasonpSE[t]])T(-10,10)
   }

      for (s in 1:4){
         mean.phiSE[s] ~ dunif(0, 1)
         mu.phiSE[s] <- logit(mean.phiSE[s])
         sigma.phiSE[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiSE[s] <- pow(sigma.phiSE[s], -2)
         sigma2.phiSE[s] <- pow(sigma.phiSE[s], 2)
         
         mean.pSE[s] ~ dunif(0, 1)
         mu.pSE[s] <- logit(mean.pSE[s])
         sigma.pSE[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pSE[s] <- pow(sigma.pSE[s], -2)
         sigma2.pSE[s] <- pow(sigma.pSE[s], 2)
      }

phitSE[1] <- phiSE[1]
phitSE[2] <- phiSE[2]
phitSE[3] <- phiSE[3]
phitSE[4] <- phiSE[4]
phitSE[5] <- phiSE[5]
phitSE[6] <- phiSE[6]
phitSE[7] <- phiSE[7]
phitSE[8] <- phiSE[8]
phitSE[9] <- phiSE[9]
phitSE[10] <- phiSE[10]
phitSE[11] <- phiSE[11]
phitSE[12] <- phiSE[12]
phitSE[13] <- phiSE[13]
phitSE[14] <- phiSE[14]
phitSE[15] <- phiSE[15]
phitSE[16] <- phiSE[16]
phitSE[17] <- phiSE[17]
phitSE[18] <- phiSE[18]
phitSE[19] <- phiSE[19]
phitSE[20] <- phiSE[20]
phitSE[21] <- phiSE[21]
phitSE[22] <- phiSE[22]
phitSE[23] <- phiSE[23]
phitSE[24] <- phiSE[24]
phitSE[25] <- phiSE[25]
phitSE[26] <- phiSE[26]
phitSE[27] <- phiSE[27]
phitSE[28] <- phiSE[28]
phitSE[29] <- phiSE[29]
phitSE[30] <- phiSE[30]
phitSE[31] <- phiSE[31]
phitSE[32] <- phiSE[32]
phitSE[33] <- phiSE[33]
phitSE[34] <- phiSE[34]
phitSE[35] <- phiSE[35]
phitSE[36] <- phiSE[36]
phitSE[37] <- phiSE[37]
phitSE[38] <- phiSE[38]
phitSE[39] <- phiSE[39]

psightSE[1] <- pSE[1]
psightSE[2] <- pSE[2]
psightSE[3] <- pSE[3]
psightSE[4] <- pSE[4]
psightSE[5] <- pSE[5]
psightSE[6] <- pSE[6]
psightSE[7] <- pSE[7]
psightSE[8] <- pSE[8]
psightSE[9] <- pSE[9]
psightSE[10] <- pSE[10]
psightSE[11] <- pSE[11]
psightSE[12] <- pSE[12]
psightSE[13] <- pSE[13]
psightSE[14] <- pSE[14]
psightSE[15] <- pSE[15]
psightSE[16] <- pSE[16]
psightSE[17] <- pSE[17]
psightSE[18] <- pSE[18]
psightSE[19] <- pSE[19]
psightSE[20] <- pSE[20]
psightSE[21] <- pSE[21]
psightSE[22] <- pSE[22]
psightSE[23] <- pSE[23]
psightSE[24] <- pSE[24]
psightSE[25] <- pSE[25]
psightSE[26] <- pSE[26]
psightSE[27] <- pSE[27]
psightSE[28] <- pSE[28]
psightSE[29] <- pSE[29]
psightSE[30] <- pSE[30]
psightSE[31] <- pSE[31]
psightSE[32] <- pSE[32]
psightSE[33] <- pSE[33]
psightSE[34] <- pSE[34]
psightSE[35] <- pSE[35]
psightSE[36] <- pSE[36]
psightSE[37] <- pSE[37]
psightSE[38] <- pSE[38]
psightSE[39] <- pSE[39]
   
   # BR
   for (t in 1:19){
         logit(phiBR[t]) <- mu.phiBR[seasonphiBR[t]] + eps.phiBR[t]
         eps.phiBR[t] ~ dnorm(0, tau.phiBR[seasonphiBR[t]])T(-10,10)
   }
   
   for (t in 1:20){
         logit(pBR[t]) <- mu.pBR[seasonpBR[t]] + eps.pBR[t]
         eps.pBR[t] ~ dnorm(0, tau.pBR[seasonpBR[t]])T(-10,10)
   }

      for (s in 1:2){
         mean.phiBR[s] ~ dunif(0, 1)
         mu.phiBR[s] <- logit(mean.phiBR[s])
         sigma.phiBR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiBR[s] <- pow(sigma.phiBR[s], -2)
         sigma2.phiBR[s] <- pow(sigma.phiBR[s], 2)
      }
      
      for (s in 1:2){
         mean.pBR[s] ~ dunif(0, 1)
         mu.pBR[s] <- logit(mean.pBR[s])
         sigma.pBR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.pBR[s] <- pow(sigma.pBR[s], -2)
         sigma2.pBR[s] <- pow(sigma.pBR[s], 2)
      }

phitBR[1] <- 0
phitBR[2] <- 0
phitBR[3] <- phiBR[1]
phitBR[4] <- phiBR[2]
phitBR[5] <- 0
phitBR[6] <- 0
phitBR[7] <- phiBR[3]
phitBR[8] <- phiBR[4]
phitBR[9] <- 0
phitBR[10] <- 0
phitBR[11] <- phiBR[5]
phitBR[12] <- phiBR[6]
phitBR[13] <- 0
phitBR[14] <- 0
phitBR[15] <- phiBR[7]
phitBR[16] <- phiBR[8]
phitBR[17] <- 0
phitBR[18] <- 0
phitBR[19] <- phiBR[9]
phitBR[20] <- phiBR[10]
phitBR[21] <- 0
phitBR[22] <- 0
phitBR[23] <- phiBR[11]
phitBR[24] <- phiBR[12]
phitBR[25] <- 0
phitBR[26] <- 0
phitBR[27] <- phiBR[13]
phitBR[28] <- phiBR[14]
phitBR[29] <- 0
phitBR[30] <- 0
phitBR[31] <- phiBR[15]
phitBR[32] <- phiBR[16]
phitBR[33] <- 0
phitBR[34] <- 0
phitBR[35] <- phiBR[17]
phitBR[36] <- phiBR[18]
phitBR[37] <- 0
phitBR[38] <- 0
phitBR[39] <- phiBR[19]

psightBR[1] <- 0
psightBR[2] <- pBR[1]
psightBR[3] <- pBR[2]
psightBR[4] <- 0
psightBR[5] <- 0
psightBR[6] <- pBR[3]
psightBR[7] <- pBR[4]
psightBR[8] <- 0
psightBR[9] <- 0
psightBR[10] <- pBR[5]
psightBR[11] <- pBR[6]
psightBR[12] <- 0
psightBR[13] <- 0
psightBR[14] <- pBR[7]
psightBR[15] <- pBR[8]
psightBR[16] <- 0
psightBR[17] <- 0
psightBR[18] <- pBR[9]
psightBR[19] <- pBR[10]
psightBR[20] <- 0
psightBR[21] <- 0
psightBR[22] <- pBR[11]
psightBR[23] <- pBR[12]
psightBR[24] <- 0
psightBR[25] <- 0
psightBR[26] <- pBR[13]
psightBR[27] <- pBR[14]
psightBR[28] <- 0
psightBR[29] <- 0
psightBR[30] <- pBR[15]
psightBR[31] <- pBR[16]
psightBR[32] <- 0
psightBR[33] <- 0
psightBR[34] <- pBR[17]
psightBR[35] <- pBR[18]
psightBR[36] <- 0
psightBR[37] <- 0
psightBR[38] <- pBR[19]
psightBR[39] <- pBR[20]

      # AR
      for (t in 1:9){
         logit(phiAR[t]) <- mu.phiAR + eps.phiAR[t]
         eps.phiAR[t] ~ dnorm(0, tau.phiAR)T(-10,10)
      }
      
      for (t in 1:10){
         logit(pAR[t]) <- mu.pAR + eps.pAR[t]
         eps.pAR[t] ~ dnorm(0, tau.pAR)T(-10,10)
   }

         mean.phiAR ~ dunif(0, 1)
         mu.phiAR <- logit(mean.phiAR)
         sigma.phiAR ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiAR <- pow(sigma.phiAR, -2)
         sigma2.phiAR <- pow(sigma.phiAR, 2)

         mean.pAR ~ dunif(0, 1)
         mu.pAR <- logit(mean.pAR)
         sigma.pAR ~ dunif(0, 10)               # Prior for standard deviation
         tau.pAR <- pow(sigma.pAR, -2)
         sigma2.pAR <- pow(sigma.pAR, 2)

phitAR[1] <- 0
phitAR[2] <- 0
phitAR[3] <- 0
phitAR[4] <- phiAR[1]
phitAR[5] <- 0
phitAR[6] <- 0
phitAR[7] <- 0
phitAR[8] <- phiAR[2]
phitAR[9] <- 0
phitAR[10] <- 0
phitAR[11] <- 0
phitAR[12] <- phiAR[3]
phitAR[13] <- 0
phitAR[14] <- 0
phitAR[15] <- 0
phitAR[16] <- phiAR[4]
phitAR[17] <- 0
phitAR[18] <- 0
phitAR[19] <- 0
phitAR[20] <- phiAR[5]
phitAR[21] <- 0
phitAR[22] <- 0
phitAR[23] <- 0
phitAR[24] <- phiAR[6]
phitAR[25] <- 0
phitAR[26] <- 0
phitAR[27] <- 0
phitAR[28] <- phiAR[7]
phitAR[29] <- 0
phitAR[30] <- 0
phitAR[31] <- 0
phitAR[32] <- phiAR[8]
phitAR[33] <- 0
phitAR[34] <- 0
phitAR[35] <- 0
phitAR[36] <- phiAR[9]
phitAR[37] <- 0
phitAR[38] <- 0
phitAR[39] <- 0
         
psightAR[1] <- 0
psightAR[2] <- 0
psightAR[3] <- pAR[1]
psightAR[4] <- 0
psightAR[5] <- 0
psightAR[6] <- 0
psightAR[7] <- pAR[2]
psightAR[8] <- 0
psightAR[9] <- 0
psightAR[10] <- 0
psightAR[11] <- pAR[3]
psightAR[12] <- 0
psightAR[13] <- 0
psightAR[14] <- 0
psightAR[15] <- pAR[4]
psightAR[16] <- 0
psightAR[17] <- 0
psightAR[18] <- 0
psightAR[19] <- pAR[5]
psightAR[20] <- 0
psightAR[21] <- 0
psightAR[22] <- 0
psightAR[23] <- pAR[6]
psightAR[24] <- 0
psightAR[25] <- 0
psightAR[26] <- 0
psightAR[27] <- pAR[7]
psightAR[28] <- 0
psightAR[29] <- 0
psightAR[30] <- 0
psightAR[31] <- pAR[8]
psightAR[32] <- 0
psightAR[33] <- 0
psightAR[34] <- 0
psightAR[35] <- pAR[9]
psightAR[36] <- 0
psightAR[37] <- 0
psightAR[38] <- 0
psightAR[39] <- pAR[10]

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

# Define state-transition and reencounter probabilities - note no i index - is no longer individual 
   for (t in 1:(n.occasions-1)){
      psi[1,t,1] <- 0
      psi[1,t,2] <- phitDB[t] * psiDB[2,t]
      psi[1,t,3] <- phitDB[t] * psiDB[3,t]
      psi[1,t,4] <- phitDB[t] * psiDB[4,t]
      psi[1,t,5] <- phitDB[t] * psiDB[5,t]
      psi[1,t,6] <- 0
      psi[1,t,7] <- 0
      
      psi[2,t,1] <- 0
      psi[2,t,2] <- 0
      psi[2,t,3] <- 0
      psi[2,t,4] <- phitJB[t] * psiJB[4,t]
      psi[2,t,5] <- phitJB[t] * psiJB[5,t]
      psi[2,t,6] <- phitJB[t] * psiJB[6,t]
      psi[2,t,7] <- 0

      psi[3,t,1] <- 0
      psi[3,t,2] <- 0
      psi[3,t,3] <- 0
      psi[3,t,4] <- phitMI[t] * psiMI[4,t]
      psi[3,t,5] <- phitMI[t] * psiMI[5,t]
      psi[3,t,6] <- phitMI[t] * psiMI[6,t]
      psi[3,t,7] <- 0
      
      psi[4,t,1] <- 0
      psi[4,t,2] <- 0
      psi[4,t,3] <- 0
      psi[4,t,4] <- phitCC[t] * psiCC[4,t]
      psi[4,t,5] <- phitCC[t] * psiCC[5,t]
      psi[4,t,6] <- phitCC[t] * psiCC[6,t]
      psi[4,t,7] <- phitCC[t] * psiCC[7,t]
      
      psi[5,t,1] <- phitSE[t] * psiSE[1,t]
      psi[5,t,2] <- phitSE[t] * psiSE[2,t]
      psi[5,t,3] <- phitSE[t] * psiSE[3,t]
      psi[5,t,4] <- phitSE[t] * psiSE[4,t]
      psi[5,t,5] <- phitSE[t] * psiSE[5,t]
      psi[5,t,6] <- phitSE[t] * psiSE[6,t]
      psi[5,t,7] <- phitSE[t] * psiSE[7,t]
      
      psi[6,t,1] <- phitBR[t] * psiBR[1,t]
      psi[6,t,2] <- 0
      psi[6,t,3] <- 0
      psi[6,t,4] <- 0
      psi[6,t,5] <- phitBR[t] * psiBR[5,t]
      psi[6,t,6] <- phitBR[t] * psiBR[6,t]
      psi[6,t,7] <- phitBR[t] * psiBR[7,t]
      
      psi[7,t,1] <- phitAR[t] * psiAR[1,t]
      psi[7,t,2] <- 0
      psi[7,t,3] <- 0
      psi[7,t,4] <- 0
      psi[7,t,5] <- phitAR[t] * psiAR[5,t]
      psi[7,t,6] <- 0
      psi[7,t,7] <- 0
      
      po[1,t] <- psightDB[t]
      po[2,t] <- psightJB[t]
      po[3,t] <- psightMI[t]
      po[4,t] <- psightCC[t]
      po[5,t] <- psightSE[t]
      po[6,t] <- psightBR[t]
      po[7,t] <- psightAR[t]


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
                  transAR = rep(c(0, 0, 0, 1), length.out = n.occasions-1),
                  seasonphiCC = rep(c(1, 2), length.out = 20),
                  seasonphiSE = rep(c(1, 2, 3, 4), length.out = n.occasions-1),
                  seasonphiBR = rep(c(1, 2), length.out = 19),
                  seasonpCC = rep(c(1, 2), length.out = 20),
                  seasonpSE = rep(c(1, 2, 3, 4), length.out = n.occasions-1),
                  seasonpBR = rep(c(1, 2), length.out = 20))

# seasonDBphi = rep(c(1,0,0,0), length.out = n.occasions-1),
# seasonJBphi = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonMIphi = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonCCphi = rep(c(0,2,3,0), length.out = n.occasions-1),
# seasonSEphi = rep(c(1,2,3,4), length.out = n.occasions-1),
# seasonBRphi = rep(c(0,0,3,4), length.out = n.occasions-1),
# seasonARphi = rep(c(0,0,0,4), length.out = n.occasions-1),
# seasonDBp = rep(c(1,0,0,0), length.out = n.occasions-1),
# seasonJBp = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonMIp = rep(c(0,2,0,0), length.out = n.occasions-1),
# seasonCCp = rep(c(0,2,3,0), length.out = n.occasions-1),
# seasonSEp = rep(c(1,2,3,4), length.out = n.occasions-1),
# seasonBRp = rep(c(0,0,3,4), length.out = n.occasions-1),
# seasonARp = rep(c(0,0,0,4), length.out = n.occasions-1)

# Initial values 
inits <- function(){list(mean.phiDB = runif(1, 0, 1), mean.phiJB = runif(1, 0, 1), mean.phiMI = runif(1, 0, 1),
                         mean.phiCC = runif(2, 0, 1), mean.phiSE = runif(4, 0, 1), mean.phiBR = runif(2, 0, 1), mean.phiAR = runif(1, 0, 1),
                         sigma.phiDB = runif(1, 0, 10), sigma.phiJB = runif(1, 0, 10), sigma.phiMI = runif(1, 0, 10),
                         sigma.phiCC = runif(2, 0, 10), sigma.phiSE = runif(4, 0, 10), sigma.phiBR = runif(2, 0, 10), sigma.phiAR = runif(1, 0, 10),
                         mean.pDB = runif(1, 0, 1), mean.pJB = runif(1, 0, 1), mean.pMI = runif(1, 0, 1),
                         mean.pCC = runif(2, 0, 1), mean.pSE = runif(4, 0, 1), mean.pBR = runif(2, 0, 1), mean.pAR = runif(1, 0, 1),
                         sigma.pDB = runif(1, 0, 10), sigma.pJB = runif(1, 0, 10), sigma.pMI = runif(1, 0, 10),
                         sigma.pCC = runif(2, 0, 10), sigma.pSE = runif(4, 0, 10), sigma.pBR = runif(2, 0, 10), sigma.pAR = runif(1, 0, 10))}  


# Parameters monitored
parameters <- c("mean.phiDB", "mean.phiJB", "mean.phiMI", "mean.phiCC", "mean.phiSE", "mean.phiBR", "mean.phiAR",
                "mean.pDB", "mean.pJB", "mean.pMI", "mean.pCC", "mean.pSE", "mean.pBR", "mean.pAR",
                "psiDB", "psiJB", "psiMI", "psiCC", "psiSE", "psiBR", "psiAR",
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

saveRDS(sim.marr$summary, file = paste0("./analysis-output/seasonal-scenario-3-summary", Sys.Date(), ".rds"))

saveRDS(sim.marr$sims.list, file = paste0("./analysis-output/seasonal-scenario-3-simslist", Sys.Date(), ".rds"))

# parameter identifiability checks
sims.list <- sim.marr$sims.list
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
   
   png(filename = paste0("analysis-output/parameter-identifiability/seasonal-scenario-2/",
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

png(filename = "figures/seasonal-scenario-2.ppcheck.png", width = 8, height = 6,
    units = "in", res = 600)

print(ppcheck)

dev.off()

# create dataframes of set values
# survival
phiPreB.sim <- tibble(site = c("DB-PreB", "JB-PreB", "MI-PreB"),
                      value = c(mean.phiDB1, mean.phiJB1, mean.phiMI1))

phiPostB1.sim <- tibble(site = c("DB-PostB1", "JB-PostB1", "MI-PostB1"),
                        value = c(mean.phiDB2, mean.phiJB2, mean.phiMI2))

phiPostB2.sim <- tibble(site = c("DB-PostB2", "JB-PostB2", "MI-PostB2"),
                        value = c(mean.phiDB3, mean.phiJB3, mean.phiMI3))

phiNonB.sim <- tibble(site = c("DB-NonB", "JB-NonB", "MI-NonB"),
                      value = c(mean.phiDB4, mean.phiJB4, mean.phiMI4))

# resighting
pPreB.sim <- tibble(site = c("DB-PreB", "JB-PreB", "MI-PreB"),
                    value = c(mean.pDB1, mean.pJB1, mean.pMI1))

pPostB1.sim <- tibble(site = c("DB-PostB1", "JB-PostB1", "MI-PostB1"),
                      value = c(mean.pDB2, mean.pJB2, mean.pMI2))

pPostB2.sim <- tibble(site = c("DB-PostB2", "JB-PostB2", "MI-PostB2"),
                      value = c(mean.pDB3, mean.pJB3, mean.pMI3))

pNonB.sim <- tibble(site = c("DB-NonB", "JB-NonB", "MI-NonB"),
                    value = c(mean.pDB4, mean.pJB4, mean.pMI4))

# transition
psiDB.sim <- tibble(site = c("DB-DB", "DB-JB", "DB-MI"),
                    value = c(psiDB.DB, psiDB.JB, psiDB.MI))

psiJB.sim <- tibble(site = c("JB-DB", "JB-JB", "JB-MI"),
                    value = c(psiJB.DB, psiJB.JB, psiJB.MI))

psiMI.sim <- tibble(site = c("MI-DB", "MI-JB", "MI-MI"),
                    value = c(psiMI.DB, psiMI.JB, psiMI.MI))

# format survival
phiPreB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.1, mean.phi.2, mean.phi.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.1", "DB-PreB")) %>%
   mutate(site = str_replace(site, "mean.phi.2", "JB-PreB")) %>%
   mutate(site = str_replace(site, "mean.phi.3", "MI-PreB"))

phiPostB1.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.4, mean.phi.5, mean.phi.6) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.4", "DB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.phi.5", "JB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.phi.6", "MI-PostB1"))

phiPostB2.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.7, mean.phi.8, mean.phi.9) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.7", "DB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.phi.8", "JB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.phi.9", "MI-PostB2"))

phiNonB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.phi.10, mean.phi.11, mean.phi.12) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.phi.10", "DB-NonB")) %>%
   mutate(site = str_replace(site, "mean.phi.11", "JB-NonB")) %>%
   mutate(site = str_replace(site, "mean.phi.12", "MI-NonB"))

# plot survival
phiPreB.plot <- ggplot() +
   geom_violin(phiPreB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPreB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiPostB1.plot <- ggplot() +
   geom_violin(phiPostB1.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPostB1.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiPostB2.plot <- ggplot() +
   geom_violin(phiPostB2.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiPostB2.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

phiNonB.plot <- ggplot() +
   geom_violin(phiNonB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(phiNonB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

# format resighting
pPreB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.1, mean.p.2, mean.p.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.1", "DB-PreB")) %>%
   mutate(site = str_replace(site, "mean.p.2", "JB-PreB")) %>%
   mutate(site = str_replace(site, "mean.p.3", "MI-PreB"))

pPostB1.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.4, mean.p.5, mean.p.6) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.4", "DB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.p.5", "JB-PostB1")) %>%
   mutate(site = str_replace(site, "mean.p.6", "MI-PostB1"))

pPostB2.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.7, mean.p.8, mean.p.9) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.7", "DB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.p.8", "JB-PostB2")) %>%
   mutate(site = str_replace(site, "mean.p.9", "MI-PostB2"))

pNonB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(mean.p.10, mean.p.11, mean.p.12) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>% 
   mutate(site = str_replace(site, "mean.p.10", "DB-NonB")) %>%
   mutate(site = str_replace(site, "mean.p.11", "JB-NonB")) %>%
   mutate(site = str_replace(site, "mean.p.12", "MI-NonB"))

# plot resighting
pPreB.plot <- ggplot() +
   geom_violin(pPreB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPreB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pPostB1.plot <- ggplot() +
   geom_violin(pPostB1.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPostB1.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pPostB2.plot <- ggplot() +
   geom_violin(pPostB2.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pPostB2.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

pNonB.plot <- ggplot() +
   geom_violin(pNonB.mod, mapping = aes(x = site, y = estimate, group = site, fill = site), alpha = 0.6) +
   geom_point(pNonB.sim, mapping = aes(x = site, y = value), size = 2, colour = "black") +
   scale_fill_viridis(discrete = TRUE) +
   #coord_flip() +
   xlab("Site") +
   ylab("Survival probability") +
   theme(legend.position = "none")

# format transitions
psiDB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiDB.1, psiDB.2, psiDB.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiDB.1", "DB-DB")) %>%
   mutate(site = str_replace(site, "psiDB.2", "DB-JB")) %>%
   mutate(site = str_replace(site, "psiDB.3", "DB-MI"))

psiJB.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiJB.1, psiJB.2, psiJB.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiJB.1", "JB-DB")) %>%
   mutate(site = str_replace(site, "psiJB.2", "JB-JB")) %>%
   mutate(site = str_replace(site, "psiJB.3", "JB-MI"))

psiMI.mod <- sims.list %>% 
   as.data.frame() %>% 
   select(psiMI.1, psiMI.2, psiMI.3) %>% 
   pivot_longer(cols = 1:3, names_to = "site", values_to = "estimate") %>%
   mutate(site = str_replace(site, "psiMI.1", "MI-DB")) %>%
   mutate(site = str_replace(site, "psiMI.2", "MI-JB")) %>%
   mutate(site = str_replace(site, "psiMI.3", "MI-MI"))

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

# plot all
# survival
png(filename = "figures/seasonal-scenario-1.phiPreB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPreB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiPostB1.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPostB1.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiPostB2.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiPostB2.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.phiNonB.png", width = 8, height = 8,
    units = "in", res = 600)

print(phiNonB.plot)

dev.off()

# resighting
png(filename = "figures/seasonal-scenario-1.pPreB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPreB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pPostB1.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPostB1.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pPostB2.png", width = 8, height = 8,
    units = "in", res = 600)

print(pPostB2.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.pNonB.png", width = 8, height = 8,
    units = "in", res = 600)

print(pNonB.plot)

dev.off()

# transition
png(filename = "figures/seasonal-scenario-1.psiDB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiDB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.psiJB.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiJB.plot)

dev.off()

png(filename = "figures/seasonal-scenario-1.psiMI.png", width = 8, height = 8,
    units = "in", res = 600)

print(psiMI.plot)

dev.off()