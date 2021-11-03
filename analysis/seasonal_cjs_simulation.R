################################################################################
# SEASONAL SURVIVAL SIMULATION SCENARIOS - CJS MODEL
#
# using m-array format
#
# single state seasonal survival before moving to multistate
################################################################################

library(rjags)
library(jagsUI)
library(tidyverse)
library(cowplot)
library(viridis)

set.seed(42)
theme_set(theme_bw())

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

# Random time effects
# Define parameter values
n.occasions <- 30                  # Number of capture occasions
n.years <- 10
marked <- rep(c(1000, 400, 50), length.out = n.occasions-1)   # Annual number of newly marked individuals
mean.phi1 <- 0.8
mean.phi2 <- 0.9
mean.phi3 <- 0.95
var.phi <- 0.3                       # Temporal variance of survival

mean.p1 <- 0.7
mean.p2 <- 0.5
mean.p3 <- 0.4
var.p <- 0.5

# Determine annual survival probabilities
logit.phi1 <- rnorm(n.years, qlogis(mean.phi1), var.phi^0.5)
phi1 <- plogis(logit.phi1)

logit.phi2 <- rnorm(n.years, qlogis(mean.phi2), var.phi^0.5)
phi2 <- plogis(logit.phi2)

logit.phi3 <- rnorm(n.years, qlogis(mean.phi3), var.phi^0.5)
phi3 <- plogis(logit.phi3)

phi <- c(rbind(phi1, phi2, phi3))
phi <- phi[1:29]

# Determine annual resighting probabilities
logit.p1 <- rnorm(n.years, qlogis(mean.p1), var.p^0.5)
p1 <- plogis(logit.p1)

logit.p2 <- rnorm(n.years, qlogis(mean.p2), var.p^0.5)
p2 <- plogis(logit.p2)

logit.p3 <- rnorm(n.years, qlogis(mean.p3), var.p^0.5)
p3 <- plogis(logit.p3)


p <- c(rbind(p1, p2, p3))
p <- p[1:29]

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

CH.marr <- marray(CH)

# # Create vector with occasion of marking
# get.first <- function(x) min(which(x!=0))
# f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-seasonal-simulation.jags")
cat("
model {

# Priors and constraints
for (t in 1:(n.occasions-1)){
   logit(phi[t]) <- mu.phi[season[t]] + eps.phi[t]
   eps.phi[t] ~ dnorm(0, tau.phi[season[t]])T(-10,10)
   logit(p[t]) <- mu.p[season[t]] + eps.p[t] 
   eps.p[t] ~ dnorm(0, tau.p[season[t]])T(-10,10)
   }

for(s in 1:3){
   mean.phi[s] ~ dunif(0, 1)             # Prior for mean survival
   mu.phi[s] <- logit(mean.phi[s])       # Logit transformation
   sigma.phi[s] ~ dunif(0, 10)               # Prior for standard deviation
   tau.phi[s] <- pow(sigma.phi[s], -2)
   sigma2.phi[s] <- pow(sigma.phi[s], 2)
   
   mean.p[s] ~ dunif(0, 1)             # Prior for mean resighting
   mu.p[s] <- logit(mean.p[s])         # Logit transformation
   sigma.p[s] ~ dunif(0, 10)               # Prior for standard deviation
   tau.p[s] <- pow(sigma.p[s], -2)
   sigma2.p[s] <- pow(sigma.p[s], 2)
}

# Temporal variance on real scale
# sigma2.real <- sigma2 * pow(mean.phi, 2) * pow((1-mean.phi), 2) 


# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
   marr[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
   }

# Define the cell probabilities of the m-array:
# Main diagonal
for (t in 1:(n.occasions-1)){
   q[t] <- 1-p[t]
   pr[t,t] <- phi[t]*p[t]	
   # Above main diagonal
   for (j in (t+1):(n.occasions-1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j	
   # Below main diagonal
   for (j in 1:(t-1)){
      pr[t,j]<-0
      } #j
   } #t

# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
   pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
   } # t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      }
   }
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      }
   }
fit <- sum(E.org[,])
fit.new <- sum(E.new[,])

}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr = CH.marr, n.occasions = dim(CH.marr)[2], rel = rowSums(CH.marr),
                  season = rep(c(1,2,3), length.out = n.occasions-1))

# Initial values
inits <- function(){list(mean.phi = runif(3, 0, 1), sigma.phi = runif(3, 0, 5),
                         mean.p = runif(3, 0, 1), sigma.p = runif(3, 0, 5))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "phi", "p", "sigma2.phi", "sigma2.p", "fit", "fit.new")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 7000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.seasonal <- jags(jags.data, inits, parameters, "cjs-seasonal-simulation.jags", n.chains = nc, n.thin = nt,
                     n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.seasonal, digits = 3)

sims.list <- cjs.seasonal$sims.list
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

################################################################################
# SEASONAL SURVIVAL SIMULATION SCENARIOS - CJS MODEL WITH TRANSIENCE
#
# using m-array format
#
# single state seasonal survival before moving to multistate
################################################################################

# Random time effects
# Define parameter values
n.occasions <- 30                  # Number of capture occasions
n.years <- 10
marked <- rep(c(1000, 400, 50), length.out = n.occasions-1)   # Seasonal number of first resights

mean.phi1 <- 0.8
mean.phi2 <- 0.9
mean.phi3 <- 0.95
var.phi <- 0.3                       # Temporal variance of survival

mean.phi.t1 <- 0.5
mean.phi.t2 <- 0.6
mean.phi.t3 <- 0.65
var.phi.t <- 0.3

mean.p1 <- 0.7
mean.p2 <- 0.5
mean.p3 <- 0.4
var.p <- 0.5

# Determine annual survival probabilities
logit.phi1 <- rnorm(n.years, qlogis(mean.phi1), var.phi^0.5)
phi1 <- plogis(logit.phi1)

logit.phi2 <- rnorm(n.years, qlogis(mean.phi2), var.phi^0.5)
phi2 <- plogis(logit.phi2)

logit.phi3 <- rnorm(n.years, qlogis(mean.phi3), var.phi^0.5)
phi3 <- plogis(logit.phi3)

phi <- c(rbind(phi1, phi2, phi3))
phi <- phi[1:29]

logit.phi.t1 <- rnorm(n.years, qlogis(mean.phi.t1), var.phi.t^0.5)
phi.t1 <- plogis(logit.phi.t1)

logit.phi.t2 <- rnorm(n.years, qlogis(mean.phi.t2), var.phi.t^0.5)
phi.t2 <- plogis(logit.phi.t2)

logit.phi.t3 <- rnorm(n.years, qlogis(mean.phi.t3), var.phi.t^0.5)
phi.t3 <- plogis(logit.phi.t3)

phi.t <- c(rbind(phi.t1, phi.t2, phi.t3))
phi.t <- phi.t[1:29]

# Determine annual resighting probabilities
logit.p1 <- rnorm(n.years, qlogis(mean.p1), var.p^0.5)
p1 <- plogis(logit.p1)

logit.p2 <- rnorm(n.years, qlogis(mean.p2), var.p^0.5)
p2 <- plogis(logit.p2)

logit.p3 <- rnorm(n.years, qlogis(mean.p3), var.p^0.5)
p3 <- plogis(logit.p3)

p <- c(rbind(p1, p2, p3))
p <- p[1:29]

# Define matrices with survival and recapture probabilities
PHI <- matrix(0, ncol = n.occasions-1, nrow = sum(marked))
for (i in 1:(length(marked)-1)){
  #PHI[(sum(marked[1:i])-marked[i]+1):sum(marked[1:i]),i:(n.occasions-1)] <- matrix(rep(phi.t[1:(n.occasions-i)],marked[i]), ncol = n.occasions-i, byrow = TRUE)
  PHI[(sum(marked[1:i])-marked[i]+1):sum(marked[1:i]),i] <- matrix(rep(phi.t[i],marked[i]), ncol = 1, byrow = TRUE)
  PHI[(sum(marked[1:i])-marked[i]+1):sum(marked[1:i]),(i+1):(n.occasions-1)] <- matrix(rep(phi[i+1:(n.occasions-(i+1))],marked[i]), ncol = n.occasions-(i+1), byrow = TRUE)
}
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)

# Apply simulation function
CH <- simul.cjs(PHI, P, marked) 

# Create separate m-arrays for first capture and subsequent recaptures
cap <- apply(CH, 1, sum)
ind <- which(cap >= 2)
CH.R <- CH[ind,] # First period CH recaptured at least once
CH.N <- CH[-ind,] # First period CH never recaptured
# Remove first capture
first <- numeric()
for (i in 1:dim(CH.R)[1]){
  first[i] <- min(which(CH.R[i,]==1))
}
CH.R1 <- CH.R
for (i in 1:dim(CH.R)[1]){
  CH.R1[i,first[i]] <- 0
}
# Create m-array of those recaptured at least once
CH.marray <- marray(CH.R1)
# Create CH matrix for first period, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.R1)[1]){
  second[i] <- min(which(CH.R1[i,]==1))
}
CH.R2 <- matrix(0, nrow = dim(CH.R)[1], ncol = dim(CH.R)[2])
for (i in 1:dim(CH.R)[1]){
  CH.R2[i,first[i]] <- 1
  CH.R2[i,second[i]] <- 1
}
# Create m-array for these
CH.R.marray <- marray(CH.R2)
# The last column ought to show the number of transients not recaptured
# again and should all be zeros, since all of them are released as "residents"
CH.R.marray[,dim(CH)[2]] <- 0
# Create the m-array for transients never recaptured and add it to the
# previous m-array
CH.N.marray <- marray(CH.N)
CH.T.marray <- CH.R.marray + CH.N.marray

# Specify model in BUGS language

sink("cjs-seasonal-sim.jags")
cat("
model {
# Priors and constraints
for (t in 1:(n.occasions-1)){
logit(phi.trans[t]) <- mu.t[season[t]] + epsilon.t[t]
epsilon.t[t] ~ dnorm(0, tau.t[season[t]])T(-15,15) # Range restriction
phi.t[t] <- 1/(1+exp(-mu.t[season[t]]-epsilon.t[t]))

logit(phi.app[t]) <- mu.phi[season[t]] + epsilon.phi[t]
epsilon.phi[t] ~ dnorm(0, tau.phi[season[t]])T(-15,15) # Range restriction
phi[t] <- 1/(1+exp(-mu.phi[season[t]]-epsilon.phi[t]))

logit(psight[t]) <- mu.p[season[t]] + epsilon.p[t]
epsilon.p[t] ~ dnorm(0, tau.p[season[t]])T(-15,15) # Range restriction
p[t] <- 1/(1+exp(-mu.p[season[t]]-epsilon.p[t]))

}

for (s in 1:3){

mu.t[s] <- log(mean.phi.t[s] / (1-mean.phi.t[s]))
mean.phi.t[s] ~ dunif(0, 1) # Prior for mean transience
sigma.t[s] ~ dunif(0, 5) # Prior on sd of temp. var
tau.t[s] <- pow(sigma.t[s], -2)
sigma.t2[s] <- pow(sigma.t[s], 2)

mu.phi[s] <- log(mean.phi[s] / (1-mean.phi[s]))
mean.phi[s] ~ dunif(0, 1) # Prior for mean survival
sigma.phi[s] ~ dunif(0, 5) # Prior on sd of temp. var
tau.phi[s] <- pow(sigma.phi[s], -2)
sigma2.phi[s] <- pow(sigma.phi[s], 2)

mu.p[s] <- log(mean.p[s] / (1-mean.p[s]))
mean.p[s] ~ dunif(0, 1) # Prior for mean survival
sigma.p[s] ~ dunif(0, 5) # Prior on sd of temp. var
tau.p[s] <- pow(sigma.p[s], -2)
sigma2.p[s] <- pow(sigma.p[s], 2)

}

# Define the multinomial likelihood
for (t in 1:(n.occasions-1)){
marr.t[t,1:n.occasions] ~ dmulti(pr.t[t,], rel.t[t])
marr[t,1:n.occasions] ~ dmulti(pr[t,], rel[t])
}

# Define the cell probabilities of the m-arrays
# Main diagonal
for (t in 1:(n.occasions-1)){
q[t] <- 1-p[t] # Probability of non-recapture
pr.t[t,t] <- phi.trans[t]*p[t]
pr[t,t] <- phi[t]*p[t]

# Above main diagonal
for (j in (t+1):(n.occasions-1)){
pr.t[t,j] <- phi.trans[t]*prod(phi[(t+1):j])*prod(q[t:(j-1)])*p[j]
pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
} # j

# Below main diagonal
for (j in 1:(t-1)){
pr.t[t,j] <- 0
pr[t,j] <- 0
} # j
} # t

# Last column: probability of non-recapture
for (t in 1:(n.occasions-1)){
pr.t[t,n.occasions] <- 1-sum(pr.t[t,1:(n.occasions-1)])
pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
} # t

# Assess model fit using Freeman-Tukey statistic
# Compute fit statistics for observed data
for (t in 1:(n.occasions-1)){
   for (j in 1:n.occasions){
      expmarr[t,j] <- rel[t]*pr[t,j]
      expmarr.t[t,j] <- rel.t[t]*pr.t[t,j]
      E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      E.org.t[t,j] <- pow((pow(marr.t[t,j], 0.5)-pow(expmarr.t[t,j], 0.5)), 2)
      } #j
   } #t
   
# Generate replicate data and compute fit stats from them
for (t in 1:(n.occasions-1)){
   marr.new[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
   marr.t.new[t,1:n.occasions] ~ dmulti(pr.t[t, ], rel.t[t])
   for (j in 1:n.occasions){
      E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      E.t.new[t,j] <- pow((pow(marr.t.new[t,j], 0.5)-pow(expmarr.t[t,j], 0.5)), 2)
      } #j
   } #t
   
fit <- sum(E.org[,])
fit.t <- sum(E.org.t[,])
fit.new <- sum(E.new[,])
fit.t.new <- sum(E.t.new[,])
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(marr.t = CH.T.marray, marr = CH.marray, n.occasions = dim(CH)[2], rel.t = rowSums(CH.T.marray),
                  rel = rowSums(CH.marray), season = rep(c(1,2,3), length.out = dim(CH.marray)[2]-1))

# Initial values
inits <- function(){list(mean.phi = runif(3, 0, 1), sigma.phi = runif(3, 0, 5),
                         mean.phi.t = runif(3, 0, 1), sigma.phi.t = runif(3, 0, 5),
                         mean.p = runif(3, 0, 1), sigma.p = runif(3, 0, 5))}  

# Parameters monitored
parameters <- c("mean.phi.t", "mean.phi", "mean.p", "phi.t", "phi", "p",
                "sigma.2t", "sigma2.phi", "sigma2.p", "fit.t", "fit.t.new", "fit", "fit.new")

# MCMC settings
ni <- 500000
nt <- 30
nb <- 200000
nc <- 3

# Call JAGS from R (BRT 3 min)
cjs.seasonal.sim <- jags(jags.data, inits, parameters, "cjs-seasonal-sim.jags", n.chains = nc, n.thin = nt,
                     n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(cjs.seasonal.sim, digits = 3)

sims.list <- cjs.seasonal.sim$sims.list
sims.list <- as.data.frame(sims.list)

# evaluation of fit
mean(sims.list$fit.new > sims.list$fit)

ppcheck <- ggplot(sims.list, aes(x = fit, y = fit.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-test/ppcheck.png", width = 8, height = 8,
    units = "in", res = 600)

print(ppcheck)

dev.off()

mean(sims.list$fit.t.new > sims.list$fit.t)

ppcheck.t <- ggplot(sims.list, aes(x = fit.t, y = fit.t.new)) +
  geom_point(alpha = 0.3, pch = 16) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Discrepancy actual data") +
  ylab("Discrepancy ideal data") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-test/ppcheck.t.png", width = 8, height = 8,
    units = "in", res = 600)

print(ppcheck.t)

dev.off()

saveRDS(cjs.seasonal.sim$summary, "./analysis-output/cjs-seasonal-sim-summary.rds")
write.csv(cjs.seasonal.sim$summary, "./analysis-output/cjs-seasonal-sim-summary.csv")

saveRDS(cjs.seasonal.sim$sims.list, "./analysis-output/cjs-seasonal-sim-simslist.rds")

################################################################################

sims.list <- readRDS("./analysis-output/cjs-seasonal-sim-simslist.rds")
sims.list <- as.data.frame(sims.list)

# create dataframes of set values
mean.phi.sim <- tibble(season = c("prebreeding", "postbreeding", "nonbreeding"),
                       value = c(mean.phi1, mean.phi2, mean.phi3))

phi.sim <- tibble(prebreeding = phi1,
                  postbreeding = phi2,
                  nonbreeding = phi3) %>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "value") %>% 
  mutate(year = rep(1:10, each = 3))

mean.t.sim <- tibble(season = c("prebreeding", "postbreeding", "nonbreeding"),
                     value = c(mean.phi.t1, mean.phi.t2, mean.phi.t3))

t.sim <- tibble(prebreeding = phi.t1,
                postbreeding = phi.t2,
                nonbreeding = phi.t3)%>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "value") %>% 
  mutate(year = rep(1:10, each = 3))

mean.p.sim <- tibble(season = c("prebreeding", "postbreeding", "nonbreeding"),
                     value = c(mean.p1, mean.p2, mean.p3))

p.sim <- tibble(prebreeding = p1,
                postbreeding = p2,
                nonbreeding = p3)%>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "value") %>% 
  mutate(year = rep(1:10, each = 3))

# format output
mean.phi.mod <- sims.list %>% 
  select(mean.phi.1, mean.phi.2, mean.phi.3) %>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "estimate") %>% 
  mutate(season = str_replace(season, "mean.phi.1", "prebreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.2", "postbreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.3", "nonbreeding"))

mean.t.mod <- sims.list %>% 
  select(mean.phi.t.1, mean.phi.t.2, mean.phi.t.3) %>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "estimate") %>% 
  mutate(season = str_replace(season, "mean.phi.t.1", "prebreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.t.2", "postbreeding")) %>%
  mutate(season = str_replace(season, "mean.phi.t.3", "nonbreeding"))

mean.p.mod <- sims.list %>% 
  select(mean.p.1, mean.p.2, mean.p.3) %>% 
  pivot_longer(cols = 1:3, names_to = "season", values_to = "estimate") %>% 
  mutate(season = str_replace(season, "mean.p.1", "prebreeding")) %>%
  mutate(season = str_replace(season, "mean.p.2", "postbreeding")) %>%
  mutate(season = str_replace(season, "mean.p.3", "nonbreeding"))

phi.mod <- sims.list %>% 
  select(phi.1, phi.2, phi.3, phi.4, phi.5, phi.6, phi.7, phi.8, phi.9, phi.10,
         phi.11, phi.12, phi.13, phi.14, phi.15, phi.16, phi.17, phi.18, phi.19, phi.20,
         phi.21, phi.22, phi.23, phi.24, phi.25, phi.26, phi.27, phi.28, phi.29) %>% 
  pivot_longer(cols = 1:29, names_to = "parameter", values_to = "estimate") %>% 
  arrange(match(parameter,
                c("phi.1", "phi.2", "phi.3", "phi.4", "phi.5", "phi.6", "phi.7", "phi.8", "phi.9", "phi.10",
                  "phi.11", "phi.12", "phi.13", "phi.14", "phi.15", "phi.16", "phi.17", "phi.18", "phi.19", "phi.20",
                  "phi.21", "phi.22", "phi.23", "phi.24", "phi.25", "phi.26", "phi.27", "phi.28", "phi.29")),
          desc(estimate)) %>% 
  mutate(season = rep(c("prebreeding", "postbreeding", "nonbreeding"), each = 30000, length.out = 870000)) %>%   
  mutate(year = rep(1:10, each = 90000, length.out = 870000))

# make plots
mean.phi.plot <- ggplot() +
  geom_violin(mean.phi.mod, mapping = aes(x = season, y = estimate, group = season, fill = season), alpha = 0.6) +
  geom_boxplot(mean.phi.mod, mapping = aes(x = season, y = estimate, group = season),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(mean.phi.mod, mapping = aes(x = season, y = estimate, group = season),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(mean.phi.sim, mapping = aes(x = season, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Season") +
  ylab("Survival probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-test/mean.phi.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.phi.plot)

dev.off()

mean.t.plot <- ggplot() +
  geom_violin(mean.t.mod, mapping = aes(x = season, y = estimate, group = season, fill = season), alpha = 0.6) +
  geom_boxplot(mean.t.mod, mapping = aes(x = season, y = estimate, group = season),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(mean.t.mod, mapping = aes(x = season, y = estimate, group = season),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(mean.t.sim, mapping = aes(x = season, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Season") +
  ylab("Survival/transience probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-test/mean.t.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.t.plot)

dev.off()

mean.p.plot <- ggplot() +
  geom_violin(mean.p.mod, mapping = aes(x = season, y = estimate, group = season, fill = season), alpha = 0.6) +
  geom_boxplot(mean.p.mod, mapping = aes(x = season, y = estimate, group = season),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(mean.p.mod, mapping = aes(x = season, y = estimate, group = season),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(mean.p.sim, mapping = aes(x = season, y = value), size = 2, colour = "black") +
  scale_fill_viridis(discrete = TRUE) +
  #coord_flip() +
  xlab("Season") +
  ylab("Resighting probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-test/mean.p.png", width = 8, height = 8,
    units = "in", res = 600)

print(mean.p.plot)

dev.off()

phi.plot <- ggplot() +
  geom_violin(phi.mod, mapping = aes(x = year, y = estimate, group = year, fill = season), alpha = 0.6) +
  geom_boxplot(phi.mod, mapping = aes(x = year, y = estimate, group = year),
               width = .05, fill = NA, outlier.colour = NA) +
  stat_summary(phi.mod, mapping = aes(x = year, y = estimate, group = year),
               fun = median, geom = "point", fill = "white", shape = 21, size = 2.5) +
  geom_point(phi.sim, mapping = aes(x = year, y = value), size = 2, colour = "black") +
  facet_wrap(. ~ season, nrow = 3) +
  scale_fill_viridis(discrete = TRUE) +
  scale_x_continuous(breaks = 1:10) +
  #coord_flip() +
  xlab("Year") +
  ylab("Survival probability") +
  theme(legend.position = "none")

png(filename = "figures/seasonal-cjs-test/phi.png", width = 8, height = 5,
    units = "in", res = 600)

print(phi.plot)

dev.off()
