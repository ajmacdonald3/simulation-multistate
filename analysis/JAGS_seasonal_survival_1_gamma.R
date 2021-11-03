sink("seasonal_survival_1.jags")
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
   phiDB ~ dunif(0, 1)
   phiJB ~ dunif(0, 1)
   phiMI ~ dunif(0, 1)
   phiAR ~ dunif(0, 1)
   phiUN ~ dunif(0, 1)
   
   pDB[1] <- 0
   pDB[2] <- 0
   pDB[3] <- mean.pDB
   pDB[4] <- 0
   pDB[5] <- 0
   pDB[6] <- mean.pDB
   pDB[7] <- 0
   pDB[8] <- 0
   pDB[9] <- mean.pDB
   pDB[10] <- 0
   pDB[11] <- 0
   pDB[12] <- mean.pDB
   pDB[13] <- 0
   pDB[14] <- 0
   pDB[15] <- mean.pDB
   
   mean.pDB ~ dunif(0, 1)
   
   pJB[1] ~ dunif(0, 1)
   pJB[2] <- 0
   pJB[3] <- 0
   pJB[4] ~ dunif(0, 1)
   pJB[5] <- 0
   pJB[6] <- 0
   pJB[7] ~ dunif(0, 1)
   pJB[8] <- 0
   pJB[9] <- 0
   pJB[10] ~ dunif(0, 1)
   pJB[11] <- 0
   pJB[12] <- 0
   pJB[13] ~ dunif(0, 1)
   pJB[14] <- 0
   pJB[15] <- 0
   
   pMI[1] <- mean.pMI
   pMI[2] <- 0
   pMI[3] <- 0
   pMI[4] <- mean.pMI
   pMI[5] <- 0
   pMI[6] <- 0
   pMI[7] <- mean.pMI
   pMI[8] <- 0
   pMI[9] <- 0
   pMI[10] <- mean.pMI
   pMI[11] <- 0
   pMI[12] <- 0
   pMI[13] <- mean.pMI
   pMI[14] <- 0
   pMI[15] <- 0
   
   mean.pMI ~ dunif(0, 1)
   
   pAR[1] <- 0
   pAR[2] <- mean.pAR
   pAR[3] <- 0
   pAR[4] <- 0
   pAR[5] <- mean.pAR
   pAR[6] <- 0
   pAR[7] <- 0
   pAR[8] <- mean.pAR
   pAR[9] <- 0
   pAR[10] <- 0
   pAR[11] <- mean.pAR
   pAR[12] <- 0
   pAR[13] <- 0
   pAR[14] <- mean.pAR
   pAR[15] <- 0
   
   mean.pAR ~ dunif(0, 1)

   # Transitions: gamma priors
   for (i in 1:5){
      a[i] ~ dgamma(1, 1)
      psiDB[i] <- a[i]/sum(a[])
      b[i] ~ dgamma(1, 1)
      psiJB[i] <- b[i]/sum(b[])
      c[i] ~ dgamma(1, 1)
      psiMI[i] <- c[i]/sum(c[])
      d[i] ~ dgamma(1, 1)
      psiAR[i] <- d[i]/sum(d[])
      
      for (t in (n.occasions-1)){
      e[t,i] ~ dgamma(1, 1)
      psiUN[i,t] <- e[i,t]/sum(e[])
      }
   }    

# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phiDB * psiDB[2]
      ps[1,i,t,3] <- phiDB * psiDB[3]
      ps[1,i,t,4] <- 0
      ps[1,i,t,5] <- phiDB * psiDB[5]
      ps[1,i,t,6] <- 1-phiDB
      
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 0
      ps[2,i,t,4] <- phiJB * psiJB[4]
      ps[2,i,t,5] <- phiJB * psiJB[5]
      ps[2,i,t,6] <- 1-phiJB
      
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phiMI * psiMI[4]
      ps[3,i,t,5] <- phiMI * psiMI[5]
      ps[3,i,t,6] <- 1-phiMI
      
      ps[4,i,t,1] <- phiAR * psiAR[1]
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- phiAR * psiAR[5]
      ps[4,i,t,6] <- 1-phiAR
      
      ps[5,i,t,1] <- phiUN * psiUN[t,1]
      ps[5,i,t,2] <- phiUN * psiUN[t,2]
      ps[5,i,t,3] <- phiUN * psiUN[t,3]
      ps[5,i,t,4] <- phiUN * psiUN[t,4]
      ps[5,i,t,5] <- phiUN * psiUN[t,5]
      ps[5,i,t,6] <- 1-phiUN     
      
      ps[6,i,t,1] <- 0
      ps[6,i,t,2] <- 0
      ps[6,i,t,3] <- 0
      ps[6,i,t,4] <- 0
      ps[6,i,t,5] <- 0
      ps[6,i,t,6] <- 1

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
      
      po[6,i,t,1] <- 0
      po[6,i,t,2] <- 0
      po[6,i,t,3] <- 0
      po[6,i,t,4] <- 0
      po[6,i,t,5] <- 1
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
