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
   for (t in 1:(n.occasions-1)){
   phiDB[t] <- mean.phiDB
   phiJB[t] <- mean.phiJB
   phiMI[t] <- mean.phiMI
   phiAR[t] <- mean.phiAR
   phiUN[t] <- mean.phiUN
   
   pDB[t] <- mean.pDB
   pJB[t] ~ dunif(0, 1)
   pMI[t] <- mean.pDB
   pAR[t] <- mean.pAR
   
   

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probs
      for (i in 1:4){
         lpsiDB[t,i] ~ dnorm(0, 0.001)
         lpsiJB[t,i] ~ dnorm(0, 0.001)
         lpsiMI[t,i] ~ dnorm(0, 0.001)
         lpsiAR[t,i] ~ dnorm(0, 0.001)
         lpsiUN[t,i] ~ dnorm(0, 0.001)
      }
      
      # Constrain the transitions such that their sum is < 1
      for (i in 1:4){
         psiDB[t,i] <- exp(lpsiDB[t,i]) / (1 + exp(lpsiDB[t,1]) + exp(lpsiDB[t,2]) + exp(lpsiDB[t,3]) + exp(lpsiDB[t,4]))
         psiJB[t,i] <- exp(lpsiJB[t,i]) / (1 + exp(lpsiJB[t,1]) + exp(lpsiJB[t,2]) + exp(lpsiJB[t,3]) + exp(lpsiJB[t,4]))
         psiMI[t,i] <- exp(lpsiMI[t,i]) / (1 + exp(lpsiMI[t,1]) + exp(lpsiMI[t,2]) + exp(lpsiMI[t,3]) + exp(lpsiMI[t,4]))
         psiAR[t,i] <- exp(lpsiAR[t,i]) / (1 + exp(lpsiAR[t,1]) + exp(lpsiAR[t,2]) + exp(lpsiAR[t,3]) + exp(lpsiAR[t,4]))
         psiUN[t,i] <- exp(lpsiUN[t,i]) / (1 + exp(lpsiUN[t,1]) + exp(lpsiUN[t,2]) + exp(lpsiUN[t,3]) + exp(lpsiUN[t,4]))
      }
         
    # Calculate the last transition probability
      psiDB[t,5] <- 1-psiDB[t,1]-psiDB[t,2]-psiDB[t,3]-psiDB[t,4]
      psiJB[t,5] <- 1-psiJB[t,1]-psiJB[t,2]-psiJB[t,3]-psiJB[t,4]
      psiMI[t,5] <- 1-psiMI[t,1]-psiMI[t,2]-psiMI[t,3]-psiMI[t,4]
      psiAR[t,5] <- 1-psiAR[t,1]-psiAR[t,2]-psiAR[t,3]-psiAR[t,4]
      psiUN[t,5] <- 1-psiUN[t,1]-psiUN[t,2]-psiUN[t,3]-psiUN[t,4]
}
   mean.phiDB ~ dunif(0, 1)
   mean.phiJB ~ dunif(0, 1)
   mean.phiMI ~ dunif(0, 1)
   mean.phiAR ~ dunif(0, 1)
   mean.phiUN ~ dunif(0, 1)
   
   mean.pDB ~ dunif(0, 1)
   mean.pMI ~ dunif(0, 1)
   mean.pAR ~ dunif(0, 1)
   
# Define state-transition and observation matrices 	
for (i in 1:nind){
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- phiDB[t] * psiDB[t,2]
      ps[1,i,t,3] <- phiDB[t] * psiDB[t,3]
      ps[1,i,t,4] <- 0
      ps[1,i,t,5] <- phiDB[t] * psiDB[t,5]
      ps[1,i,t,6] <- 1-phiDB[t]
      
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 0
      ps[2,i,t,3] <- 0
      ps[2,i,t,4] <- phiJB[t] * psiJB[t,4]
      ps[2,i,t,5] <- phiJB[t] * psiJB[t,5]
      ps[2,i,t,6] <- 1-phiJB[t]
      
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- phiMI[t] * psiMI[t,4]
      ps[3,i,t,5] <- phiMI[t] * psiMI[t,5]
      ps[3,i,t,6] <- 1-phiMI[t]
      
      ps[4,i,t,1] <- phiAR[t] * psiAR[t,1]
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 0
      ps[4,i,t,5] <- phiAR[t] * psiAR[t,5]
      ps[4,i,t,6] <- 1-phiAR[t]
      
      ps[5,i,t,1] <- phiUN[t] * psiUN[t,1]
      ps[5,i,t,2] <- phiUN[t] * psiUN[t,2]
      ps[5,i,t,3] <- phiUN[t] * psiUN[t,3]
      ps[5,i,t,4] <- phiUN[t] * psiUN[t,4]
      ps[5,i,t,5] <- phiUN[t] * psiUN[t,5]
      ps[5,i,t,6] <- 1-phiUN[t]     
      
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
