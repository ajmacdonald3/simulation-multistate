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

   # Transitions: multinomial logit
      # Normal priors on logit of all but one transition probs
      for (i in 1:4){
         lpsiDB[i] ~ dnorm(0, 0.001)
         lpsiJB[i] ~ dnorm(0, 0.001)
         lpsiMI[i] ~ dnorm(0, 0.001)
         lpsiAR[i] ~ dnorm(0, 0.001)
      }
      
      for (i in 1:4){
      for (t in (n.occasions-1)){
         lpsiUN[i,t] ~ dnorm(0, 0.001)
        }
      }
         
      # Constrain the transitions such that their sum is < 1
      for (i in 1:4){
         psiDB[i] <- exp(lpsiDB[i]) / (1 + exp(lpsiDB[1]) + exp(lpsiDB[2]) + exp(lpsiDB[3]) + exp(lpsiDB[4]))
         psiJB[i] <- exp(lpsiJB[i]) / (1 + exp(lpsiJB[1]) + exp(lpsiJB[2]) + exp(lpsiJB[3]) + exp(lpsiJB[4]))
         psiMI[i] <- exp(lpsiMI[i]) / (1 + exp(lpsiMI[1]) + exp(lpsiMI[2]) + exp(lpsiMI[3]) + exp(lpsiMI[4]))
         psiAR[i] <- exp(lpsiAR[i]) / (1 + exp(lpsiAR[1]) + exp(lpsiAR[2]) + exp(lpsiAR[3]) + exp(lpsiAR[4]))
      }
       
      for (i in 1:4){
      for (t in (n.occasions-1)){
         psiUN[i,t] <- exp(lpsiUN[i,t]) / (1 + exp(lpsiUN[1,t]) + exp(lpsiUN[2,t]) + exp(lpsiUN[3,t]) + exp(lpsiUN[4,t]))
        }
      }
         
    # Calculate the last transition probability
      psiDB[5] <- 1-psiDB[1]-psiDB[2]-psiDB[3]-psiDB[4]
      psiJB[5] <- 1-psiJB[1]-psiJB[2]-psiJB[3]-psiJB[4]
      psiMI[5] <- 1-psiMI[1]-psiMI[2]-psiMI[3]-psiMI[4]
      psiAR[5] <- 1-psiAR[1]-psiAR[2]-psiAR[3]-psiAR[4]
      
    for (t in (n.occasions-1)){
      psiUN[5,t] <- 1-psiUN[1,t]-psiUN[2,t]-psiUN[3,t]-psiUN[4,t]
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
