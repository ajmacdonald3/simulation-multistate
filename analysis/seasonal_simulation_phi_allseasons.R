sink("simulation_scenarios.jags")
cat("
model {

# Priors and constraints
   # Survival and recapture: uniform
   
   # DB
   for (t in 1:(n.occasions-1)){
         logit(phiDB[t]) <- mu.phiDB[season[t]] + eps.phiDB[t]
         eps.phiDB[t] ~ dnorm(0, tau.phiDB[season[t]])T(-10,10)
   }
   
   for (t in 1:9){
         logit(pDB[t]) <- mu.pDB + eps.pDB[t]
         eps.pDB[t] ~ dnorm(0, tau.pDB)T(-10,10)
   }

      for (s in 1:4){
         mean.phiDB[s] ~ dunif(0, 1)
         mu.phiDB[s] <- logit(mean.phiDB[s])
         sigma.phiDB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiDB[s] <- pow(sigma.phiDB[s], -2)
         sigma2.phiDB[s] <- pow(sigma.phiDB[s], 2)
      }
      
         mean.pDB ~ dunif(0, 1)
         mu.pDB <- logit(mean.pDB)
         sigma.pDB ~ dunif(0, 10)               # Prior for standard deviation
         tau.pDB <- pow(sigma.pDB, -2)
         sigma2.pDB <- pow(sigma.pDB, 2)
         
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
   for (t in 1:(n.occasions-1)){
         logit(phiJB[t]) <- mu.phiJB[season[t]] + eps.phiJB[t]
         eps.phiJB[t] ~ dnorm(0, tau.phiJB[season[t]])T(-10,10)
   }
   
   for(t in 1:10){
         logit(pJB[t]) <- mu.pJB + eps.pJB[t]
         eps.pJB[t] ~ dnorm(0, tau.pJB)T(-10,10)
   }

      for (s in 1:4){
         mean.phiJB[s] ~ dunif(0, 1)
         mu.phiJB[s] <- logit(mean.phiJB[s])
         sigma.phiJB[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiJB[s] <- pow(sigma.phiJB[s], -2)
         sigma2.phiJB[s] <- pow(sigma.phiJB[s], 2)
      }
      
         mean.pJB ~ dunif(0, 1)
         mu.pJB <- logit(mean.pJB)
         sigma.pJB ~ dunif(0, 10)               # Prior for standard deviation
         tau.pJB <- pow(sigma.pJB, -2)
         sigma2.pJB <- pow(sigma.pJB, 2)
         
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
   for (t in 1:(n.occasions-1)){
         logit(phiMI[t]) <- mu.phiMI[season[t]] + eps.phiMI[t]
         eps.phiMI[t] ~ dnorm(0, tau.phiMI[season[t]])T(-10,10)
   }
   
   for (t in 1:10){
         logit(pMI[t]) <- mu.pMI + eps.pMI[t]
         eps.pMI[t] ~ dnorm(0, tau.pMI)T(-10,10)
   }

      for (s in 1:4){
         mean.phiMI[s] ~ dunif(0, 1)
         mu.phiMI[s] <- logit(mean.phiMI[s])
         sigma.phiMI[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiMI[s] <- pow(sigma.phiMI[s], -2)
         sigma2.phiMI[s] <- pow(sigma.phiMI[s], 2)
      }
      
         mean.pMI ~ dunif(0, 1)
         mu.pMI <- logit(mean.pMI)
         sigma.pMI ~ dunif(0, 10)               # Prior for standard deviation
         tau.pMI <- pow(sigma.pMI, -2)
         sigma2.pMI <- pow(sigma.pMI, 2)

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
   for (t in 1:(n.occasions-1)){
         logit(phiCC[t]) <- mu.phiCC[season[t]] + eps.phiCC[t]
         eps.phiCC[t] ~ dnorm(0, tau.phiCC[season[t]])T(-10,10)
   }
   
   for (t in 1:20){
         logit(pCC[t]) <- mu.pCC[seasonpCC[t]] + eps.pCC[t]
         eps.pCC[t] ~ dnorm(0, tau.pCC[seasonpCC[t]])T(-10,10)
   }

      for (s in 1:4){
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
         logit(phiSE[t]) <- mu.phiSE[season[t]] + eps.phiSE[t]
         eps.phiSE[t] ~ dnorm(0, tau.phiSE[season[t]])T(-10,10)
         
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
   for (t in 1:(n.occasions-1)){
         logit(phiBR[t]) <- mu.phiBR[season[t]] + eps.phiBR[t]
         eps.phiBR[t] ~ dnorm(0, tau.phiBR[season[t]])T(-10,10)
   }
   
   for (t in 1:20){
         logit(pBR[t]) <- mu.pBR[seasonpBR[t]] + eps.pBR[t]
         eps.pBR[t] ~ dnorm(0, tau.pBR[seasonpBR[t]])T(-10,10)
   }

      for (s in 1:4){
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
      for (t in 1:(n.occasions-1)){
         logit(phiAR[t]) <- mu.phiAR[season[t]] + eps.phiAR[t]
         eps.phiAR[t] ~ dnorm(0, tau.phiAR[season[t]])T(-10,10)
      }
      
      for (t in 1:10){
         logit(pAR[t]) <- mu.pAR + eps.pAR[t]
         eps.pAR[t] ~ dnorm(0, tau.pAR)T(-10,10)
   }

      for (s in 1:4){
         mean.phiAR[s] ~ dunif(0, 1)
         mu.phiAR[s] <- logit(mean.phiAR[s])
         sigma.phiAR[s] ~ dunif(0, 10)               # Prior for standard deviation
         tau.phiAR[s] <- pow(sigma.phiAR[s], -2)
         sigma2.phiAR[s] <- pow(sigma.phiAR[s], 2)
      }
      
         mean.pAR ~ dunif(0, 1)
         mu.pAR <- logit(mean.pAR)
         sigma.pAR ~ dunif(0, 10)               # Prior for standard deviation
         tau.pAR <- pow(sigma.pAR, -2)
         sigma2.pAR <- pow(sigma.pAR, 2)
         
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
