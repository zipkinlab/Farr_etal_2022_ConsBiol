#----------------#
#-Load Libraries-#
#----------------#

library(nimble)
library(coda)

#-----------#
#-Load Data-#
#-----------#

load(file = "~/HMSNO/DataFormat/HMSNO.data.Rdata")
load(file = "~/HMSNO/DataFormat/HMSNO.con.Rdata")

nspecs <- HMSNO.con$nspecs
parkS <- HMSNO.con$parkS
parkE <- HMSNO.con$parkE
nsite <- HMSNO.con$nsite
nstart <- HMSNO.con$nstart
nend <- HMSNO.con$nend
nparks <- max(parkE)
nyrs <- max(nend)

#--------------#
#-NIMBLE model-#
#--------------#

model.code <- nimbleCode({
  
  #-Priors-#
  
  # Detection hyperparameters
  mu.a0 ~ dunif(0, 1) # Community-level intercept on probability scale
  mu.a0L <- logit(mu.a0) # Community-level intercept on logit scale
  tau.a0 ~ dgamma(0.1, 0.1) # Community-level precision
  
  # Effect of effort (days sampled)
  alpha1 ~ dnorm(0, 0.1)
  
  # Initial abundance hyperparameters
  mu.b0 ~ dnorm(0, 0.1) # Community-level intercept on log scale
  tau.b0 ~ dgamma(0.1, 0.1) # Community-level precision
  
  # Apparent survival hyperparameters
  mu.o0 ~ dunif(0, 1) # Community-level intercept on probability scale
  mu.o0L <- logit(mu.o0) # Community-level intercept on logit scale
  tau.o0 ~ dgamma(0.1, 0.1) # Community-level precision
  
  # Gains hyperparameters
  mu.g0 ~ dnorm(0, 0.1) # Community-level intercept on log scale
  tau.g0 ~ dgamma(0.1, 0.1) # Community-level precision
  
  # Precison of park-level random effect on initial abundance
  tau.eps.l ~ dgamma(0.1, 0.1)
  
  # Precison of park-level random effect on apparent survival
  tau.eps.o ~ dgamma(0.1, 0.1)
  
  for(r in 1:nparks){
    
    # Park-level random effect on initial abundance
    eps.l[r] ~ dnorm(0, tau.eps.l)
    
    # Park-level random effect on apparent survival
    eps.o[r] ~ dnorm(0, tau.eps.o)
    
    # Predicted community-level apparent survival for each park 
    logit(park.surv[r]) <- mu.o0L + eps.o[r]
    
  } # End r
  
  for(i in 1:nspecs){
    
    # Species-specific intercept on detection
    alpha0[i] ~ dnorm(mu.a0L, tau.a0)
    
    # Species-specific intercept on initial abundance
    beta0[i] ~ dnorm(mu.b0, tau.b0)
    
    # Species-specific intercept on apparent survival
    omega0[i] ~ dnorm(mu.o0L, tau.o0)
    
    # Species-specific intercept on gains
    gamma0[i] ~ dnorm(mu.g0, tau.g0)
    
    for(r in parkS[i]:parkE[i]){
      
      # Predicted population-level apparent survival
      logit(pop.surv[r,i]) <- omega0[i] + eps.o[r]
      
      for(j in 1:nsite[r]){
        
        # Linear predictor for initial abundance
        log(lambda[j,r,i]) <- beta0[i] + eps.l[r]
        
        # Initial abundance
        N[nstart[r],j,r,i] ~ dpois(lambda[j,r,i])
        
        for(k in 1:nreps){
          
          # Linear predictor for detection probability (in year 1)
          logit(r[k,nstart[r],j,r,i]) <- alpha0[i] + alpha1 * days[k,nstart[r],j,i]
          
          # Observation process (in year 1)
          y[k,nstart[r],j,r,i] ~ dbern(p[k,nstart[r],j,r,i])
          
          # Site level detection (N-occupancy parameterization) (in year 1)
          p[k,nstart[r],j,r,i] <- 1 - pow((1 - r[k,nstart[r],j,r,i]), N[nstart[r],j,r,i])
          
        } # End k
        
        for(t in (nstart[r]+1):nend[r]){
          
          for(k in 1:nreps){
            
            # Linear predictor for detection probability (in year t+1)
            logit(r[k,t,j,r,i]) <- alpha0[i] + alpha1 * days[k,t,j,i]
            
            # Observation process (in year t+1)
            y[k,t,j,r,i] ~ dbern(p[k,t,j,r,i])
            
            # Observation process (in year t+1)
            p[k,t,j,r,i] <- 1 - pow((1 - r[k,t,j,r,i]), N[t,j,r,i])
            
          } # End k
          
          # Linear predictor for apparent survival
          logit(omega[t-1,j,r,i]) <- omega0[i] + eps.o[r]
          
          # Biological process
          N[t,j,r,i] <- S[t-1,j,r,i] + G[t-1,j,r,i]
          
          # Apparent survival
          S[t-1,j,r,i] ~ dbin(omega[t-1,j,r,i], N[t-1,j,r,i])
          
          # Gains
          G[t-1,j,r,i] ~ dpois(gamma[t-1,j,r,i])
          
          # Linear predictor for gains
          log(gamma[t-1,j,r,i]) <- gamma0[i]
          
        } # End t
        
      } # End j
      
      for(t in nstart[r]:nend[r]){
        
        # Population (species-park) abundance per year 
        Nhat[t,r,i] <- sum(N[t,1:nsite[r],r,i])
        
      } # End t
      
    } # End r
    
  } # End s
  
})

#--------------#
#-Compile data-#
#--------------#

Data <- list(y = HMSNO.data$y, days = HMSNO.data$days.scaled)

#----------------#
#-Initial values-#
#----------------#

dim2 <- dim1 <- dim(HMSNO.data$y)[-1]
dim2[1] <- dim2[1] - 1

Nst <- array(NA, dim = dim1)
Sst <- array(NA, dim = dim2)
Gst <- array(NA, dim = dim2)

for(i in 1:nspecs){
  for(r in parkS[i]:parkE[i]){
    for(j in 1:nsite[r]){
      Nst[nstart[r],j,r,i] <- 10
      Sst[nstart[r]:(nend[r]-1),j,r,i] <- 5
      Gst[nstart[r]:(nend[r]-1),j,r,i] <- 2
    }
  }
}

alpha0.fun <- function(){
  alpha0 <- NULL
  alpha0[1] <- runif(1, -1, 0)
  alpha0[2] <- runif(1, -1.5, -1)
  alpha0[3] <- runif(1, -1.5, -1)
  alpha0[4] <- runif(1, -2.5, -1.5)
  alpha0[5] <- runif(1, -1, -0.5)
  alpha0[6] <- runif(1, -1.5, -0.5)
  alpha0[7] <- runif(1, -0.7, -0.5) #Issues with this node
  alpha0[8] <- runif(1, -2, -1.5)
  alpha0[9] <- runif(1, -1, 0)
  alpha0[10] <- runif(1, -3, -2)
  alpha0[11] <- runif(1, -3, -2)
  alpha0[12] <- runif(1, -1.1, -0.9) #Issues with this node
  return(alpha0)
}

beta0.fun <- function(){
  beta0 <- NULL
  beta0[1] <- runif(1, 0, 1)
  beta0[2] <- runif(1, -1, 1)
  beta0[3] <- runif(1, 1, 3)
  beta0[4] <- runif(1, -8, -4)
  beta0[5] <- runif(1, 0, 2)
  beta0[6] <- runif(1, 0, 2)
  beta0[7] <- runif(1, -3, -1) #Issues with this node
  beta0[8] <- runif(1, 1, 2)
  beta0[9] <- runif(1, -3, -1)
  beta0[10] <- runif(1, 0, 2)
  beta0[11] <- runif(1, -5, -2)
  beta0[12] <- runif(1, -2, 1) #Issues with this node
  return(beta0)
}

omega0.fun <- function(){
  omega0 <- NULL
  omega0[1] <- runif(1, -0.5, 0.5)
  omega0[2] <- runif(1, -2, 0)
  omega0[3] <- runif(1, 0, 2)
  omega0[4] <- runif(1, -1, 1)
  omega0[5] <- runif(1, -1, 1)
  omega0[6] <- runif(1, 0, 2)
  omega0[7] <- runif(1, 1, 2) #Issues with this node
  omega0[8] <- runif(1, 0, 2)
  omega0[9] <- runif(1, -2, -1)
  omega0[10] <- runif(1, -1, 1)
  omega0[11] <- runif(1, -2, 2)
  omega0[12] <- runif(1, 1, 3)
  return(omega0)
}

gamma0.fun <- function(){
  gamma0 <- NULL
  gamma0[1] <- runif(1, 0, 0.5)
  gamma0[2] <- runif(1, -0.2, 0.2)
  gamma0[3] <- runif(1, -1.5, -0.5)
  gamma0[4] <- runif(1, -2.5, -2)
  gamma0[5] <- runif(1, -0.2, 0.2)
  gamma0[6] <- runif(1, -2.5, -1.5)
  gamma0[7] <- runif(1, -1.8, -1.4)
  gamma0[8] <- runif(1, -1, -0.5)
  gamma0[9] <- runif(1, -3, 2.5)
  gamma0[10] <- runif(1, -3, -1)
  gamma0[11] <- runif(1, -5, -3)
  gamma0[12] <- runif(1, -2.5, -1.5)
  return(gamma0)
}

eps.l.fun <- function(){
  eps.l <- NULL
  eps.l[1] <- runif(1, -2, 0)
  eps.l[2] <- runif(1, 2, 4)
  eps.l[3] <- runif(1, -2, 0)
  eps.l[4] <- runif(1, 0, 2)
  eps.l[5] <- runif(1, 0, 2)
  eps.l[6] <- runif(1, -2, 0)
  return(eps.l)
}

eps.o.fun <- function(){
  eps.o <- NULL
  eps.o[1] <- runif(1, -2, 0)
  eps.o[2] <- runif(1, 0, 2)
  eps.o[3] <- runif(1, -2, 0)
  eps.o[4] <- runif(1, -1, 1)
  eps.o[5] <- runif(1, -1, 1)
  eps.o[6] <- runif(1, -2, 0)
  return(eps.o)
}

inits <- function()list(N=Nst, S=Sst, G=Gst,
                        mu.a0 = runif(1, 0.15, 0.3), tau.a0 = runif(1, 1, 3),
                        mu.b0 = runif(1, -2, 1), tau.b0 = runif(1, 0.2, 0.4),
                        mu.o0 = runif(1, 0.4, 0.7), tau.o0 = runif(1, 0, 1.5), 
                        mu.g0 = runif(1, -2, -1), tau.g0 = runif(1, 0.25, 1), 
                        gamma0 = gamma0.fun(),
                        alpha0 = alpha0.fun(), alpha1 = runif(1, 0.2, 0.4),
                        beta0 = beta0.fun(), 
                        omega0 = omega0.fun(), 
                        eps.l = eps.l.fun(), 
                        eps.o = eps.o.fun(),
                        tau.eps.l = runif(1, 0, 1),
                        tau.eps.o = runif(1, 0, 1)
)

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c("mu.a0", "tau.a0", "alpha0", "alpha1",
            "mu.b0", "tau.b0", "beta0", "eps.l", "tau.eps.l",
            "mu.o0", "tau.o0", "omega0", "eps.o", "tau.eps.o",
            "mu.g0", "tau.g0", "gamma0",
            "park.surv", "pop.surv", "Npark") 

#---------------#
#-MCMC settings-#
#---------------#

model <- nimbleModel(model.code, constants = HMSNO.con, data = Data, inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)


MCMCconf$addSampler(target = c("beta0[3]", "beta0[6]", "beta0[9]", "beta0[10]", "eps.l[1]"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("beta0[1]", "beta0[2]", "beta0[4]", "beta0[5]", "beta0[11]", "eps.l[5]"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("beta0[5]", "beta0[8]", "beta0[12]", "eps.l[6]"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("omega0[3]", "omega0[6]", "omega0[9]", "omega0[10]", "eps.o[1]"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("omega0[1]", "omega0[2]", "omega0[4]", "omega0[5]", "omega0[11]", "eps.o[5]"),
                    type = "AF_slice")

MCMCconf$addSampler(target = c("omega0[5]", "omega0[8]", "omega0[12]", "eps.o[6]"),
                    type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

model.comp <- compileNimble(model, MCMC)

ni <- 100000
nb <- 75000
nc <- 1
nt <- 25

#Run NIMBLE model
out <- runMCMC(model.comp$MCMC, niter = ni, nburnin = nb, nchains = nc, thin = nt, samplesAsCodaMCMC = TRUE)

#-Save output-#
ID <- paste("chain", length(list.files(path = "~/HMSNO/DataAnalysis/CaseStudy", pattern = "chain", full.names = FALSE)) + 1, sep="")
assign(ID, out)
save(list = ID, file = paste0(ID, ".Rds"))