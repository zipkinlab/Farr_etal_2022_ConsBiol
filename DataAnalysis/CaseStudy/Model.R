#----------------#
#-Load Libraries-#
#----------------#

library(nimble)
library(coda)

#-----------#
#-Load Data-#
#-----------#

load(file = "../../DataFormat/HMSNO.Adata.Rdata")
load(file = "../../DataFormat/HMSNO.Acon.Rdata")

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

MSdyn.c <- nimbleCode({
  
  mu.b0 ~ dnorm(0, 0.1)
  tau.b0 ~ dgamma(0.1, 0.1)
  mu.o0L <- logit(mu.o0)
  mu.o0 ~ dunif(0, 1)
  tau.o0 ~ dgamma(0.1, 0.1)
  mu.g0 ~ dnorm(0, 0.1)
  tau.g0 ~ dgamma(0.1, 0.1)
  mu.a0 ~ dunif(0, 1)
  mu.a0L <- logit(mu.a0)
  tau.a0 ~ dgamma(0.1, 0.1)
  alpha1 ~ dnorm(0, 0.1)
  tau.eps.o ~ dgamma(0.1, 0.1)
  tau.eps.l ~ dgamma(0.1, 0.1)
  # tau.eps.g ~ dgamma(0.1, 0.1)
  
  for(i in 1:nparks){
    eps.o[i] ~ dnorm(0, tau.eps.o)
    eps.l[i] ~ dnorm(0, tau.eps.l)
    # eps.g[i] ~ dnorm(0, tau.eps.g)
    
    logit(park.surv[i]) <- mu.o0L + eps.o[i]
    # log(park.gain[i]) <- mu.g0 + eps.g[i]
  }
  
  for(s in 1:nspecs){
    
    alpha0[s] ~ dnorm(mu.a0L, tau.a0)
    beta0[s] ~ dnorm(mu.b0, tau.b0)
    omega0[s] ~ dnorm(mu.o0L, tau.o0)
    gamma0[s] ~ dnorm(mu.g0, tau.g0)
    
    for(i in parkS[s]:parkE[s]){
      
      logit(pop.surv[i,s]) <- omega0[s] + eps.o[i]
      # log(pop.gain[i,s]) <- gamma0[s] + eps.g[i]
      
      for(j in 1:nsite[i]){
        
        log(lambda[j,i,s]) <- beta0[s] + eps.l[i]
        
        N[nstart[i],j,i,s] ~ dpois(lambda[j,i,s])
        
        for(k in 1:nreps){
          
          logit(r[k,nstart[i],j,i,s]) <- alpha0[s] + alpha1 * days[k,nstart[i],j,i]
          
          y[k,nstart[i],j,i,s] ~ dbern(p[k,nstart[i],j,i,s])
          p[k,nstart[i],j,i,s] <- 1 - pow((1 - r[k,nstart[i],j,i,s]), N[nstart[i],j,i,s])
          
        }#end k
        
        for(t in (nstart[i]+1):nend[i]){
          
          for(k in 1:nreps){
            
            logit(r[k,t,j,i,s]) <- alpha0[s] + alpha1 * days[k,t,j,i]
            
            y[k,t,j,i,s] ~ dbern(p[k,t,j,i,s])
            p[k,t,j,i,s] <- 1 - pow((1 - r[k,t,j,i,s]), N[t,j,i,s])
            
          }#end k
          
          
          logit(omega[t-1,j,i,s]) <- omega0[s] + eps.o[i]
          
          N[t,j,i,s] <- S[t-1,j,i,s] + G[t-1,j,i,s]
          S[t-1,j,i,s] ~ dbin(omega[t-1,j,i,s], N[t-1,j,i,s])
          G[t-1,j,i,s] ~ dpois(gamma[t-1,j,i,s])
          
          log(gamma[t-1,j,i,s]) <- gamma0[s] #+ eps.g[i]
          
        }#end t
      }#end j
      
      for(t in nstart[i]:nend[i]){
        
        Npark[t,i,s] <- sum(N[t,1:nsite[i],i,s])
        
      }#end t
    }#end i
  }#end s
  
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

for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    for(j in 1:nsite[i]){
      Nst[nstart[i],j,i,s] <- 10
      Sst[nstart[i]:(nend[i]-1),j,i,s] <- 5
      Gst[nstart[i]:(nend[i]-1),j,i,s] <- 2
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

# eps.g.fun <- function(){
#   eps.g <- NULL
#   eps.g[1] <- runif(1, -1, 1)
#   eps.g[2] <- runif(1, -1, 1)
#   eps.g[3] <- runif(1, -1, 1)
#   eps.g[4] <- runif(1, -1, 1)
#   eps.g[5] <- runif(1, -1, 1)
#   eps.g[6] <- runif(1, -1, 1)
#   return(eps.g)
# }

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
                        # eps.g = eps.g.fun(),
                        tau.eps.l = runif(1, 0, 1),
                        tau.eps.o = runif(1, 0, 1)
                        # tau.eps.g = runif(1, 0, 1) 
)

#Parameters to save
params <- c("mu.a0",  
            "tau.a0", 
            "mu.b0",  
            "tau.b0", 
            "mu.o0",
            "tau.o0",
            "mu.g0",
            "tau.g0",
            "alpha0",
            "alpha1",
            "beta0",
            "omega0",
            "gamma0",
            "eps.l",
            "eps.o",
            # "eps.g",
            "tau.eps.l",
            "tau.eps.o",
            # "tau.eps.g",
            "park.surv",
            # "park.gain",
            "pop.surv",
            # "pop.gain",
            "Npark") 

#MCMC settings
MSdyn.m <- nimbleModel(MSdyn.c, constants = HMSNO.con, data = Data, inits = inits())

MCMCconf <- configureMCMC(MSdyn.m, monitors = params)


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

# MCMCconf$addSampler(target = c("gamma0[3]", "gamma0[6]", "gamma0[9]", "gamma0[10]", "eps.g[1]"),
#                     type = "AF_slice")
# 
# MCMCconf$addSampler(target = c("gamma0[1]", "gamma0[2]", "gamma0[4]", "gamma0[5]", "gamma0[11]", "eps.g[5]"),
#                     type = "AF_slice")
# 
# MCMCconf$addSampler(target = c("gamma0[5]", "gamma0[8]", "gamma0[12]", "eps.g[6]"),
#                     type = "AF_slice")


MCMC <- buildMCMC(MCMCconf)

MSdyn.comp <- compileNimble(MSdyn.m, MCMC)

ni <- 100000
nb <- 75000
nc <- 1
nt <- 25

#Run NIMBLE model
# invisible(capture.output(MSdyn.o <- runMCMC(MSdyn.comp$MCMC, niter = ni, nburnin = nb, nchains = nc, thin = nt, samplesAsCodaMCMC = TRUE)))
MSdyn.o <- runMCMC(MSdyn.comp$MCMC, niter = ni, nburnin = nb, nchains = nc, thin = nt, samplesAsCodaMCMC = TRUE)

#-Save output-#
ID <- paste("chain", length(list.files(pattern = "chain", full.names = FALSE)) + 1, sep="")
assign(ID, MSdyn.o)
save(list = ID, file = paste0(ID, ".Rds"))