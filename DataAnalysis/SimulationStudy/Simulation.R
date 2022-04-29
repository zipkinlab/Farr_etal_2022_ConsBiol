#-----------#
#-Libraries-#
#-----------#

library(nimble)
library(coda)

#-----------#
#-Functions-#
#-----------#

#Logit function
logit <- function(pp) 
{ 
  log(pp) - log(1-pp) 
}

#Inverse logit
expit <- function(eta) 
{
  1/(1+exp(-eta))
}

#Species ID
spp.fun <- function(x,y){
  out <- c(
    which.max(apply(y, MARGIN = 1, sum)),
    which.min(apply(y, MARGIN = 1, sum)),
    which.min(apply(x, MARGIN = 1, sum)))
  return(out)
}


# "True" values
nYears <- 10
nReps <- 3
nSites <- 75
nSpecies <- 30

mu.lambda <- runif(1,0.1,1.5)
mu.lambdaL <- log(mu.lambda)
mu.omega <- runif(1,0,1)
mu.omegaL <- logit(mu.omega)
mu.gamma <- runif(1,0,1)
mu.gammaL <- log(mu.gamma)
mu.r <- runif(1,0,1)
mu.rL <- logit(mu.r)

sd.lambda <- 0.5
sd.omega <- 0.5
sd.gamma <- 0.5
sd.r <- 0.5

lambda0 <- rnorm(nSpecies, mu.lambdaL, sd.lambda)
omega0 <- rnorm(nSpecies, mu.omegaL, sd.omega)
gamma0 <- rnorm(nSpecies, mu.gammaL, sd.gamma)
r0 <- rnorm(nSpecies, mu.rL, sd.r)

lambda <- exp(lambda0)
# lambda <- lambda[order(lambda)]
omega <- expit(omega0)
# omega <- omega[order(omega)]
gamma <- exp(gamma0)
# gamma <- gamma[order(gamma)]
r <- expit(r0)
# r <- r[order(r)]

#Simulate true abundances, N, for each location
N <- array(NA, dim = c(nSpecies, nSites, nYears))
S <- G <- array(NA, dim = c(nSpecies, nSites, nYears-1))

#subsequent years follow the birth-death-immigration process
for(i in 1:nSpecies) {
  N[i,,1] <- rpois(nSites, lambda[i])
  for(t in 2:nYears) {
    S[i,,t-1] <- rbinom(nSites, N[i,,t-1], omega[i])
    G[i,,t-1] <- rpois(nSites, gamma[i])
    N[i,,t] <- S[i,,t-1] + G[i,,t-1] 
  }
}

# Generate data vector y for the counts
y <- array(NA, c(nSpecies, nSites, nYears, nReps))
for(i in 1:nSpecies){
  for(t in 1:nYears) {
    for(k in 1:nReps) {
      y[i,,t,k] <- rbinom(nSites, N[i,,t], r[i])
    }
  }
}

# Convert to occupancy data
y[which(y != 0)] <- 1


# Common, rare, elusive
spp.ID <- spp.fun(N,y)

#-------------------------#
#-Multispecies model code-#
#-------------------------#

code <- nimbleCode({
  
  mu.lambdaL ~ dnorm(0, 0.1)
  log(mu.lambda) <- mu.lambdaL
  tau.lambda ~ dgamma(0.1, 0.1)
  sd.lambda <- 1/sqrt(tau.lambda)
  mu.omegaL <- logit(mu.omega)
  mu.omega ~ dunif(0, 1)
  tau.omega ~ dgamma(0.1, 0.1)
  sd.omega <- 1/sqrt(tau.omega)
  mu.gammaL ~ dnorm(0, 0.1)
  log(mu.gamma) <- mu.gammaL
  tau.gamma ~ dgamma(0.1, 0.1)
  sd.gamma <- 1/sqrt(tau.gamma)
  mu.rL <- logit(mu.r)
  mu.r ~ dunif(0, 1)
  tau.r ~ dgamma(0.1, 0.1)
  sd.r <- 1/sqrt(tau.r)
  
  for(i in 1:nSpecies){
    lambda0[i] ~ dnorm(mu.lambdaL, tau.lambda)
    omega0[i] ~ dnorm(mu.omegaL, tau.omega)
    gamma0[i] ~ dnorm(mu.gammaL, tau.gamma)
    r0[i] ~ dnorm(mu.r, tau.r)
    
    log(lambda[i]) <- lambda0[i]
    logit(omega[i]) <- omega0[i]
    log(gamma[i]) <- gamma0[i]
    logit(r[i]) <- r0[i]
    
    for(j in 1:nSites){
      
      N[i,j,1] ~ dpois(lambda[i])
      
      for(t in 2:nYears){
        
        S[i,j,t-1] ~ dbin(omega[i], N[i,j,t-1])
        G[i,j,t-1] ~ dpois(gamma[i])
        N[i,j,t] <- S[i,j,t-1] + G[i,j,t-1]
        
      }#end t
      
      for(t in 1:nYears){
        
        p[i,j,t] <- 1 - pow((1 - r[i]), N[i,j,t])
        
        for(k in 1:nReps){
          
          y[i,j,t,k] ~ dbern(p[i,j,t])
          
        }#end k
      }#end t
    }#end j
  }#end i
})

#--------------#
#-Compile data-#
#--------------#

data <- list(y = y)

constants <- list(nSpecies = nSpecies, nSites = nSites, nYears = nYears, nReps = nReps)

#----------------#
#-Initial values-#
#----------------#

inits <- function()list(N = N, S = S, G = G,
                        mu.lambdaL = mu.lambdaL,
                        # tau.lambda,
                        mu.omega = mu.omega,
                        # tau.omega,
                        mu.gammaL = mu.gammaL,
                        # tau.gamma,
                        mu.r = mu.r,
                        # tau.r,
                        lambda0 = lambda0,
                        omega0 = omega0,
                        gamma0 = gamma0,
                        r0 = r0
                        )

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c("mu.lambda", "sd.lambda",
            "mu.omega", "sd.omega",
            "mu.gamma", "sd.gamma",
            "mu.r", "sd.r",
            "lambda", "omega", "gamma", "r")

#---------------#
#-MCMC settings-#
#---------------#

model <- nimbleModel(code = code, 
                     constants = constants, 
                     data = data, 
                     inits =  inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

ni <- 35000
nb <- 10000
nc <- 3
nt <- 25

# ni <- 350
# nb <- 100
# nc <- 3
# nt <- 1

#-----------#
#-Run model-#
#-----------#
  
out <- runMCMC(compiled.model$MCMC, 
               niter = ni, nburnin = nb, 
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#----------------#
#-Compile output-#
#----------------#

coefs <- c(outer(paste0(c("gamma", "lambda", "omega", "r"), "["), paste0(spp.ID, "]"), paste0),
            "mu.gamma", "mu.lambda", "mu.omega", "mu.r", "sd.gamma", "sd.lambda", "sd.omega", "sd.r")

coefs.names <- c(outer(c("gamma.", "lambda.", "omega.", "r."), c("common", "elusive", "rare"), paste0),
                 "mu.gamma", "mu.lambda", "mu.omega", "mu.r", "sd.gamma", "sd.lambda", "sd.omega", "sd.r")

output <- cbind(sapply(coefs, FUN = function(x){eval(parse(text = x))}),
                summary(out)$statistics[coefs,"Mean"],
                summary(out)$statistics[coefs,"SD"],
                gelman.diag(out)$psrf[coefs,1])

# output <- cbind(c(gamma, lambda, mu.gamma, mu.lambda, 
#                   mu.omega, mu.r, omega, r,
#                   sd.gamma, sd.lambda, sd.omega, sd.r),
#                 summary(out)$statistics[,"Mean"],
#                 gelman.diag(out)$psrf[,1])

output <- as.data.frame(output)
# output$Params <- rownames(output)
output$Params <- coefs.names
rownames(output) <- NULL
colnames(output)[1:4] <- c("Truth", "Mean_MS", "SD_MS", "Rhat_MS")
output <- output[,c(5,1,2,3,4)]
output$Rhat_SS <- output$SD_SS <- output$Mean_SS <- NA

#---------------------------#
#-Single-species model code-#
#---------------------------#

code <- nimbleCode({
  
  lambda0 ~ dnorm(0, 0.01)
  omega ~ dunif(0, 1)
  gamma0 ~ dnorm(0, 0.01)
  r ~ dunif(0, 1)
  
  log(lambda) <- lambda0
  log(gamma) <- gamma0
  
  for(j in 1:nSites){
    
    N[j,1] ~ dpois(lambda)
    
    for(t in 2:nYears){
      
      S[j,t-1] ~ dbin(omega, N[j,t-1])
      G[j,t-1] ~ dpois(gamma)
      N[j,t] <- S[j,t-1] + G[j,t-1]
      
    }#end t
    
    for(t in 1:nYears){
      
      p[j,t] <- 1 - pow((1 - r), N[j,t])
      
      for(k in 1:nReps){
        
        y[j,t,k] ~ dbern(p[j,t])
        
      }#end k
    }#end t
  }#end j
})

#--------------#
#-Compile data-#
#--------------#

iter <- 1

# for(i in 1:nSpecies){
for(i in unique(spp.ID)){

data <- list(y = y[i,,,])

constants <- list(nSites = nSites, nYears = nYears, nReps = nReps)

#----------------#
#-Initial values-#
#----------------#

inits <- function()list(N = N[i,,], S = S[i,,], G = G[i,,],
                        lambda0 = lambda0[i],
                        omega = omega[i],
                        gamma0 = gamma0[i],
                        r = r[i]
)

#--------------------#
#-Parameters to save-#
#--------------------#

params <- c("gamma", "lambda", "omega", "r")

#---------------#
#-MCMC settings-#
#---------------#

model <- nimbleModel(code = code, 
                     constants = constants, 
                     data = data, 
                     inits =  inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

#-----------#
#-Run model-#
#-----------#

out <- runMCMC(compiled.model$MCMC, 
               niter = ni, nburnin = nb, 
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#----------------#
#-Compile output-#
#----------------#


if(iter == 1){
  
  output[1:4,6] <- summary(out)$statistics[,"Mean"]
  output[1:4,7] <- summary(out)$statistics[,"SD"]
  output[1:4,8] <- gelman.diag(out)$psrf[,1]
  
}else{
  
  if(length(unique(spp.ID) == 2)){
    
    output[9:12,6] <- output[5:8,6] <- summary(out)$statistics[,"Mean"]
    output[9:12,7] <- output[5:8,7] <- summary(out)$statistics[,"SD"]
    output[9:12,8] <- output[5:8,8] <- gelman.diag(out)$psrf[,1]
    
  }else{
    
    if(iter == 2){
      
      output[5:8,6] <- summary(out)$statistics[,"Mean"]
      output[5:8,7] <- summary(out)$statistics[,"SD"]      
      output[5:8,8] <- gelman.diag(out)$psrf[,1]
      
    }else{
      
      output[9:12,6] <- summary(out)$statistics[,"Mean"]
      output[9:12,7] <- summary(out)$statistics[,"SD"]
      output[9:12,8] <- gelman.diag(out)$psrf[,1]
      
    }
  }
}

# output[c(i,i+nSpecies,i+4+nSpecies*2,i+4+nSpecies*3),5] <- summary(out)$statistics[,"Mean"]
# 
# output[c(i,i+nSpecies,i+4+nSpecies*2,i+4+nSpecies*3),6] <- gelman.diag(out)$psrf[,1]
# 
# print(i)

iter <- iter + 1

print(iter)

}#end i

#-------------#
#-Save output-#
#-------------#

ID <- length(list.files("./output_75_final/")) + 1
save(output, file = paste("./output_75_final/output", ID, ".Rds", sep=""))
