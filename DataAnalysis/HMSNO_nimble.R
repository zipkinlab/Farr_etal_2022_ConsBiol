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

#---------#
#-Indices-#
#---------#

nspec <- HMSNO.con$nspec
nparks <- max(HMSNO.data$park)
npark <- HMSNO.con$npark
nsites <- max(HMSNO.con$site)
nsite <- HMSNO.con$nsite
nyrs <- max(HMSNO.con$yr)
nstart <- HMSNO.data$nstart
nend <- HMSNO.data$nend
parkspec <- HMSNO.con$parkspec
sitespec <- HMSNO.con$sitespec

#--------------#
#-NIMBLE model-#
#--------------#

MSdyn.c <- nimbleCode({
  
  mu.b0 ~ dnorm(0, 0.1)
  tau.b0 ~ dgamma(0.1, 0.1)
  mu.o0L <- logit(mu.o0)
  mu.o0 ~ dunif(0, 1)
  tau.o0 ~ dgamma(0.1, 0.1)
  mu.o1 ~ dnorm(0, 0.1)
  tau.o1 ~ dgamma(0.1, 0.1)
  mu.o2 ~ dnorm(0, 0.1)
  tau.o2 ~ dgamma(0.1, 0.1)
  
  log(gamma) <- gamma0
  gamma0 ~ dnorm(0, 0.1)
  
  alpha0 ~ dunif(0, 1)
  alpha1 ~ dnorm(0, 0.1)
  
  for(s in 1:nspec){
    
    omega0[s] ~ dnorm(mu.o0L, tau.o0)
    omega1[s] ~ dnorm(mu.o1, tau.o1)
    omega2[s] ~ dnorm(mu.o2, tau.o2)
    
    for(j in sitespec[1,s]:sitespec[nsite[s],s]){
      logit(r[ns[j],j,s]) <- logit(alpha0) + alpha1 * days[ns[j],j]
      N[ns[j],j,s] ~ dpois(lambda[park[j],s])
      for(t in (ns[j]+1):ne[j]){
        
        logit(r[t,j,s]) <- alpha0 + alpha1 * days[t,j]
        
        logit(omega[t-1,j,s]) <- omega0[s] + omega1[s] * density[j] + omega2[s] * edge[t,j]
        
        N[t,j,s] <- S[t-1,j,s] + G[t-1,j,s]
        S[t-1,j,s] ~ dbin(omega[t-1,j,s], N[t-1,j,s])
        #G[t-1,j,s] ~ dpois(gamma[t-1,s])
        G[t-1,j,s] ~ dpois(gamma)
        
      }#end t
    }#end j
    
    for(i in 1:npark[s]){
      beta0[parkspec[i,s],s] ~ dnorm(mu.b0, tau.b0)
      log(lambda[parkspec[i,s],s]) <- beta0[parkspec[i,s],s]
      for(t in nyrstr[parkspec[i,s]]:nyrend[parkspec[i,s]]){
        Npark[t,parkspec[i,s],s] <- sum(N[t,nsitestr[parkspec[i,s]]:nsiteend[t,parkspec[i,s]],s])
      }#end t
    }#end i
    
  }#end s
  
  
  for(k in 1:nobs){
    y[k] ~ dbern(p[k])
    p[k] <- 1 - pow((1 - r[yr[k], site[k], spec[k]]), N[yr[k], site[k], spec[k]])
  }#end k
  
})

#---------------#
#-Inital values-#
#---------------#

Nst <- array(NA, dim = c(nyrs, nsites, nspec))
Sst <- array(NA, dim = c(nyrs-1, nsites, nspec))
Gst <- array(NA, dim = c(nyrs-1, nsites, nspec))
for(s in 1:nspec){
  for(j in 1:nsite[s]){
    jj <- sitespec[j,s]
    Nst[nstart[jj],jj,s] <- 10
    Sst[nstart[jj]:(nend[jj]-1),jj,s] <- 5
    Gst[nstart[jj]:(nend[jj]-1),jj,s] <- 2
  }
}
# for(j in 1:nsites){
#   Nst[nstart[j],j,] <- 10
#   Sst[nstart[j]:(nend[j]-1),j,] <- 5
#   Gst[nstart[j]:(nend[j]-1),j,] <- 2
# }

# alpha0.fun <- function(){
#   alpha0 <- NULL
#   alpha0[1] <- runif(1, -1.5, -1)
#   alpha0[2] <- runif(1, -1.5, -1)
#   alpha0[3] <- runif(1, -2, -1.5)
#   alpha0[4] <- runif(1, -2.5, -1.5)
#   alpha0[5] <- runif(1, -1.5, -1)
#   alpha0[6] <- runif(1, -1.5, -1)
#   alpha0[7] <- runif(1, -2, -1.5)
#   alpha0[8] <- runif(1, -4, -2)
#   alpha0[9] <- runif(1, -2, -1.5)
#   return(alpha0)
# }

# beta0.fun <- function(){
#   beta0 <- NULL
#   beta0[1] <- runif(1, 1.5, 2.5)
#   beta0[2] <- runif(1, 0.5, 1.25)
#   beta0[3] <- runif(1, 0, 0.75)
#   beta0[4] <- runif(1, -5, -3)
#   beta0[5] <- runif(1, 1.5, 2.25)
#   beta0[6] <- runif(1, -0.25, 0.5)
#   beta0[7] <- runif(1, 0.75, 1.5)
#   beta0[8] <- runif(1, -0.75, 0)
#   beta0[9] <- runif(1, -0.25, 0.5)
#   return(beta0)
# }

beta0.fun <- function(){
  beta0 <- matrix(NA, nrow = nparks, ncol = nspec)
  # for(s in 1:nspec){
  #   for(i in specPark[1:specNpark[s],s]){
  #     beta0[i,s] <- runif(1, -1, 1)
  #   }
  # }
  for(s in 1:nspec){
    for(i in 1:npark[s]){
      ii <- parkspec[i,s]
      beta0[ii,s] <- runif(1, -1, 2)
    }
  }
  return(beta0)
}

# beta4.fun <- function(){
#   beta4 <- NULL
#   beta4[1] <- runif(1, 0, 1)
#   beta4[2] <- runif(1, 0, 2)
#   beta4[3] <- runif(1, -1, 0)
#   beta4[4] <- runif(1, 0, 4)
#   beta4[5] <- runif(1, 2, 3)
#   beta4[6] <- runif(1, 1, 1.5)
#   beta4[7] <- runif(1, 0, 2)
#   beta4[8] <- runif(1, 1, 2)
#   beta4[9] <- runif(1, 0.4, 0.6)
#   return(beta4)
# }

omega0.fun <- function(){
  omega0 <- NULL
  omega0[1] <- runif(1, 0.2, 0.4)
  omega0[2] <- runif(1, 0, 0.25)
  omega0[3] <- runif(1, 1, 1.5)
  omega0[4] <- runif(1, -1.5, -0.5)
  omega0[5] <- runif(1, 0.4, 0.6)
  omega0[6] <- runif(1, 1.9, 2.2)
  omega0[7] <- runif(1, -0.25, 0)
  omega0[8] <- runif(1, -0.5, 0)
  omega0[9] <- runif(1, 2, 2.25)
  return(omega0)
}

# gamma0.fun <- function(){
#   gamma0 <- NULL
#   gamma0[1] <- runif(1, 0, 1)
#   gamma0[2] <- runif(1, -0.2, 0.2)
#   gamma0[3] <- runif(1, -1.5, -0.5)
#   gamma0[4] <- runif(1, -2.5, -1.5)
#   gamma0[5] <- runif(1, -0.6, -0.2)
#   gamma0[6] <- runif(1, -3, -2)
#   gamma0[7] <- runif(1, -1.5, -0.5)
#   gamma0[8] <- runif(1, -3, -1)
#   gamma0[9] <- runif(1, -3, -2)
#   return(gamma0)
# }

# inits <- function()list(N=Nst, S=Sst, G=Gst,
#                         mu.a0 = runif(1, 0.1, 0.25), tau.a0 = runif(1, 0, 5), alpha0 = alpha0.fun(), alpha1 = runif(1, -0.25, 0),
#                         mu.b0 = runif(1, 0, 3), tau.b0 = runif(1, 0, 1), beta0 = beta0.fun(),
#                         mu.b1 = runif(1, -1, 1), tau.b1 = runif(1, 0, 10), beta1 = runif(nspec, -1, 1),
#                         mu.b2 = runif(1, -1, 1), tau.b2 = runif(1, 0, 10), beta2 = runif(nspec, -0.5, 0.5),
#                         mu.b3 = runif(1, -1, 1), tau.b3 = runif(1, 0, 10), beta3 = runif(nspec, -1, 1),
#                         mu.b4 = runif(1, 0, 2), tau.b4 = runif(1, 0, 2), beta4 = beta4.fun(),
#                         mu.o0 = runif(1, 0.5, 0.9), tau.o0 = runif(1, 0, 1), omega0 = omega0.fun(),
#                         mu.o1 = runif(1, 0, 4), tau.o1 = runif(1, 0, 1), omega1 = runif(nspec, -5, 5),
#                         mu.o2 = runif(1, -1, 1), tau.o2 = runif(1, 0, 10), omega2 = runif(nspec, -1, 1),
#                         mu.o3 = runif(1, -1, 1), tau.o3 = runif(1, 0, 10), omega3 = runif(nspec, -1, 1),
#                         mu.g0 = runif(1, -2, 0), tau.g0 = runif(1, 0, 1), gamma0 = gamma0.fun()
#                         )

inits <- function()list(N=Nst, S=Sst, G=Gst,
                        alpha0 = runif(1, 0.15, 0.2), alpha1 = runif(1, -0.5, -0.45),
                        mu.b0 = runif(1, -0.3, 0.75), #sig.b0 = runif(1, 1, 2), beta0 = beta0.fun(),
                        beta0 = beta0.fun(), tau.b0 = runif(1, 0, 1),
                        mu.o0 = runif(1, 0.5, 0.75), tau.o0 = runif(1, 0.5, 1.25), omega0 = omega0.fun(),
                        mu.o1 = runif(1, 0, 4), tau.o1 = runif(1, 0, 1), omega1 = runif(nspec, -5, 5),
                        mu.o2 = runif(1, -1, 1), tau.o2 = runif(1, 0, 10), omega2 = runif(nspec, -1, 1),
                        gamma0 = runif(1, -2, -1.75) #, tau.p = runif(1, 0.5, 2)
)

#Parameters to save
# params <- c("mu.a0", "tau.a0", 
#             "mu.b0", "tau.b0", "mu.b1", "tau.b1", "mu.b2", "tau.b2", "mu.b3", "tau.b3", "mu.b4", "tau.b4",
#             "mu.o0", "tau.o0", "mu.o1", "tau.o1", "mu.o2", "tau.o2", "mu.o3", "tau.o3",
#             "mu.g0", "tau.g0",
#             "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "beta4",
#             "omega0", "omega1", "omega2", "omega3", "gamma0")

params <- c("mu.b0", "tau.b0", "mu.o0", "tau.o0", 
            "mu.o1", "tau.o1", "mu.o2", "tau.o2",
            "alpha0", "alpha1", "beta0",
            "omega0", "omega1", "omega2",
            "gamma0", "Npark")

#MCMC settings
MSdyn.m <- nimbleModel(MSdyn.c, constants = HMSNO.con, data = HMSNO.data, inits = inits())

MCMCconf <- configureMCMC(MSdyn.m, monitors = params)
#MCMCconf$printSamplers(1:122)

# MCMCconf$removeSampler(c('alpha0', 'alpha1', 'gamma0',
#                          'beta0', 'beta1', 'beta2', 'beta3', 'beta4',
#                          'omega0', 'omega1', 'omega2', 'omega3'))
# 
# MCMCconf$addSampler(target = c('mu.a0', 'alpha0', 'alpha1'),
#                     type = "AF_slice")
# 
# MCMCconf$addSampler(target = c('mu.b0', 'beta0'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.b1', 'beta1'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.b2', 'beta2'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.b3', 'beta3'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.b4', 'beta4'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.o0', 'omega0'),
#                     type = "AF_slice")
# 
# MCMCconf$addSampler(target = c('mu.o1', 'omega1'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.o2', 'omega2'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.o3', 'omega3'),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c('mu.g0', 'gamma0'),
#                     type = "RW_block")
# 
# nFun <- function(node, i){
#   paste(node, "[", i, "]", sep = "")
# }
# 
# for(i in 1:nspec){
#   MCMCconf$addSampler(target = c(nFun("beta0", i), nFun("beta1", i), 
#                                  nFun("beta2", i), nFun("beta3", i),
#                                  nFun("beta4", i)),
#                       type = "RW_block")
#   
#   MCMCconf$addSampler(target = c(nFun("omega0", i), nFun("omega1", i), 
#                                  nFun("omega2", i), nFun("omega3", i)),
#                       type = "AF_slice")  
# }
# 
# MCMCconf$samplerExecutionOrder <- c(seq(1,22,1), seq(52, 41641, 1), seq(23, 51, 1))

MCMCconf$addSampler(target = c('alpha0', 'alpha1'), type = "AF_slice")

species <- seq(1, nspec, 1)

nFun <- function(node, i){
  paste(node, "[", i, "]", sep = "")
}

omit <- NULL

# MCMCconf$addSampler(target = c("beta0[1]", "beta0[2]", "beta0[3]", "beta0[4]", "beta0[5]",
#                                "beta0[6]", "beta0[7]", "beta0[8]", "beta0[9]"),
#                     type = "RW_block")

MCMCconf$addSampler(target = c("omega0[1]", "omega0[2]", "omega0[3]", "omega0[4]", "omega0[5]",
                               "omega0[6]", "omega0[7]", "omega0[8]", "omega0[9]"),
                    type = "RW_block")

for(i in 1:nspec){
  omit <- c(omit, i)
  if(i != nspec){
    for(j in species[-omit]){
      # MCMCconf$addSampler(target = c(nFun("beta0", i), nFun("beta0", j)),
      #                     type = "RW_block") 
      MCMCconf$addSampler(target = c(nFun("omega0", i), nFun("omega0", j)),
                          type = "RW_block")
    } 
  }
  # MCMCconf$addSampler(target = c("mu.b0", nFun("beta0", i)),
  #                     type = "RW_block") 
  MCMCconf$addSampler(target = c("mu.o0", nFun("omega0", i)),
                      type = "AF_slice")  
  MCMCconf$addSampler(target = c(nFun("omega0", i), nFun("omega1", i),
                                 nFun("omega2", i)),
                      type = "AF_slice")
}


#MCMCconf$addSampler(target = c('mu.b0', 'beta0'), type = "RW_block")

#MCMCconf$addSampler(target = c('mu.o0', 'omega0'), type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

MSdyn.comp <- compileNimble(MSdyn.m, MCMC)

ni <- 10000
nb <- 5000
nc <- 3
nt <- 5

#Run NIMBLE model
MSdyn.o <- runMCMC(MSdyn.comp$MCMC, niter = ni, nburnin = nb, nchains = nc, thin = nt, samplesAsCodaMCMC = TRUE)

#-Save output-#
ID <- gsub(" ","_",Sys.time())
ID <- gsub(":", "-", ID)
save(MSdyn.o, file = paste("output", ID, ".Rdata", sep=""))

#traceplot(out[c(1:3)][,"Npark[1, 1, 1]"])
#plot(unlist(output[c(1:5)][,"beta0[1]"]), unlist(output[c(1:5)][,"beta0[4]"])

for(i in 1:npark){
  plot(nyrstr[i], mean(unlist(output[c(1:5)][,paste("Npark[", nyrstr[i], ", ", i, ", 8]", sep = "")])),
       xlim = c(1,9), ylim = c(0,500), xlab = "Abundance", ylab = "Year")
  for(t in (nyrstr[i] + 1):nyrend[i]){
    points(t, mean(unlist(output[c(1:5)][,paste("Npark[", t, ", ", i, ", 8]", sep = "")])))
  }
}
