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

HMSNO.con$nyrs <- nyrs

#--------------#
#-NIMBLE model-#
#--------------#

MSdyn.c <- nimbleCode({
  
  mu.b0 ~ dnorm(0, 0.1)
  tau.b0 ~ dgamma(0.1, 0.1)
  mu.a0L <- logit(mu.a0)
  mu.a0 ~ dunif(0, 1)
  tau.a0 ~ dgamma(0.1, 0.1)
  mu.o0L <- logit(mu.o0)
  mu.o0 ~ dunif(0, 1)
  tau.o0 ~ dgamma(0.1, 0.1)
  mu.o1 ~ dnorm(0, 0.1)
  tau.o1 ~ dgamma(0.1, 0.1)
  mu.g0 ~ dnorm(0, 0.1)
  tau.g0 ~ dgamma(0.1, 0.1)
  tau.phi ~ dgamma(0.1, 0.1)
  alpha1 ~ dnorm(0, 0.1)
  
  for(t in 1:nyrs){
    eps.o[t] ~ T(dnorm(0, tau.phi), -5, 5)
  }#end t
  
  for(s in 1:nspec){
    
    omega1[s] ~ dnorm(mu.o1, tau.o1)
    
    for(j in sitespec[1,s]:sitespec[nsite[s],s]){
      logit(r[ns[j],j,s]) <- alpha0[park[j],s] + alpha1 * days[ns[j],j]
      N[ns[j],j,s] ~ dpois(lambda[park[j],s])
      for(t in (ns[j]+1):ne[j]){
        
        logit(r[t,j,s]) <- alpha0[park[j],s] + alpha1 * days[t,j]
        logit(omega[t-1,j,s]) <- omega0[park[j],s] + omega1[s] * edge[t-1,j] + eps.o[t]
        
        N[t,j,s] <- S[t-1,j,s] + G[t-1,j,s]
        S[t-1,j,s] ~ dbin(omega[t-1,j,s], N[t-1,j,s])
        G[t-1,j,s] ~ dpois(gamma[park[j],s])
        
      }#end t
    }#end j
    
    for(i in 1:npark[s]){
      alpha0[parkspec[i,s],s] ~ dnorm(mu.a0L, tau.a0)
      beta0[parkspec[i,s],s] ~ dnorm(mu.b0, tau.b0)
      omega0[parkspec[i,s],s] ~ dnorm(mu.o0, tau.o0)
      gamma0[parkspec[i,s],s] ~ dnorm(mu.g0, tau.g0)
      
      log(lambda[parkspec[i,s],s]) <- beta0[parkspec[i,s],s]
      log(gamma[parkspec[i,s],s]) <- gamma0[parkspec[i,s],s]
      
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

alpha0.fun <- function(){
  alpha0 <- matrix(NA, nrow = nparks, ncol = nspec)
  for(s in 1:nspec){
    for(i in 1:npark[s]){
      ii <- parkspec[i,s]
      alpha0[ii,s] <- runif(1, -3, -1)
    }
  }
  return(alpha0)
}

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

# omega0.fun <- function(){
#   omega0 <- NULL
#   omega0[1] <- runif(1, 0.2, 0.4)
#   omega0[2] <- runif(1, 0, 0.25)
#   omega0[3] <- runif(1, 1, 1.5)
#   omega0[4] <- runif(1, -1.5, -0.5)
#   omega0[5] <- runif(1, 0.4, 0.6)
#   omega0[6] <- runif(1, 1.9, 2.2)
#   omega0[7] <- runif(1, -0.25, 0)
#   omega0[8] <- runif(1, -0.5, 0)
#   omega0[9] <- runif(1, 2, 2.25)
#   return(omega0)
# }

omega0.fun <- function(){
  omega0 <- matrix(NA, nrow = nparks, ncol = nspec)
  for(s in 1:nspec){
    for(i in 1:npark[s]){
      ii <- parkspec[i,s]
      omega0[ii,s] <- runif(1, -1, 2)
    }
  }
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

gamma0.fun <- function(){
  gamma0 <- matrix(NA, nrow = nparks, ncol = nspec)
  for(s in 1:nspec){
    for(i in 1:npark[s]){
      ii <- parkspec[i,s]
      gamma0[ii,s] <- runif(1, -1, 2)
    }
  }
  return(gamma0)
}

inits <- function()list(N=Nst, S=Sst, G=Gst,
                        mu.a0 = runif(1, 0.1, 0.25), tau.a0 = runif(1, 0, 5),
                        mu.b0 = runif(1, -0.3, 0.75), tau.b0 = runif(1, 0, 1),
                        mu.o0 = runif(1, 0.5, 0.75), tau.o0 = runif(1, 0.5, 1.25),
                        mu.o1 = runif(1, 0, 4), tau.o1 = runif(1, 0, 1),
                        mu.g0 = runif(1, -2, 0), tau.g0 = runif(1, 0, 1),
                        alpha0 = alpha0.fun(), alpha1 = runif(1, -0.5, -0.45),
                        beta0 = beta0.fun(), omega0 = omega0.fun(), 
                        omega1 = runif(nspec, -5, 5), gamma0 = gamma0.fun(),
                        tau.phi = runif(1, 0, 1)
)

#Parameters to save
params <- c("mu.a0", "tau.a0", "mu.b0", "tau.b0", 
            "mu.o0", "tau.o0", "mu.o1", "tau.o1",
            "mu.g0", "tau.g0", "tau.phi",
            "alpha0", "alpha1", "beta0", "omega0", 
            "omega1", "gamma0")

#MCMC settings
MSdyn.m <- nimbleModel(MSdyn.c, constants = HMSNO.con, data = HMSNO.data, inits = inits())

MCMCconf <- configureMCMC(MSdyn.m, monitors = params)

MCMCconf$addSampler(target = c("alpha0[2,1]", "alpha0[2,2]", "alpha0[6,3]", 
                               "alpha0[2,4]", "alpha0[1,5]", "alpha0[2,5]",
                               "alpha0[2,6]", "alpha0[3,6]", "alpha0[4,6]",
                               "alpha0[5,6]", "alpha0[1,7]", "alpha0[6,8]", 
                               "alpha0[1,9]", "alpha0[2,9]", "alpha0[3,9]",
                               "alpha0[4,9]"),
                    type = "RW_block")

MCMCconf$addSampler(target = c("beta0[2,1]", "beta0[2,2]", "beta0[6,3]", 
                               "beta0[2,4]", "beta0[1,5]", "beta0[2,5]",
                               "beta0[2,6]", "beta0[3,6]", "beta0[4,6]",
                               "beta0[5,6]", "beta0[1,7]", "beta0[6,8]", 
                               "beta0[1,9]", "beta0[2,9]", "beta0[3,9]",
                               "beta0[4,9]"),
                    type = "RW_block")

MCMCconf$addSampler(target = c("omega0[2,1]", "omega0[2,2]", "omega0[6,3]", 
                               "omega0[2,4]", "omega0[1,5]", "omega0[2,5]",
                               "omega0[2,6]", "omega0[3,6]", "omega0[4,6]",
                               "omega0[5,6]", "omega0[1,7]", "omega0[6,8]", 
                               "omega0[1,9]", "omega0[2,9]", "omega0[3,9]",
                               "omega0[4,9]"),
                    type = "RW_block")

MCMCconf$addSampler(target = c("gamma0[2,1]", "gamma0[2,2]", "gamma0[6,3]", 
                               "gamma0[2,4]", "gamma0[1,5]", "gamma0[2,5]",
                               "gamma0[2,6]", "gamma0[3,6]", "gamma0[4,6]",
                               "gamma0[5,6]", "gamma0[1,7]", "gamma0[6,8]", 
                               "gamma0[1,9]", "gamma0[2,9]", "gamma0[3,9]",
                               "gamma0[4,9]"),
                    type = "RW_block")

species <- seq(1, nspec, 1)

nFun <- function(node, i, s){
  paste(node, "[", i, ",", s, "]", sep = "")
}

omit <- NULL

for(s in 1:nspec){
  omit <- c(omit, s)
  include <- parkspec[1:npark[s],s]
  for(i in parkspec[1:npark[s],s]){
    include <- include[-1]
    if(i != parkspec[npark[s],s]){
      for(ii in include){
        # print(c(nFun("alpha0", i, s), nFun("alpha0", ii, s)))
        # print(c(nFun("beta0", i, s), nFun("beta0", ii, s)))
        # print(c(nFun("omega0", i, s), nFun("omega0", ii, s)))
        # print(c(nFun("gamma0", i, s), nFun("gamma0", ii, s)))
        MCMCconf$addSampler(target = c(nFun("alpha0", i, s), nFun("alpha0", ii, s)),
                            type = "RW_block")
        MCMCconf$addSampler(target = c(nFun("beta0", i, s), nFun("beta0", ii, s)),
                            type = "RW_block")
        MCMCconf$addSampler(target = c(nFun("omega0", i, s), nFun("omega0", ii, s)),
                            type = "RW_block")
        MCMCconf$addSampler(target = c(nFun("gamma0", i, s), nFun("gamma0", ii, s)),
                            type = "RW_block")
      }#end ii
    }#end if
    if(s != nspec){
      for(ss in species[-omit]){
        for(ii in parkspec[1:npark[ss],ss]){
          # print(c(nFun("alpha0", i, s), nFun("alpha0", ii, ss)))
          # print(c(nFun("beta0", i, s), nFun("beta0", ii, ss)))
          # print(c(nFun("omega0", i, s), nFun("omega0", ii, ss)))
          # print(c(nFun("gamma0", i, s), nFun("gamma0", ii, ss)))
          MCMCconf$addSampler(target = c(nFun("alpha0", i, s), nFun("alpha0", ii, ss)),
                              type = "RW_block")
          MCMCconf$addSampler(target = c(nFun("beta0", i, s), nFun("beta0", ii, ss)),
                              type = "RW_block")
          MCMCconf$addSampler(target = c(nFun("omega0", i, s), nFun("omega0", ii, ss)),
                              type = "RW_block")
          MCMCconf$addSampler(target = c(nFun("gamma0", i, s), nFun("gamma0", ii, ss)),
                              type = "RW_block")
        }#end ii
      }#end ss
    }#end if
    # print(c("mu.a0", nFun("alpha0", i, s)))
    # print(c("mu.b0", nFun("beta0", i, s)))
    # print(c("mu.o0", nFun("omega0", i, s)))
    # print(c("mu.g0", nFun("gamma0", i, s)))
    MCMCconf$addSampler(target = c("mu.a0", nFun("alpha0", i, s)),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c("mu.b0", nFun("beta0", i, s)),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c("mu.o0", nFun("omega0", i, s)),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c("mu.g0", nFun("gamma0", i, s)),
                        type = "AF_slice")
    MCMCconf$addSampler(target = c(nFun("omega0", i, s), paste("omega1", "[", s, "]", sep = "")),
                        type = "AF_slice")
  }#end i
}#end s

MCMC <- buildMCMC(MCMCconf)

MSdyn.comp <- compileNimble(MSdyn.m, MCMC)

ni <- 20000
nb <- 10000
nc <- 3
nt <- 5

#Run NIMBLE model
MSdyn.o <- runMCMC(MSdyn.comp$MCMC, niter = ni, nburnin = nb, nchains = nc, thin = nt, 
                   setSeed = c(1,2,3), samplesAsCodaMCMC = TRUE)

#-Save output-#
ID <- gsub(" ","_",Sys.time())
ID <- gsub(":", "-", ID)
save(MSdyn.o, file = paste("output", ID, ".Rdata", sep=""))

traceplot(MSdyn.o[c(1:3)][,attr(MSdyn.o$chain1, "dimnames")[[2]][543]])

param.list <- c("mu.a0", "tau.a0", "mu.b0", "tau.b0", 
                "mu.o0", "tau.o0", "mu.o1", "tau.o1",
                "mu.g0", "tau.g0", "tau.phi", "alpha1",
  "alpha0[2, 1]", "alpha0[2, 2]", "alpha0[6, 3]", 
  "alpha0[2, 4]", "alpha0[1, 5]", "alpha0[2, 5]",
  "alpha0[2, 6]", "alpha0[3, 6]", "alpha0[4, 6]",
  "alpha0[5, 6]", "alpha0[1, 7]", "alpha0[6, 8]", 
  "alpha0[1, 9]", "alpha0[2, 9]", "alpha0[3, 9]",
  "alpha0[4, 9]", "beta0[2, 1]", "beta0[2, 2]", "beta0[6, 3]", 
  "beta0[2, 4]", "beta0[1, 5]", "beta0[2, 5]",
  "beta0[2, 6]", "beta0[3, 6]", "beta0[4, 6]",
  "beta0[5, 6]", "beta0[1, 7]", "beta0[6, 8]", 
  "beta0[1, 9]", "beta0[2, 9]", "beta0[3, 9]",
  "beta0[4, 9]", "omega0[2, 1]", "omega0[2, 2]", "omega0[6, 3]", 
  "omega0[2, 4]", "omega0[1, 5]", "omega0[2, 5]",
  "omega0[2, 6]", "omega0[3, 6]", "omega0[4, 6]",
  "omega0[5, 6]", "omega0[1, 7]", "omega0[6, 8]", 
  "omega0[1, 9]", "omega0[2, 9]", "omega0[3, 9]",
  "omega0[4, 9]", "omega1[1]", "omega1[2]", "omega1[3]",
  "omega1[4]", "omega1[5]", "omega1[6]", "omega1[7]",
  "omega1[8]", "omega1[9]",
  "gamma0[2, 1]", "gamma0[2, 2]", "gamma0[6, 3]", 
  "gamma0[2, 4]", "gamma0[1, 5]", "gamma0[2, 5]",
  "gamma0[2, 6]", "gamma0[3, 6]", "gamma0[4, 6]",
  "gamma0[5, 6]", "gamma0[1, 7]", "gamma0[6, 8]", 
  "gamma0[1, 9]", "gamma0[2, 9]", "gamma0[3, 9]",
  "gamma0[4, 9]"
)

traceplot(MSdyn.o[c(1:3)][,param.list])
gelman.diag(MSdyn.o[c(1:3)][,param.list])
effectiveSize(MSdyn.o[c(1:3)][,param.list])
summary(MSdyn.o[c(1:3)][,param.list])
#traceplot(out[c(1:3)][,"Npark[1, 1, 1]"])
#plot(unlist(output[c(1:5)][,"beta0[1]"]), unlist(output[c(1:5)][,"beta0[4]"])

for(i in 1:npark){
  plot(nyrstr[i], mean(unlist(output[c(1:5)][,paste("Npark[", nyrstr[i], ", ", i, ", 8]", sep = "")])),
       xlim = c(1,9), ylim = c(0,500), xlab = "Abundance", ylab = "Year")
  for(t in (nyrstr[i] + 1):nyrend[i]){
    points(t, mean(unlist(output[c(1:5)][,paste("Npark[", t, ", ", i, ", 8]", sep = "")])))
  }
}
