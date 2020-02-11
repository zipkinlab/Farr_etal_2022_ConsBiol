#-------------------#
#-Working Directory-#
#-------------------#

#setwd("./DataFormat")

#----------------#
#-Load Libraries-#
#----------------#

library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(nimble)
#library(compareMCMCs)
library(coda)

#-----------#
#-Load Data-#
#-----------#

Data <- NULL
files <- list.files(path = "./EditedData", pattern = "Data", full.names = TRUE)
for(i in 1:length(files)){
  sheets <- excel_sheets(files[i])
  for(j in 1:length(sheets)){
    Temp <- read_excel(files[i], sheet = sheets[j])
    Temp <- Temp %>% mutate(species = as.character(sheets[j]), park = as.factor(i)) %>% 
      select(park, `deployment ID`, Year, species, ls(Temp, pattern = "R"), density)
    Data <- bind_rows(Data, Temp)
  }
}

covfile <- "./RawData/camera trap metadata.xlsx"
covsheets <- excel_sheets(covfile)
covsheets <- covsheets[-3] #Remove Nyungwe

CovData <- NULL
for(i in 1:5){
  Temp <- read_excel(covfile, sheet = covsheets[i])
  Temp <- Temp %>% select(`deployment ID`, days, elevation, edge, Year)
  CovData <- bind_rows(CovData, Temp)
}

CovData <- bind_rows(CovData, read_excel("./EditedData/Data_Nyungwe_nig_sylv_April_1.xlsx", sheet = 1) %>%
                       select(`deployment ID`, days, elevation, edge, Year)) %>%
  bind_rows(., read_excel("./EditedData/Data_Nyungwe_nig_sylv_April_1.xlsx", sheet = 2) %>%
              select(`deployment ID`, days, elevation, edge, Year)) %>%
  distinct()

Data <- full_join(Data, CovData, by = c("deployment ID", "Year")) %>% drop_na()

rm(Temp, CovData)

#-------------#
#-Format Data-#
#-------------#

#Convert to NAs
Data[,5:10] <- apply(Data[,5:10], 2, function(y) as.numeric(gsub("-", NA, y)))

#Remove cameras that were sampled in only 1 year
Data <- Data %>% group_by(`deployment ID`) %>% filter(n() > 2) %>% ungroup()

#Arrange dataframe
Data <- Data %>% arrange(park, Year, `deployment ID`)

#Convert to single vector format
Data2 <- melt(Data, id = c("park", "deployment ID", "Year", "species", "days", "elevation", "edge"))

#Occupancy data
y <- Data2$value

#Indices
df <- Data %>% group_by(`deployment ID`, park) %>% summarize(minYear = min(Year), maxYear = max(Year))
df$parkID <- df %>% select(park) %>% .$park %>% as.numeric(.)
df$siteID <- df %>% select(`deployment ID`) %>% .$`deployment ID` %>% as.factor(.) %>% as.numeric(.)
df$nstart <- as.numeric(df$minYear) - as.numeric(min(df$minYear)) + 1
df$nend <- 1 + as.numeric(max(df$maxYear)) - as.numeric(min(df$minYear)) - 
  (as.numeric(max(df$maxYear)) - as.numeric(df$maxYear))

nyrs <- as.numeric(max(df$maxYear)) - as.numeric(min(df$minYear)) + 1
nsites <- as.numeric(Data %>% summarize(n_distinct(`deployment ID`)))
nspec <- as.numeric(Data %>% summarize (n_distinct(species)))
nobs <- length(y)
npark <- max(as.numeric(df$park))
nstart <- df$nstart
nend <- df$nend

parkdf <- Data %>% mutate(siteID = as.numeric(as.factor(`deployment ID`))) %>% group_by(park) %>% 
  summarize(minYear = min(Year), maxYear = max(Year), 
            minSite = min(siteID),
            maxSite = max(siteID))
parkdf$nstart <- as.numeric(parkdf$minYear) - as.numeric(min(parkdf$minYear)) + 1
parkdf$nend <- 1 + as.numeric(max(parkdf$maxYear)) - as.numeric(min(parkdf$minYear)) - 
  (as.numeric(max(parkdf$maxYear)) - as.numeric(parkdf$maxYear))
nyrstr <- parkdf$nstart
nyrend <- parkdf$nend
nsitestr <- parkdf$minSite
nsiteend <- parkdf$maxSite

#Nested indices
Data2 <- Data2 %>% left_join(., df %>% select(`deployment ID`, park, parkID, siteID), by = c("deployment ID", "park"))

yr <- as.numeric(as.factor(Data2$Year))
site <- Data2$siteID
spec <- as.numeric(as.factor(Data2$species))
park <- as.numeric(df$park)

# 
# #Site indices
# site1 <- cbind(as.factor(unique(Data[,c("deployment ID", "park", "Year")])$Year), 
#               as.factor(unique(Data[,c("deployment ID", "park", "Year")])$park), 
#               as.factor(unique(Data[,c("deployment ID", "park", "Year")])$`deployment ID`))
# site2 <- list()
# for(i in 1:npark){
#   for(t in 1:nyrs){
#     print(c(i,t))
#     print(site1[which(site1[,1]==t&site1[,2]==i),3])
#   }
# }

#Covariates
Cov <- Data %>% group_by(`deployment ID`, Year) %>%
  distinct(`deployment ID`, Year, .keep_all = TRUE)
Cov$elevation <- scale(as.numeric(Cov$elevation))
Cov$edge <- scale(as.numeric(Cov$edge))
Cov$days <- scale(as.numeric(Cov$days))

elev <- dcast(Cov %>% select(`deployment ID`, Year, elevation), `deployment ID` ~ Year)
elev <- t(as.matrix(elev[,-1]))
edge <- dcast(Cov %>% select(`deployment ID`, Year, edge), `deployment ID` ~ Year)
edge <- t(as.matrix(edge[,-1]))
days <- dcast(Cov %>% select(`deployment ID`, Year, days), `deployment ID` ~ Year)
days <- t(as.matrix(days[,-1]))
density <- Data %>% distinct(`deployment ID`, density)
density <- as.numeric(density$density)

#--------------#
#-NIMBLE model-#
#--------------#

MSdyn.c <- nimbleCode({
  
  mu.b0 ~ dnorm(0, 0.1)
  # #tau.b0 <- 1/(sig.b0 * sig.b0)
  # #sig.b0 ~ dt(0, pow(2.5,-2), 1) T(0,)
  sig.b0 ~ dunif(0, 3)
  # tau.b0 ~ dgamma(0.1, 0.1)
  # mu.b1 ~ dnorm(0, 0.1)
  # tau.b1 ~ dgamma(0.1, 0.1)
  # mu.b2 ~ dnorm(0, 0.1)
  # tau.b2 ~ dgamma(0.1, 0.1)
  # mu.b3 ~ dnorm(0, 0.1)
  # tau.b3 ~ dgamma(0.1, 0.1)
  # mu.b4 ~ dnorm(0, 0.1)
  # tau.b4 ~ dgamma(0.1, 0.1)
  # mu.a0 ~ dnorm(0, 0.368) #dunif(0, 1)
  # mu.a0L <- logit(mu.a0)
  # mu.a0 ~ dunif(0, 1)
  # tau.a0 ~ dgamma(0.1, 0.1)
  #mu.o0 ~ dnorm(0, 0.368) #dunif(0, 1)
  mu.o0L <- logit(mu.o0)
  mu.o0 ~ dunif(0, 1)
  #tau.o0 <- 1/(sig.o0 * sig.o0)
  #sig.o0 ~ dt(0, pow(2.5,-2), 1) T(0,)
  tau.o0 ~ dgamma(0.1, 0.1)
  # #mu.o1 ~ dnorm(0, 0.368)
  # mu.o1 ~ dnorm(0, 0.1)
  # #tau.o1 <- 1/(sig.o1 * sig.o1)
  # #sig.o1 ~ dt(0, pow(2.5,-2), 1) T(0,)
  # tau.o1 ~ dgamma(0.1, 0.1)
  # #mu.o2 ~ dnorm(0, 0.368)
  # mu.o2 ~ dnorm(0, 0.1)
  # #tau.o2 <- 1/(sig.o2 * sig.o2)
  # #sig.o2 ~ dt(0, pow(2.5,-2), 1) T(0,)
  # tau.o2 ~ dgamma(0.1, 0.1)
  # #mu.o3 ~ dnorm(0, 0.368)
  # mu.o3 ~ dnorm(0, 0.1)
  # #tau.o3 <- 1/(sig.o3 * sig.o3)
  # #sig.o3 ~ dt(0, pow(2.5,-2), 1) T(0,)
  # tau.o3 ~ dgamma(0.1, 0.1)
  # #mu.g0 ~ dnorm(0, 0.1)
  # #tau.g0 <- 1/(sig.g0 * sig.g0)
  # #sig.g0 ~ dt(0, pow(2.5,-2), 1) T(0,)
  # #tau.g0 ~ dgamma(0.1, 0.1)
  
  log(gamma) <- gamma0
  gamma0 ~ dnorm(0, 0.1)
  
  #alpha1 ~ dnorm(0, 0.368)
  alpha0 ~ dunif(0, 1)
  alpha1 ~ dnorm(0, 0.1)
  
  
  for(i in 1:npark){
    eps[i] ~ dnorm(0, tau.p)
  }
  
  tau.p ~ dgamma(0.1, 0.1)
  
  for(s in 1:nspec){
    
    beta0[s] ~ dnorm(mu.b0, sd = sig.b0) 
    # beta1[s] ~ dnorm(mu.b1, tau.b1)
    # beta2[s] ~ dnorm(mu.b2, tau.b2)
    # beta3[s] ~ dnorm(mu.b3, tau.b3)
    # beta4[s] ~ dnorm(mu.b4, tau.b4)
    # alpha0[s] ~ dnorm(mu.a0L, tau.a0)
    omega0[s] ~ dnorm(mu.o0L, tau.o0)
    # omega1[s] ~ dnorm(mu.o1, tau.o1)
    # omega2[s] ~ dnorm(mu.o2, tau.o2)
    # omega3[s] ~ dnorm(mu.o3, tau.o3)
    #gamma0[s] ~ dnorm(mu.g0, tau.g0)
    
    for(j in 1:nsites){
      
      logit(r[nstart[j],j,s]) <- logit(alpha0) + alpha1 * days[nstart[j],j]
      
      N[nstart[j],j,s] ~ dpois(lambda[j,s])
      log(lambda[j,s]) <- beta0[s] + eps[park[j]]
      #+ beta1[s] * density[j] + beta2[s] * edge[nstart[j],j] +
      #   beta3[s] * density[j] * edge[nstart[j],j] + beta4[s] * elev[nstart[j],j]
      
      for(t in (nstart[j]+1):nend[j]){
        
        logit(r[t,j,s]) <- alpha0 + alpha1 * days[t,j]
        
        # logit(omega[t-1,j,s]) <- omega0[s] + omega1[s] * density[j] + omega2[s] * edge[t,j] +
        #   omega3[s] * density[j] * edge[t,j]
        
        logit(omega[t-1,j,s]) <- omega0[s] #+ omega1[s] * density[j] + omega2[s] * edge[t,j] +
        #omega3[s] * density[j] * edge[t,j]
        
        S[t-1,j,s] ~ dbin(omega[t-1,j,s], N[t-1,j,s])
        #G[t-1,j,s] ~ dpois(gamma[t-1,s])
        G[t-1,j,s] ~ dpois(gamma)
        N[t,j,s] <- S[t-1,j,s] + G[t-1,j,s]
      }#end t
    }#end j
    
    #for(t in 2:nyrs){ #check 9
    #  log(gamma[t-1,s]) <- gamma0[s]
    #}#end t
    for(i in 1:npark){
      for(t in nyrstr[i]:nyrend[i]){
        Npark[t,i,s] <- sum(N[t,nsitestr[i]:nsiteend[i],s])
      }#end t
    }#end i
  }#end s
  
  for(k in 1:nobs){
    y[k] ~ dbern(p[k])
    p[k] <- 1 - pow((1 - r[yr[k], site[k], spec[k]]), N[yr[k], site[k], spec[k]])
  }#end k
  
})

#--------------#
#-Compile Data-#
#--------------#

#Data
MSdyn.data <- list(y = y)
MSdyn.con <- list(nobs = nobs, nstart = nstart, nend = nend, nsites = nsites, nspec = nspec, nyrs = nyrs, npark = npark,
                  yr = yr, site = site, spec = spec, park = park, elev = elev, edge = edge, density = density, days = days,
                  nyrstr = nyrstr, nyrend = nyrend, nsitestr = nsitestr, nsiteend = nsiteend)

#Initial values
Nst <- array(NA, dim = c(nyrs, nsites, nspec))
Sst <- array(NA, dim = c(nyrs-1, nsites, nspec))
Gst <- array(NA, dim = c(nyrs-1, nsites, nspec))
for(j in 1:nsites){
  Nst[nstart[j],j,] <- 10
  Sst[nstart[j]:(nend[j]-1),j,] <- 5
  Gst[nstart[j]:(nend[j]-1),j,] <- 2
}

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

beta0.fun <- function(){
  beta0 <- NULL
  beta0[1] <- runif(1, 1.5, 2.5)
  beta0[2] <- runif(1, 0.5, 1.25)
  beta0[3] <- runif(1, 0, 0.75)
  beta0[4] <- runif(1, -5, -3)
  beta0[5] <- runif(1, 1.5, 2.25)
  beta0[6] <- runif(1, -0.25, 0.5)
  beta0[7] <- runif(1, 0.75, 1.5)
  beta0[8] <- runif(1, -0.75, 0)
  beta0[9] <- runif(1, -0.25, 0.5)
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
                        mu.b0 = runif(1, -0.3, 0.75), sig.b0 = runif(1, 1, 2), beta0 = beta0.fun(),
                        mu.o0 = runif(1, 0.5, 0.75), tau.o0 = runif(1, 0.5, 1.25), omega0 = omega0.fun(),
                        gamma0 = runif(1, -2, -1.75), tau.p = runif(1, 0.5, 2)
)

#Parameters to save
# params <- c("mu.a0", "tau.a0", 
#             "mu.b0", "tau.b0", "mu.b1", "tau.b1", "mu.b2", "tau.b2", "mu.b3", "tau.b3", "mu.b4", "tau.b4",
#             "mu.o0", "tau.o0", "mu.o1", "tau.o1", "mu.o2", "tau.o2", "mu.o3", "tau.o3",
#             "mu.g0", "tau.g0",
#             "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "beta4",
#             "omega0", "omega1", "omega2", "omega3", "gamma0")

params <- c("mu.b0", "sig.b0", "mu.o0", "tau.o0",
            "alpha0", "alpha1", "beta0",
            "omega0", "gamma0", "tau.p", "Npark")

#MCMC settings
MSdyn.m <- nimbleModel(MSdyn.c, constants = MSdyn.con, data = MSdyn.data, inits = inits())

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

MCMCconf$addSampler(target = c("beta0[1]", "beta0[2]", "beta0[3]", "beta0[4]", "beta0[5]",
                               "beta0[6]", "beta0[7]", "beta0[8]", "beta0[9]"),
                    type = "RW_block")

MCMCconf$addSampler(target = c("omega0[1]", "omega0[2]", "omega0[3]", "omega0[4]", "omega0[5]",
                               "omega0[6]", "omega0[7]", "omega0[8]", "omega0[9]"),
                    type = "RW_block")

for(i in 1:nspec){
  omit <- c(omit, i)
  if(i != nspec){
    for(j in species[-omit]){
      MCMCconf$addSampler(target = c(nFun("beta0", i), nFun("beta0", j)),
                          type = "RW_block") 
      MCMCconf$addSampler(target = c(nFun("omega0", i), nFun("omega0", j)),
                          type = "RW_block")
    } 
  }
  MCMCconf$addSampler(target = c("mu.b0", nFun("beta0", i)),
                      type = "RW_block") 
  MCMCconf$addSampler(target = c("mu.o0", nFun("omega0", i)),
                      type = "AF_slice")  
}


#MCMCconf$addSampler(target = c('mu.b0', 'beta0'), type = "RW_block")

#MCMCconf$addSampler(target = c('mu.o0', 'omega0'), type = "AF_slice")

MCMC <- buildMCMC(MCMCconf)

MSdyn.comp <- compileNimble(MSdyn.m, MCMC)

ni <- 200000
nb <- 100000
nc <- 1
nt <- 50

#Run NIMBLE model
MSdyn.o <- runMCMC(MSdyn.comp$MCMC, niter = ni, nburnin = nb, nchains = nc, thin = nt, samplesAsCodaMCMC = TRUE)

#-Save output-#
ID <- gsub(" ","_",Sys.time())
ID <- gsub(":", "-", ID)
save(MSdyn.o, file = paste("output", ID, ".Rdata", sep=""))

#traceplot(out[c(1:3)][,"Npark[1, 1, 1]"])
#plot(unlist(output[c(1:5)][,"beta0[1]"]), unlist(output[c(1:5)][,"beta0[4]"])

# for(i in 1:npark){
#   plot(nyrstr[i], mean(unlist(output[c(1:5)][,paste("Npark[", nyrstr[i], ", ", i, ", 8]", sep = "")])), 
#        xlim = c(1,9), ylim = c(0,500), xlab = "Abundance", ylab = "Year")
#   for(t in (nyrstr[i] + 1):nyrend[i]){
#     points(t, mean(unlist(output[c(1:5)][,paste("Npark[", t, ", ", i, ", 8]", sep = "")])))
#   }
# }
