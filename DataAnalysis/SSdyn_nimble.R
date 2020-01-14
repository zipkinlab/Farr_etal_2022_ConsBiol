#-------------------#
#-Working Directory-#
#-------------------#

setwd("./DataFormat")

#----------------#
#-Load Libraries-#
#----------------#

library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(nimble)
library(compareMCMCs)
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

#-------------#
#-Format Data-#
#-------------#

#Convert to NAs
Data[,5:10] <- apply(Data[,5:10], 2, function(y) as.numeric(gsub("-", NA, y)))

#Remove cameras that were sampled in only 1 year
Data <- Data %>% group_by(`deployment ID`) %>% filter(n() > 2) %>% ungroup()

#Filter for single-species
Data <- Data %>% filter(species == "nigrifrons")

#Convert to single vector format
Data2 <- melt(Data, id = c("park", "deployment ID", "Year", "species", "days", "elevation", "edge"))

#Occupancy data
y <- Data2$value

#Indices
nyrs <- Data %>% group_by(`deployment ID`) %>% summarize(min(Year), max(Year))
nstart <- as.numeric(nyrs$`min(Year)`) - as.numeric(min(nyrs$`min(Year)`)) + 1
nend <- 1 + as.numeric(max(nyrs$`max(Year)`)) - as.numeric(min(nyrs$`min(Year)`)) - 
  (as.numeric(max(nyrs$`max(Year)`)) - as.numeric(nyrs$`max(Year)`))

nyrs <- as.numeric(max(nyrs$`max(Year)`)) - as.numeric(min(nyrs$`min(Year)`)) + 1

nsites <- as.numeric(Data %>% summarize(n_distinct(`deployment ID`)))
nobs <- length(y)

#Nested indices
yr <- as.numeric(as.factor(Data2$Year))
site <- as.numeric(as.factor(Data2$`deployment ID`))

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

SSdyn.c <- nimbleCode({
  
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  beta2 ~ dnorm(0, 0.1)
  beta3 ~ dnorm(0, 0.1)
  beta4 ~ dnorm(0, 0.1)
  
  alpha0 ~ dunif(0, 1)
  alpha1 ~ dnorm(0, 0.1)
  
  omega0 ~ dunif(0, 1)
  omega1 ~ dnorm(0, 0.1)
  omega2 ~ dnorm(0, 0.1)
  omega3 ~ dnorm(0, 0.1)
  
  gamma0 ~ dnorm(0, 0.1)
  
  for(j in 1:nsites){
    
    log(lambda[j]) <- beta0 + beta1 * density[j] + beta2 * edge[nstart[j],j] +
      beta3 * density[j] * edge[nstart[j],j] + beta4 * elev[nstart[j],j]
    
    logit(r[nstart[j],j]) <- logit(alpha0) + alpha1 * days[nstart[j],j]
    N[nstart[j],j] ~ dpois(lambda[j])
    
    for(t in (nstart[j]+1):nend[j]){
      
      logit(r[t,j]) <- alpha0 + alpha1 * days[t,j]
      
      logit(omega[t-1,j]) <- logit(omega0) + omega1 * density[j] + omega2 * edge[t,j] + 
        omega3 * density[j] * edge[t,j]
      
      S[t-1,j] ~ dbin(omega[t-1,j], N[t-1,j])
      G[t-1,j] ~ dpois(gamma[t-1])
      N[t,j] <- S[t-1,j] + G[t-1,j]
      
    }#end t
  }#end j
  
  for(t in 2:nyrs){ #check 9
    log(gamma[t-1]) <- gamma0
  }#end t
  
  for(k in 1:nobs){
    y[k] ~ dbern(p[k])
    p[k] <- 1 - pow((1 - r[yr[k], site[k]]), N[yr[k], site[k]])
  }#end k
  
})

#--------------#
#-Compile Data-#
#--------------#

#Data
SSdyn.d <- list(y = y, nobs = nobs, nstart = nstart, nend = nend, nsites = nsites, nyrs = nyrs,
                  yr = yr, site = site, elev = elev, edge = edge, density = density, days = days)

#Initial values
Nst <- array(NA, dim = c(nyrs, nsites))
Sst <- array(NA, dim = c(nyrs-1, nsites))
Gst <- array(NA, dim = c(nyrs-1, nsites))
for(j in 1:nsites){
  Nst[nstart[j],j] <- 10
  Sst[nstart[j]:(nend[j]-1),j] <- 5
  Gst[nstart[j]:(nend[j]-1),j] <- 2
}

inits <- function()list(N=Nst, S=Sst, G=Gst)

#Parameters to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "beta4",
            "omega0", "omega1", "omega2", "omega3", "gamma0")

#MCMC settings
SSdyn.m <- nimbleModel(SSdyn.c, constants = SSdyn.d, inits = inits())

MCMCconf <- configureMCMC(SSdyn.m)
MCMCconf$printSamplers(1:12)

MCMCconf$removeSampler(c('beta0', 'beta1', 'beta2', 'beta3', 'beta4',
                         'omega0', 'omega1', 'omega2', 'omega3'))

MCMCconf$addSampler(target = c('beta0', 'beta1', 'beta2', 'beta3', 'beta4'), 
                    type = "AF_slice")

MCMCconf$addSampler(target = c('omega0', 'omega1', 'omega2', 'omega3'), 
                    type = "AF_slice")

MCMCconf$samplerExecutionOrder <- c(seq(3,4166,1),1,2)

MCMC <- buildMCMC(MCMCconf)

SSdyn.comp <- compileNimble(SSdyn.m, MCMC)

ni <- 25000
nb <- 20000
nc <- 3

#Run NIMBLE model
SSdyn.o <- runMCMC(SSdyn.comp$MCMC, niter = ni, nburnin = nb, nchains = nc)
SSdyn.o2 <- nimbleMCMC(SSdyn.c, inits = inits, constants = SSdyn.d,
                       niter = ni, nburnin = nb)

#------------#
#-Post check-#
#------------#

samples <- as.mcmc.list(lapply(SSdyn.o, mcmc))

gelman.diag(samples)
effectiveSize(samples)
plot(samples)
