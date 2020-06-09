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
library(jagsUI)

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
nspec <- as.numeric(Data %>% summarize (n_distinct(species)))
nobs <- length(y)

#Nested indices
yr <- as.numeric(as.factor(Data2$Year))
site <- as.numeric(as.factor(Data2$`deployment ID`))
spec <- as.numeric(as.factor(Data2$species))

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

elev[is.na(elev)] <- 0
edge[is.na(edge)] <- 0
days[is.na(days)] <- 0

#--------------#
#-Compile Data-#
#--------------#

#Data
jags.data <- list(y = y, nobs = nobs, nstart = nstart, nend = nend, nsites = nsites, nspec = nspec,
                  yr = yr, site = site, spec = spec, elev = elev, edge = edge, density = density, days = days)

#Initial values
Nst <- array(NA, dim = c(nyrs, nsites, nspec))
Sst <- array(NA, dim = c(nyrs-1, nsites, nspec))
Gst <- array(NA, dim = c(nyrs-1, nsites, nspec))
for(j in 1:nsites){
  Nst[nstart[j],j,] <- 10
  Sst[nstart[j]:(nend[j]-1),j,] <- 5
  Gst[nstart[j]:(nend[j]-1),j,] <- 2
}

alpha0.fun <- function(){
  alpha0 <- NULL
  alpha0[1] <- runif(1, -0.2, 0.2)
  alpha0[2] <- runif(1, -4, -2.75)
  alpha0[3] <- runif(1, -2, -1.5)
  alpha0[4] <- runif(1, -5, -3)
  alpha0[5] <- runif(1, -1, -0.6)
  alpha0[6] <- runif(1, -1.5, -1)
  alpha0[7] <- runif(1, -2, -1.5)
  alpha0[8] <- runif(1, -4, -2)
  alpha0[9] <- runif(1, -2, -1.5)
  return(alpha0)
}

beta0.fun <- function(){
  beta0 <- NULL
  beta0[1] <- runif(1, 0.4, 0.8)
  beta0[2] <- runif(1, 1.5, 3.5)
  beta0[3] <- runif(1, 0, 1)
  beta0[4] <- runif(1, -7, -4)
  beta0[5] <- runif(1, 0.5, 1)
  beta0[6] <- runif(1, -0.5, 0)
  beta0[7] <- runif(1, 1, 3)
  beta0[8] <- runif(1, 0, 3)
  beta0[9] <- runif(1, 0.5, 1)
  return(beta0)
}

beta4.fun <- function(){
  beta4 <- NULL
  beta4[1] <- runif(1, -0.5, 0.5)
  beta4[2] <- runif(1, 0, 0.5)
  beta4[3] <- runif(1, -0.5, 0)
  beta4[4] <- runif(1, -2, 2)
  beta4[5] <- runif(1, 0, 0.6)
  beta4[6] <- runif(1, 1, 1.5)
  beta4[7] <- runif(1, 0, 2)
  beta4[8] <- runif(1, 1, 2)
  beta4[9] <- runif(1, 0.4, 0.6)
  return(beta4)
}

omega0.fun <- function(){
  omega0 <- NULL
  omega0[1] <- runif(1, 1.4, 2.2)
  omega0[2] <- runif(1, -5, -3)
  omega0[3] <- runif(1, 0.5, 1)
  omega0[4] <- runif(1, -2, -0.4)
  omega0[5] <- runif(1, 0.6, 1)
  omega0[6] <- runif(1, 2, 4)
  omega0[7] <- runif(1, 0, 1)
  omega0[8] <- runif(1, 0, 3)
  omega0[9] <- runif(1, 2, 3)
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

inits <- function()list(N=Nst, S=Sst, G=Gst,
                        mu.a0 = runif(1, 0.1, 0.25), tau.a0 = runif(1, 0, 5), alpha0 = alpha0.fun(), alpha1 = runif(1, -0.25, 0.25),
                        mu.b0 = runif(1, 0, 3), tau.b0 = runif(1, 0, 1), beta0 = beta0.fun(),
                        mu.b1 = runif(1, -1, 1), tau.b1 = runif(1, 0, 10), beta1 = runif(nspec, -1, 1),
                        mu.b2 = runif(1, -1, 1), tau.b2 = runif(1, 0, 10), beta2 = runif(nspec, -0.5, 0.5),
                        mu.b3 = runif(1, -1, 1), tau.b3 = runif(1, 0, 10), beta3 = runif(nspec, -1, 1),
                        mu.b4 = runif(1, 0, 2), tau.b4 = runif(1, 0, 2), beta4 = beta4.fun(),
                        mu.o0 = runif(1, 0.5, 0.9), tau.o0 = runif(1, 0, 1), omega0 = omega0.fun(),
                        mu.o1 = runif(1, 0, 4), tau.o1 = runif(1, 0, 1), omega1 = runif(nspec, -5, 5),
                        mu.o2 = runif(1, -1, 1), tau.o2 = runif(1, 0, 10), omega2 = runif(nspec, -1, 1),
                        mu.o3 = runif(1, -1, 1), tau.o3 = runif(1, 0, 10), omega3 = runif(nspec, -1, 1),
                        gamma0 = runif(1, -4, 0)
)

#Parameters to save
params <- c("mu.a0", "tau.a0", 
            "mu.b0", "tau.b0", "mu.b1", "tau.b1", "mu.b2", "tau.b2", "mu.b3", "tau.b3", "mu.b4", "tau.b4",
            "mu.o0", "tau.o0", "mu.o1", "tau.o1", "mu.o2", "tau.o2", "mu.o3", "tau.o3",
            "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "beta4",
            "omega0", "omega1", "omega2", "omega3", "gamma0")


#MCMC settings
# ni <- 115000
# na <- 1000
# nb <- 100000
# nc <- 3
# nt <- 5

na <- 1000
ni <- 100000
nb <- 20000
nc <- 1
nt <- 40

#Run JAGS model
out <- jagsUI(jags.data, inits, params, model.file = "../MSdyn_cov.txt", 
              n.adapt = na, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, 
              parallel = FALSE)
