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
nstart <- as.numeric(nyrs$`min(Year)`) - 2009 + 1
nend <- 9 - (2017 - as.numeric(nyrs$`max(Year)`))

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

#--------------#
#-Compile Data-#
#--------------#

#Data
jags.data <- list(y = y, nobs = nobs, nstart = nstart, nend = nend, nsites = nsites, nspec = nspec, 
                  yr = yr, site = site, spec = spec, elev = elev, edge = edge, density = density, days = days)

#Initial values
Nst <- array(NA, dim = c(9, nsites, nspec))
Sst <- array(NA, dim = c(9-1, nsites, nspec))
Gst <- array(NA, dim = c(9-1, nsites, nspec))
for(j in 1:nsites){
Nst[nstart[j],j,] <- 10
Sst[nstart[j]:(nend[j]-1),j,] <- 5
Gst[nstart[j]:(nend[j]-1),j,] <- 2
}

inits <- function()list(N=Nst, S=Sst, G=Gst)

#Parameters to save
params <- c("mu.a0", "tau.a0", 
            "mu.b0", "sig.b0", "mu.b1", "tau.b1", "mu.b2", "tau.b2", "mu.b3", "tau.b3", "mu.b4", "tau.b4",
            "mu.o0", "sig.o0", "mu.o1", "sig.o1", "mu.o2", "sig.o2", "mu.o3", "sig.o3",
            "mu.g0", "sig.g0",
            "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "beta4",
            "omega0", "omega1", "omega2", "omega3", "gamma0")

#MCMC settings
ni <- 115000
na <- 1000
nb <- 100000
nc <- 3
nt <- 5

#Run JAGS model
out <- jagsUI(jags.data, inits, params, model.file = "../MSdyn_cov.txt", 
              n.adapt = na, n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, 
              parallel = TRUE)
