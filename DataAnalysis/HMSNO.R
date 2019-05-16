#-------------------#
#-Working Directory-#
#-------------------#

setwd("C:/Users/farrm/Documents/GitHub/HMSNO/DataFormat")

#----------------#
#-Load Libraries-#
#----------------#

library(readxl)
library(dplyr)
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
      select(park, `deployment ID`, Year, species, ls(Temp, pattern = "R"), elevation, edge, density)
    Data <- bind_rows(Data, Temp)
  }
}

#-------------#
#-Format Data-#
#-------------#

#Convert to NAs
Data[,5:10] <- apply(Data[,5:10], 2, function(y) as.numeric(gsub("-", NA, y)))

#Convert to single vector format
Data2 <- melt(Data, id = c("park", "deployment ID", "Year", "species", "elevation", "edge"))

#Occupancy data
y <- Data2$value

#Indices
nyrs <- as.numeric(Data %>% summarize(n_distinct(Year)))
nsites <- as.numeric(Data %>% summarize(n_distinct(`deployment ID`)))
nspec <- as.numeric(Data %>% summarize (n_distinct(species)))
nobs <- length(y)

#Nested indices
yr <- as.numeric(as.factor(Data2$Year))
site <- as.numeric(as.factor(Data2$`deployment ID`))
spec <- as.numeric(as.factor(Data2$species))

#Covariates
elev <- Data %>% group_by(`deployment ID`) %>% distinct(elevation) %>% 
  distinct(`deployment ID`, .keep_all = TRUE)
elev <- as.numeric(elev$elevation)

edge <- Data %>% group_by(`deployment ID`) %>% distinct(edge) %>% 
  distinct(`deployment ID`, .keep_all = TRUE)
edge <- as.numeric(edge$edge)

density <- Data %>% distinct(`deployment ID`, density)
density <- as.numeric(density$density)

#--------------#
#-Compile Data-#
#--------------#

#Data
jags.data <- list(y = y, nobs = nobs, nyrs = nyrs, nsites = nsites, nspec = nspec, 
                  yr = yr, site = site, spec = spec, elev = elev, edge = edge, density = density)

#Initial values
Nst <- array(NA, dim = c(nyrs, nsites, nspec))
Nst[1,,] <- 10
Sst <- array(5, dim = c(nyrs-1, nsites, nspec))
Gst <- array(2, dim = c(nyrs-1, nsites, nspec))

inits <- function()list(N=Nst, S=Sst, G=Gst)

#Parameters to save
params <- c("mu.a0", "tau.a0", 
            "mu.b0", "sig.b0", "mu.b1", "tau.b1", "mu.b2", "tau.b2", "mu.b3", "tau.b3", "mu.b4", "tau.b4",
            "mu.o0", "sig.o0", "mu.o1", "sig.o1", "mu.o2", "sig.o2", "mu.o3", "sig.o3",
            "mu.g0", "sig.g0",
            "alpha0", "beta0", "beta1", "beta2", "beta3", "beta4",
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
