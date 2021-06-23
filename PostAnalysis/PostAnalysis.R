#TO DO: effort by year

#-----------#
#-Libraries-#
#-----------#

library(coda)
library(tidyverse)
library(sf)
library(ggthemes)
library(reshape2)
library(grid)
library(gridExtra)
library(readxl)
library(nimble)

#-----------#
#-Load data-#
#-----------#

load(file = "L:/Farr/HMSNO/DataFormat/HMSNO.Adata.Rdata")
load(file = "L:/Farr/HMSNO/DataFormat/HMSNO.Acon.Rdata")

attach(HMSNO.data)
attach(HMSNO.con)

#-------------------#
#-Load model output-#
#-------------------#

pattern <- "chain"

#Retrospective analysis
files <- list.files(path = "L:/Farr/HMSNO/DataAnalysis/Ignorance/Chains", pattern = pattern, full.names = TRUE)

nc <- 5

for(i in 1:nc){
  load(files[i])
}

out.retro <- mcmc.list(mget(ls()[grep(pattern, ls())]))

rm(list = ls()[grep(pattern, ls())])

# #BPVA: rcp 4.5
# files <- list.files(path = "L:/Farr/HMSNO/DataAnalysis/Annual_rcp45/Final_chains", pattern = pattern, full.names = TRUE)
# 
# for(i in 1:nc){
#   load(files[i])
# }
# 
# out.rcp45 <- mcmc.list(mget(ls()[grep(pattern, ls())]))
# 
# rm(list = ls()[grep(pattern, ls())])
# 
# #BPVA: rcp 8.5
# files <- list.files(path = "L:/Farr/HMSNO/DataAnalysis/Annual_rcp85/Final_chains", pattern = pattern, full.names = TRUE)
# 
# for(i in 1:nc){
#   load(files[i])
# }
# 
# out.rcp85 <-  mcmc.list(mget(ls()[grep(pattern, ls())]))
# 
# rm(list = ls()[grep(pattern, ls())])

#BPVA: business as usual
files <- list.files(path = "L:/Farr/HMSNO/DataAnalysis/BAU/chains", pattern = pattern, full.names = TRUE)

for(i in 1:nc){
  load(files[i])
}

out.bpva <- mcmc.list(mget(ls()[grep(pattern, ls())]))

rm(list = ls()[grep(pattern, ls())])


#------------#
#-Parameters-#
#------------#

params <- attr(out.retro[[1]], "dimnames")[[2]][!grepl("Npark|park.surv|pop.surv|park.gain|pop.gain", attr(out.retro[[1]], "dimnames")[[2]])]

#-------------#
#-Convergence-#
#-------------#

Rhat <- gelman.diag(out.retro[c(1:nc)][,params])
if(all(Rhat[[1]][,1] < 1.2)){
  print("Converged")
}else{
  tmp <- as.numeric(which(Rhat[[1]][,1] > 1.2))
  print("Not converged")
  print(params[tmp])
  traceplot(out.retro[c(1:nc)][,params[tmp]])
}

# Rhat <- gelman.diag(out.rcp45[c(1:nc)][,params])
# if(all(Rhat[[1]][,1] < 1.2)){
#   print("Converged")
# }else{
#   tmp <- as.numeric(which(Rhat[[1]][,1] > 1.2))
#   print("Not converged")
#   print(params[tmp])
#   traceplot(out.rcp45[c(1:nc)][,params[tmp]])
# }
# 
# Rhat <- gelman.diag(out.rcp85[c(1:nc)][,params])
# if(all(Rhat[[1]][,1] < 1.2)){
#   print("Converged")
# }else{
#   tmp <- as.numeric(which(Rhat[[1]][,1] > 1.2))
#   print("Not converged")
#   print(params[tmp])
#   traceplot(out.rcp85[c(1:nc)][,params[tmp]])
# }

Rhat <- gelman.diag(out.bpva[c(1:nc)][,params])
if(all(Rhat[[1]][,1] < 1.2)){
  print("Converged")
}else{
  tmp <- as.numeric(which(Rhat[[1]][,1] > 1.2))
  print("Not converged")
  print(params[tmp])
  traceplot(out.bpva[c(1:nc)][,params[tmp]])
}

#-----------#
#-Extaction-#
#-----------#

for(i in 1:6){
  params <- c(params, paste0("park.surv[", i, "]"))
}

for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    for(t in nstart[i]:nend[i]){
      params <- c(params, paste0("Npark[", t, ", ", i, ", ", s, "]"))
    }
    params <- c(params, paste0("pop.surv[", i, ", ", s, "]"))
  }
}
       
# names <- factor(c("Peters's", "Bay", "Harvey's", "White-bellied", "Blue",
#                   "Black-fronted", "Ogilby's",  "Abbott's", "Yellow-backed", "Community"),
#           levels =  c("Abbott's", "Bay", "Black-fronted", "Blue",
#           "Harvey's", "Ogilby's", "Peters's","White-bellied", 
#           "Yellow-backed", "Community"))

spec.names <- c("C. callipygus", "C. dorsalis", "C. harveyi", "C. leucogaster", 
                  "P. monticola", "N. moschatus", "C. nigrifrons", "C. ogilbyi",
                  "T. scriptus", "C. spadix", "T. spekii", "C. silvicultor")

park.names <- c("UDZ", "VNP", "NFNP", "BIF", "NNNP", "KRP")

year.names <- function(x){
  return(format(as.Date("2008", format = "%Y") + lubridate::years(x), "%Y"))
}

ni <- nc * length(out.retro[[1]][,1])

#-----------------------#
#-Retrospecitve results-#
#-----------------------#

N <- array(NA, dim = c(nspecs,max(parkE),11,ni))
park.surv <- array(NA, dim = c(max(parkE),ni))
pop.surv <- array(NA, dim = c(nspecs,max(parkE),ni))

for(i in 1:6){
  park.surv[i,] <- unlist(out.retro[c(1:nc)][,paste("park.surv[", i, "]", sep = "")])
}

for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    pop.surv[s,i,] <- unlist(out.retro[c(1:nc)][,paste("pop.surv[", i, ", ", s, "]", sep = "")])
    for(t in nstart[i]:nend[i]){
      N[s,i,t,] <- unlist(out.retro[c(1:nc)][,paste("Npark[", t, ", ", i, ", ", s, "]", sep = "")])
    }
  }
}

meanN <- apply(N, MARGIN = c(1,2,3), mean, na.rm = TRUE)
N97.5 <- apply(N, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
N2.5 <- apply(N, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)

meanpark.surv <- apply(park.surv, MARGIN = 1, mean, na.rm = TRUE)
park.surv97.5 <- apply(park.surv, MARGIN = 1, quantile, probs = 0.975, na.rm = TRUE)
park.surv75 <- apply(park.surv, MARGIN = 1, quantile, probs = 0.75, na.rm = TRUE)
park.surv25 <- apply(park.surv, MARGIN = 1, quantile, probs = 0.25, na.rm = TRUE)
park.surv2.5 <- apply(park.surv, MARGIN = 1, quantile, probs = 0.025, na.rm = TRUE)

meanpop.surv <- apply(pop.surv, MARGIN = c(1,2), mean, na.rm = TRUE)
pop.surv97.5 <- apply(pop.surv, MARGIN = c(1,2), quantile, probs = 0.975, na.rm = TRUE)
pop.surv75 <- apply(pop.surv, MARGIN = c(1,2), quantile, probs = 0.75, na.rm = TRUE)
pop.surv25 <- apply(pop.surv, MARGIN = c(1,2), quantile, probs = 0.25, na.rm = TRUE)
pop.surv2.5 <- apply(pop.surv, MARGIN = c(1,2), quantile, probs = 0.025, na.rm = TRUE)

data <- reshape2::melt(meanN, varnames = c("species", "parkID", "yearID"), value.name = "abundance")
data$species <- factor(data$species, labels = spec.names)
data$park <- factor(data$park, labels = park.names)
data$year <- as.numeric(year.names(data$yearID))
data$`97.5%` <- reshape2::melt(N97.5, varnames = c("species", "parkID", "yearID"), value.name = "97.5%")$`97.5%`
data$`2.5%` <- reshape2::melt(N2.5, varnames = c("species", "parkID", "yearID"), value.name = "2.5%")$`2.5%`
data$effort <- nsite[data$parkID]

data <- data %>%
  mutate(abundance_site = abundance/effort,
         `97.5%_site` = `97.5%`/effort,
         `2.5%_site` = `2.5%`/effort)

data$park <- factor(data$park, levels = c("KRP", "NNNP", "UDZ", "BIF", "NFNP", "VNP"))


parkdata <- data.frame(reshape2::melt(meanpark.surv, value.name = "mean"),
                       reshape2::melt(park.surv2.5, value.name = "2.5"),
                       reshape2::melt(park.surv25, value.name = "25"),
                       reshape2::melt(park.surv75, value.name = "75"),
                       reshape2::melt(park.surv97.5, value.name = "97.5"))
parkdata$park <- factor(park.names)
colnames(parkdata)[2:5] <- c("l2.5", "l25", "u75", "u97.5") 

popdata <- reshape2::melt(meanpop.surv, varnames = c("species", "parkID"), value.name = "mean")
popdata$`l2.5` <- reshape2::melt(pop.surv2.5, varnames = c("species", "parkID"), value.name = "2.5%")$`2.5`
popdata$`l25` <- reshape2::melt(pop.surv25, varnames = c("species", "parkID"), value.name = "25%")$`25`
popdata$`u75` <- reshape2::melt(pop.surv75, varnames = c("species", "parkID"), value.name = "75%")$`75`
popdata$`u97.5` <- reshape2::melt(pop.surv97.5, varnames = c("species", "parkID"), value.name = "97.5%")$`97.5`

popdata$species <- factor(popdata$species, labels = spec.names)
popdata$park <- factor(popdata$parkID, labels = park.names)
popdata <- popdata[complete.cases(popdata),]
popdata <- popdata[,-2]

#Summary statistics
measurements <- c("Initial", "Final", "PopGrowth")

#Extract summary statistics
Data.pop <- array(NA, dim = c(nspecs,max(parkE),length(measurements),ni))
for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    Data.pop[s,i,1,] <- N[s,i,nstart[i],]/nsite[i]
    Data.pop[s,i,2,] <- N[s,i,nend[i],]/nsite[i]
    tmp2 <- 1
    for(t in (nstart[i]+1):nend[i]){
      tmp1 <- N[s,i,t,]/N[s,i,t-1,]
      tmp2 <- tmp1 * tmp2
    }
    Data.pop[s,i,3,] <- tmp2^(1/(nend[i]-nstart[i]))
  }
}

#Calculate quantiles from MCMC
Data.mean <- apply(Data.pop, MARGIN = c(1,2,3), mean, na.rm = TRUE)
Data.97.5 <- apply(Data.pop, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
Data.2.5 <- apply(Data.pop, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)


#Format summary statistics into dataframe
sum.df <- full_join(full_join(dcast(melt(Data.mean, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements), 
                           dcast(melt(Data.2.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
                           by = c("species", "parkID")),
                 dcast(melt(Data.97.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
                 by = c("species", "parkID"))
colnames(sum.df) <- c("species", "park", paste0(measurements, rep(c("_mean", "_lower", "_upper"), each = length(measurements))))
sum.df$species <- factor(sum.df$species, labels = spec.names)
sum.df$park <- factor(sum.df$park, labels = park.names)
sum.df <- sum.df[complete.cases(sum.df),]
sum.df <- sum.df[,c("species", "park", paste0(rep(measurements, each = 3), c("_mean", "_lower", "_upper")))]

#Fill in missing values
for(k in 1:dim(sum.df)[1]){
  if(is.infinite(sum.df[k,"PopGrowth_mean"])){
    i <- as.numeric(sum.df[k,"species"])
    r <- as.numeric(sum.df[k,"park"])
    sum.df[k,"PopGrowth_mean"] <- mean(Data.pop[i,r,3,!is.infinite(Data.pop[i,r,3,])], na.rm = TRUE)
  }
  if(is.infinite(sum.df[k,"PopGrowth_lower"])){
    i <- as.numeric(sum.df[k,"species"])
    r <- as.numeric(sum.df[k,"park"])
    sum.df[k,"PopGrowth_lower"] <- quantile(Data.pop[i,r,3,!is.infinite(Data.pop[i,r,3,])], probs = 0.025, na.rm = TRUE)
  }
  if(is.infinite(sum.df[k,"PopGrowth_upper"])){
    i <- as.numeric(sum.df[k,"species"])
    r <- as.numeric(sum.df[k,"park"])
    sum.df[k,"PopGrowth_upper"] <- quantile(Data.pop[i,r,3,!is.infinite(Data.pop[i,r,3,])], probs = 0.975, na.rm = TRUE)
  }
}

#Initial abundances minimum
sum.df[which.min(sum.df$Initial_mean), c(1:5)]

#Initial abundances maximum
sum.df[which.max(sum.df$Initial_mean), c(1:5)]

#Final abundances minimum
sum.df[which.min(sum.df$Final_mean), c(1:2, 6:8)]

#Final abundances maximum
sum.df[which.max(sum.df$Final_mean), c(1:2, 6:8)]

#Population growth rate increasing
# quantile(Data.pop[9,1,3,!is.infinite(Data.pop[9,1,3,])], probs = c(0.025, 0.975))
sum(sum.df$PopGrowth_mean > 1)

sum(sum.df$PopGrowth_mean < 0.8)

#Population growth rate minimum
sum.df[which.min(sum.df$PopGrowth_mean), c(1:2, 9:11)]

#Population growth rate maximum
sum.df[which.max(sum.df$PopGrowth_mean), c(1:2, 9:11)]

#Population growth rate decreasing
sum.df[which(sum.df$PopGrowth_upper<1), c(1:2, 9:11)]

#NNNP
sum.df %>% filter(park == "NNNP")

sum.df %>% filter(park == "VNP")

sum.df %>% filter(park == "UDZ")

sum.df %>% filter(park == "BIF")

sum.df %>% filter(park == "NFNP")

sum.df %>% filter(park == "KRP")

#Save table
sum.df %>% arrange(park) %>% select(park, species, PopGrowth_mean, PopGrowth_lower, PopGrowth_upper) %>%
  mutate(PopGrowth_mean = round(PopGrowth_mean, digits = 2),
         PopGrowth_lower = round(PopGrowth_lower, digits = 2),
         PopGrowth_upper = round(PopGrowth_upper, digits = 2)) %>%
  write.csv(., file = "./SupportingInformation/AppS3TableS1.csv")


#---------------------------------#
#-Retrospecitve covariate effects-#
#---------------------------------#

output.retro <- summary(out.retro[c(1:nc)][,params])

# sum.fun <- function(x) {
#   c(mean(x), quantile(x, probs = c(0.025, 0.25, 0.75, 0.975)))
# }
# 
# omega0.logit <- do.call(rbind, out.retro[c(1:nc)][,params[c(grep("omega0", params)[c(7,9,11,12)], grep("mu.o0", params))]])
# omega0.logit[,5] <- logit(omega0.logit[,5])
# 
# apparent.survival <-  data.frame(species = c("nigrifrons", "scriptus", "spekii", "silvicultor", "community"),
#                                  t(apply(plogis(omega0.logit +
#                                    do.call(rbind, out.retro[c(1:nc)][,params[c(grep("omega1", params)[c(7,9,11,12)], grep("mu.o1", params))]])),
#                                    MARGIN = 2, sum.fun)),
#                                  disturbance = "high")
#   
# colnames(apparent.survival)[2:6] <- c("mean", "l2.5", "l25", "u75", "u97.5")

omega0 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(plogis(output.retro$statistics[grep("omega0", params),"Mean"]), output.retro$statistics["mu.o0","Mean"]),
                     l2.5 = c(plogis(output.retro$quantiles[grep("omega0", params),"2.5%"]), output.retro$quantiles["mu.o0","2.5%"]),
                     l25 = c(plogis(output.retro$quantiles[grep("omega0", params),"25%"]), output.retro$quantiles["mu.o0","25%"]),
                     u75 = c(plogis(output.retro$quantiles[grep("omega0", params),"75%"]), output.retro$quantiles["mu.o0","75%"]),
                     u97.5 = c(plogis(output.retro$quantiles[grep("omega0", params),"97.5%"]), output.retro$quantiles["mu.o0","97.5%"]))

# omega0$disturbance <- "low"
# 
# apparent.survival <- rbind(omega0, apparent.survival)
# 
# omega1 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("omega1", params),"Mean"], output.retro$statistics["mu.o1","Mean"]),
#                      l2.5 = c(output.retro$quantiles[grep("omega1", params),"2.5%"], output.retro$quantiles["mu.o1","2.5%"]),
#                      l25 = c(output.retro$quantiles[grep("omega1", params),"25%"], output.retro$quantiles["mu.o1","25%"]),
#                      u75 = c(output.retro$quantiles[grep("omega1", params),"75%"], output.retro$quantiles["mu.o1","75%"]),
#                      u97.5 = c(output.retro$quantiles[grep("omega1", params),"97.5%"], output.retro$quantiles["mu.o1","97.5%"]))
# 
# 
# apparent.survival <- rbind(omega0, apparent.survival)
# 
# apparent.survival$species <- factor(apparent.survival$species, levels = c("callipygus", "dorsalis", "harveyi", "leucogaster", 
#                                                                           "monticola", "moschatus",  "ogilbyi", "spadix", 
#                                                                           "nigrifrons", "scriptus", "spekii", "silvicultor",
#                                                                           "community"))
# 
# apparent.survival$disturbance <- factor(apparent.survival$disturbance, levels = c("low", "high"))
# omega0.logit <- do.call(rbind, out.retro[c(1:nc)][,params[c(grep("omega0", params)[c(7,9,11,12)], grep("mu.o0", params))]])
# omega0.logit[,5] <- logit(omega0.logit[,5])


popdata <- full_join(popdata, omega0 %>% mutate(park = "TEAM"))
popdata <- full_join(popdata, parkdata %>% mutate(species = "community"))
popdata$park <- factor(popdata$park, levels = c("KRP", "NNNP", "UDZ", "BIF", "NFNP", "VNP", "TEAM"))
popdata$species <- factor(popdata$species, levels = c(spec.names, "community"))

# omega2 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("omega2", params),"Mean"], output.retro$statistics["mu.o2","Mean"]),
#                      l2.5 = c(output.retro$quantiles[grep("omega2", params),"2.5%"], output.retro$quantiles["mu.o2","2.5%"]),
#                      l25 = c(output.retro$quantiles[grep("omega2", params),"25%"], output.retro$quantiles["mu.o2","25%"]),
#                      u75 = c(output.retro$quantiles[grep("omega2", params),"75%"], output.retro$quantiles["mu.o2","75%"]),
#                      u97.5 = c(output.retro$quantiles[grep("omega2", params),"97.5%"], output.retro$quantiles["mu.o2","97.5%"]))
# 
# omega3 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("omega3", params),"Mean"], output.retro$statistics["mu.o3","Mean"]),
#                      l2.5 = c(output.retro$quantiles[grep("omega3", params),"2.5%"], output.retro$quantiles["mu.o3","2.5%"]),
#                      l25 = c(output.retro$quantiles[grep("omega3", params),"25%"], output.retro$quantiles["mu.o3","25%"]),
#                      u75 = c(output.retro$quantiles[grep("omega3", params),"75%"], output.retro$quantiles["mu.o3","75%"]),
#                      u97.5 = c(output.retro$quantiles[grep("omega3", params),"97.5%"], output.retro$quantiles["mu.o3","97.5%"]))
# 
# 
# 
# omega4 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("omega4", params),"Mean"], output.retro$statistics["mu.o4","Mean"]),
#                      l2.5 = c(output.retro$quantiles[grep("omega4", params),"2.5%"], output.retro$quantiles["mu.o4","2.5%"]),
#                      l25 = c(output.retro$quantiles[grep("omega4", params),"25%"], output.retro$quantiles["mu.o4","25%"]),
#                      u75 = c(output.retro$quantiles[grep("omega4", params),"75%"], output.retro$quantiles["mu.o4","75%"]),
#                      u97.5 = c(output.retro$quantiles[grep("omega4", params),"97.5%"], output.retro$quantiles["mu.o4","97.5%"]))


#-Gains-#

# gamma0.log <- do.call(rbind, out.retro[c(1:nc)][,params[c(grep("gamma0", params)[c(7,9,11,12)], grep("mu.g0", params))]])
# 
# gains <- data.frame(species = c("nigrifrons", "scriptus", "spekii", "silvicultor", "community"),
#                                  t(apply(exp(gamma0.log +
#                                                   do.call(rbind, out.retro[c(1:nc)][,params[c(grep("gamma1", params)[c(7,9,11,12)], grep("mu.g1", params))]])),
#                                          MARGIN = 2, sum.fun)),
#                                  disturbance = "high")
# 
# colnames(gains)[2:6] <- c("mean", "l2.5", "l25", "u75", "u97.5")

gamma0 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(exp(output.retro$statistics[grep("gamma0", params),"Mean"]), exp(output.retro$statistics["mu.g0","Mean"])),
                     l2.5 = c(exp(output.retro$quantiles[grep("gamma0", params),"2.5%"]), exp(output.retro$quantiles["mu.g0","2.5%"])),
                     l25 = c(exp(output.retro$quantiles[grep("gamma0", params),"25%"]), exp(output.retro$quantiles["mu.g0","25%"])),
                     u75 = c(exp(output.retro$quantiles[grep("gamma0", params),"75%"]), exp(output.retro$quantiles["mu.g0","75%"])),
                     u97.5 = c(exp(output.retro$quantiles[grep("gamma0", params),"97.5%"]), exp(output.retro$quantiles["mu.g0","97.5%"])))

# gamma0$species <- factor(gamma0$species, levels = c("callipygus", "dorsalis", "harveyi", "leucogaster", 
#                                                                           "monticola", "moschatus",  "ogilbyi", "spadix", 
#                                                                           "nigrifrons", "scriptus", "spekii", "silvicultor",
#                                                                           "community"))
# 
# 
# gamma0$disturbance <- "low"
# 
# gains <- rbind(gamma0, gains)

# gamma1 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("gamma1", params),"Mean"], output.retro$statistics["mu.g1","Mean"]),
#                      l2.5 = c(output.retro$quantiles[grep("gamma1", params),"2.5%"], output.retro$quantiles["mu.g1","2.5%"]),
#                      l25 = c(output.retro$quantiles[grep("gamma1", params),"25%"], output.retro$quantiles["mu.g1","25%"]),
#                      u75 = c(output.retro$quantiles[grep("gamma1", params),"75%"], output.retro$quantiles["mu.g1","75%"]),
#                      u97.5 = c(output.retro$quantiles[grep("gamma1", params),"97.5%"], output.retro$quantiles["mu.g1","97.5%"]))
# 
# 
# beta0 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(plogis(output.retro$statistics[grep("beta0", params),"Mean"]), output.retro$statistics["mu.b0","Mean"]),
#                     l2.5 = c(plogis(output.retro$quantiles[grep("beta0", params),"2.5%"]), output.retro$quantiles["mu.b0","2.5%"]),
#                     l25 = c(plogis(output.retro$quantiles[grep("beta0", params),"25%"]), output.retro$quantiles["mu.b0","25%"]),
#                     u75 = c(plogis(output.retro$quantiles[grep("beta0", params),"75%"]), output.retro$quantiles["mu.b0","75%"]),
#                     u97.5 = c(plogis(output.retro$quantiles[grep("beta0", params),"97.5%"]), output.retro$quantiles["mu.b0","97.5%"]))
# 
# beta1 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("beta1", params),"Mean"], output.retro$statistics["mu.b1","Mean"]),
#                     l2.5 = c(output.retro$quantiles[grep("beta1", params),"2.5%"], output.retro$quantiles["mu.b1","2.5%"]),
#                     l25 = c(output.retro$quantiles[grep("beta1", params),"25%"], output.retro$quantiles["mu.b1","25%"]),
#                     u75 = c(output.retro$quantiles[grep("beta1", params),"75%"], output.retro$quantiles["mu.b1","75%"]),
#                     u97.5 = c(output.retro$quantiles[grep("beta1", params),"97.5%"], output.retro$quantiles["mu.b1","97.5%"]))
#
# beta2 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("beta2", params),"Mean"], output.retro$statistics["mu.b2","Mean"]),
#                     l2.5 = c(output.retro$quantiles[grep("beta2", params),"2.5%"], output.retro$quantiles["mu.b2","2.5%"]),
#                     l25 = c(output.retro$quantiles[grep("beta2", params),"25%"], output.retro$quantiles["mu.b2","25%"]),
#                     u75 = c(output.retro$quantiles[grep("beta2", params),"75%"], output.retro$quantiles["mu.b2","75%"]),
#                     u97.5 = c(output.retro$quantiles[grep("beta2", params),"97.5%"], output.retro$quantiles["mu.b2","97.5%"]))
# 
# beta3 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output.retro$statistics[grep("beta3", params),"Mean"], output.retro$statistics["mu.b3","Mean"]),
#                     l2.5 = c(output.retro$quantiles[grep("beta3", params),"2.5%"], output.retro$quantiles["mu.b3","2.5%"]),
#                     l25 = c(output.retro$quantiles[grep("beta3", params),"25%"], output.retro$quantiles["mu.b3","25%"]),
#                     u75 = c(output.retro$quantiles[grep("beta3", params),"75%"], output.retro$quantiles["mu.b3","75%"]),
#                     u97.5 = c(output.retro$quantiles[grep("beta3", params),"97.5%"], output.retro$quantiles["mu.b3","97.5%"]))

#FIX order of parks & sig figs
data %>% drop_na(.) %>% select(species, park, year, abundance_site, `2.5%_site`, `97.5%_site`) %>%
  arrange(park, species, year) %>%
  mutate(abundance_site = round(abundance_site, digits = 2),
         `2.5%_site` = round(`2.5%_site`, digits = 2),
         `97.5%_site` = round(`97.5%_site`, digits = 2),
         abundance = paste0(abundance_site, " (CI ",  `2.5%_site`, ", ", `97.5%_site`,  ")" )) %>%
  select(park, species, year, abundance) %>%
  write.csv(., file = "./SupportingInformation/AppS3TableS2.csv")

#Table S3
full_join(omega0 %>% mutate(mean = format(round(mean, digits=2), nsmall = 2),
                            l2.5 = format(round(l2.5, digits=2), nsmall = 2),
                            u97.5 = format(round(u97.5, digits=2), nsmall = 2), 
                            `Annual apparent survival` = paste0(mean, " (CI ",  l2.5, ", ", u97.5,  ")" )) %>% 
            select(species, `Annual apparent survival`), 
          gamma0 %>% mutate(mean = format(round(mean, digits=2), nsmall = 2),
                            l2.5 = format(round(l2.5, digits=2), nsmall = 2),
                            u97.5 = format(round(u97.5, digits=2), nsmall = 2),
                            Gains = paste0(mean, " (CI ",  l2.5, ", ", u97.5,  ")" )) %>% 
            select(species, Gains), 
          by = "species") %>%
  write.csv(., file = "./SupportingInformation/AppS3TableS3.csv")

#Table S4
popdata %>% filter(!(park == "TEAM") & species == "community") %>%
  mutate(mean = format(round(mean, digits=2), nsmall = 2),
         l2.5 = format(round(l2.5, digits=2), nsmall = 2),
         u97.5 = format(round(u97.5, digits=2), nsmall = 2), 
         `Annual apparent survival` = paste0(mean, " (CI ",  l2.5, ", ", u97.5,  ")" )) %>% 
  select(park, `Annual apparent survival`) %>%
  write.csv(., file = "./SupportingInformation/AppS3TableS4.csv")

#---------------------#
#-Retrospecitve plots-#
#---------------------#

#TEAM community mean
popdata[popdata$species == "community" & popdata$park == "TEAM",]

#Range of species survival
sd.o0 <- mean(1/sqrt(unlist(out.retro[c(1:nc)][,"tau.o0"])))
sd.o0; quantile(1/sqrt(unlist(out.retro[c(1:nc)][,"tau.o0"])), c(0.025, 0.975))

plogis(logit(output.retro$statistics["mu.o0","Mean"]) + sd.o0) -
plogis(logit(output.retro$statistics["mu.o0","Mean"]) - sd.o0)

#Range of park survival
sd.eps.o <- mean(1/sqrt(unlist(out.retro[c(1:nc)][,"tau.eps.o"])))
sd.eps.o; quantile(1/sqrt(unlist(out.retro[c(1:nc)][,"tau.eps.o"])), c(0.025, 0.975))

plogis(logit(output.retro$statistics["mu.o0","Mean"]) + sd.eps.o) -
  plogis(logit(output.retro$statistics["mu.o0","Mean"]) - sd.eps.o)


#Species vs park magnitude
data.frame(species = factor(spec.names), 
           mean = output.retro$statistics[grep("omega0", params),"Mean"] - logit(output.retro$statistics["mu.o0","Mean"]),
           l2.5 = output.retro$quantiles[grep("omega0", params),"2.5%"] - logit(output.retro$quantiles["mu.o0","2.5%"]),
           l25 = output.retro$quantiles[grep("omega0", params),"25%"] - logit(output.retro$quantiles["mu.o0","25%"]),
           u75 = output.retro$quantiles[grep("omega0", params),"75%"] - logit(output.retro$quantiles["mu.o0","75%"]),
           u97.5 = output.retro$quantiles[grep("omega0", params),"97.5%"] - logit(output.retro$quantiles["mu.o0","97.5%"]))

data.frame(park = factor(park.names),
           mean = output.retro$statistics[grep("eps.o", params),"Mean"][-7],
           l2.5 = output.retro$quantiles[grep("eps.o", params),"2.5%"][-7],
           l25 = output.retro$quantiles[grep("eps.o", params),"25%"][-7],
           u75 = output.retro$quantiles[grep("eps.o", params),"75%"][-7],
           u97.5 = output.retro$quantiles[grep("eps.o", params),"97.5%"][-7])

#Maximum species survival
max(popdata[popdata$species != "community",'mean'])

popdata[popdata$species != "community",]

popdata[popdata$species == "community" & popdata$park != "TEAM",]

#TEAM community gains
gamma0[gamma0$species == 'community',]


pos <- position_dodge(0.75, preserve = "total")


Figure1 <-  ggplot(data %>% drop_na(.), aes(x = year, y = abundance_site)) +
  geom_point(col = "grey") +
  geom_line(aes(col = species)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ park, scales = 'free', ncol = 2, dir = "v") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
  # guides(colour = guide_legend(nrow = 1)) +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic"),
    # legend.key.width = unit(0.5, "lines"),
    # legend.key.height = unit(5, "pt")
    ) +
  labs(y = "Density (abundance/site)\n", x = "Year")


tiff(file = "~/HMSNO/Figure1.tiff", res = 600, width = 6, height = 6, units = "in")
Figure1
dev.off()

# Figure2 <- ggplot(data %>% 
#                     filter(species %in% c("nigrifrons", "scriptus", "spekii", "silvicultor")) %>%
#                     mutate(y_max = ifelse(species == "nigrifrons", 3.75,
#                                          ifelse(species == "scriptus", 1.75,
#                                                 ifelse(species == "spekii", 0.3, 2.5)))),
#                   aes(x = year, y = abundance_site)) +
#   geom_line(aes(col = park, linetype = park), size = 1) + 
#   geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.1) +
#   facet_wrap(. ~ species, scales = 'free', 
#              labeller = labeller(species = c("nigrifrons" = "C. nigrifrons",
#                                              "scriptus" = "T. scriptus", 
#                                              "spekii" = "T. spekii", 
#                                              "silvicultor" = "C. silvicultor"))) +
#   geom_blank(aes(y = y_max)) +
#   scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
#   scale_fill_manual(values = c("KRP" = "#2166ac", "NNNP" = "#2166ac", "UDZ" = "#2166ac",
#                     "BIF" = "#b2182b", "NFNP" = "#b2182b", "VNP" = "#b2182b")) +
#   scale_color_manual(values = c("KRP" = "#2166ac", "NNNP" = "#2166ac", "UDZ" = "#2166ac",
#                                "BIF" = "#b2182b", "NFNP" = "#b2182b", "VNP" = "#b2182b")) +
#   scale_linetype_manual(values = c("KRP" = "solid", "NNNP" = "dashed", "UDZ" = "dotted",
#                                    "BIF" = "solid", "NFNP" = "dashed", "VNP" = "dotted")) +
#   theme_few() +
#   theme(
#     strip.text = element_text(face = "italic")) +
#   labs(y = "Density (abundance/site)\n", x = "Year", fill = "Park", linetype = "Park", color = "Park")
# 
# tiff(file = "~/HMSNO/Figure2.tiff", res = 600, width = 6, height = 5, units = "in")
# Figure2
# dev.off()

# Figure3A <- ggplot(apparent.survival, aes(x = species, color = disturbance)) + 
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 1, position = pos) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25, position = pos) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3, position = pos) +
#   scale_color_manual(values = c("low" = "#2166ac", "high" = "#b2182b"), labels = expression("< 35"/km^2, "> 300"/km^2)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5)) +
#   labs(y ="Apparent survival\n", x = expression(),  color = "Human density") +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# Figure2 <- ggplotGrob(ggplot(popdata %>% filter(park == "TEAM" | species == "community") %>% mutate(param = "a.surv", type = ifelse(park == "TEAM", "team", "park")) %>% 
#                     full_join(., gamma0 %>% mutate(param = "gain", park = "TEAM", type = "team")) %>% 
#                     mutate(park = factor(park, levels = c( "TEAM", "KRP", "NNNP", "UDZ", "BIF", "NFNP", "VNP"))), 
#                   aes(x = species, color = park)) +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, position = pos) +
#   scale_color_manual(values = c("black", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")) +
#   facet_wrap(. ~ param, scales = 'free_y', nrow = 2, strip.position = "left",
#              labeller = as_labeller(c(a.surv = "Apparent survival", gain = "Individuals gained \nper site"))) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5, face = "italic"),
#         legend.title = element_blank(),
#         strip.placement = "outside",
#         strip.background = element_blank(),
#         strip.text = element_text(size = 12)) +
#   labs(y = expression(), x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash"))
# 
# Figure2$layout[grepl("guide", Figure2$layout$name),'t'] <- 9

Fig2A <- popdata %>% filter(park == "TEAM" | species == "community") %>%
  mutate(xname = ifelse(park == "TEAM", as.character(species), as.character(park))) %>%
  mutate(xname = factor(xname, levels = c(spec.names, "community", park.names))) %>%
  ggplot(., aes(x = xname, color = park)) +
  geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
               stat = "identity", size = 0.75, width = 0.5) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "black")) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        # axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5, face = "italic"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        # legend.position = c(-5,5),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  labs(y = "Apparent survival", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")

Fig2B <- gamma0 %>%
  ggplot(., aes(x = species)) +
  geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
               stat = "identity", size = 0.75, width = 0.5, col = "black") +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5, face = "italic"),
        legend.title = element_blank(),
        # legend.position = c(-5,5),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  labs(y = "Individuals gained \nper site", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")

Fig2A <- ggplotGrob(Fig2A)
Fig2B <- ggplotGrob(Fig2B)

Fig2B$heights <- Fig2A$heights

Legend <- Fig2A$grobs[[15]]

Fig2A$grobs[[15]] <- nullGrob()

Fig2B$widths <- Fig2A$widths

Fig2A$widths[9] <- unit(0, "cm")

Fig2B$widths[9] <- sum(unit(5.4, "cm"), unit(2.5, "point"))

tiff(file = "~/HMSNO/Figure2.tiff", res = 600, width = 8, height = 5, units = "in")
grid.arrange(arrangeGrob(Fig2A, Fig2B), bottom = "\n \n", left = "")
grid.draw(textGrob("A", x = unit(0.5, "cm"), y = unit(12.25, "cm"),
                   gp=grid::gpar(fontsize=14, fontface = 2)))
grid.draw(textGrob("B", x = unit(0.5, "cm"), y = unit(6.75, "cm"),
                   gp=grid::gpar(fontsize=14, fontface = 2)))
pushViewport(viewport(x = unit(17, "cm"), y = unit(4.75, "cm")))
grid.draw(Legend)
dev.off()

# Figure2A <- ggplot(popdata %>% filter(park == "TEAM" | species == "community"), aes(x = species, color = park)) +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, position = pos) +
#   scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "black")) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5, face = "italic"),
#         legend.title = element_blank()) +
#   labs(y ="Apparent survival\n", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# Figure2B <- ggplot(gamma0, aes(x = species)) + 
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, color = "black") +  
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5, face = "italic"),
#         legend.position = "none") +
#   labs(y ="Gains (individuals/site)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# Fig2A <- ggplotGrob(Figure2A)
# Fig2B <- ggplotGrob(Figure2B)
# 
# Fig2B$widths <- Fig2A$widths
# 
# tiff(file = "~/HMSNO/Figure2.tiff", res = 600, width = 8, height = 5, units = "in")
# grid.arrange(arrangeGrob(Fig2A, Fig2B, ncol = 1, nrow = 2))
# dev.off()
# 
# Figure3A_alt1 <- ggplot(popdata %>% filter(park == "TEAM" | species == "community"), aes(x = species, color = park)) +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, position = pos) +
#   # scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "black")) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5)) +
#   labs(y ="Apparent survival\n", x = expression()) +
#   scale_x_continuous(breaks = xpos) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# tiff(file = "~/HMSNO/Figure3_alt1.tiff", res = 600, width = 11, height = 5, units = "in")
# Figure3A_alt1
# dev.off()
# 
# Figure3A_alt2 <- ggplot(popdata %>% filter(park != "TEAM"), aes(x = park, color = species)) +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, position = pos) +
#   scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "grey", "#b15928", "black")) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5)) +
#   labs(y ="Apparent survival\n", x = expression())
# 
# tiff(file = "~/HMSNO/Figure3A_alt2.tiff", res = 600, width = 8, height = 5, units = "in")
# Figure3A_alt2
# dev.off()


# omega2$species <- factor(omega2$species, levels = c("callipygus", "dorsalis", "harveyi", "leucogaster", 
#                                                     "monticola", "moschatus",  "ogilbyi", "spadix", 
#                                                     "nigrifrons", "scriptus", "spekii", "silvicultor",
#                                                     "community"))
# 
# Figure4 <- ggplot(omega2, aes(x = species)) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, alpha = 0) +  
#   coord_cartesian(ylim = c(-2.75, 1.25)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Effect of annual rainfall\n(logit-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# tiff(file = "~/HMSNO/Figure4.tiff", res = 600, width = 8, height = 5, units = "in")
# Figure4
# dev.off()
# 
# Figure4B <- ggplot(omega3, aes(x = species)) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, alpha = 0) +  
#   coord_cartesian(ylim = c(-1.25, 3.5)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Effect of elevation\n(logit-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# tiff(file = "~/HMSNO/Figure4B.tiff", res = 600, width = 8, height = 5, units = "in")
# Figure4B
# dev.off()
# 
# Figure4C <- ggplot(omega4, aes(x = species)) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, alpha = 0) +  
#   coord_cartesian(ylim = c(-1, 2)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Distance from edge\n(logit-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# tiff(file = "~/HMSNO/Figure4C.tiff", res = 600, width = 8, height = 5, units = "in")
# Figure4C
# dev.off()
# 
# # beta2$species <- factor(beta2$species, levels = c("callipygus", "dorsalis", "harveyi", "leucogaster", 
# #                                                     "monticola", "moschatus",  "ogilbyi", "spadix", 
# #                                                     "nigrifrons", "scriptus", "spekii", "silvicultor",
# #                                                     "community"))
# 
# Figure5 <- ggplot(beta2, aes(x = species)) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, alpha = 0) +  
#   coord_cartesian(ylim = c(-0.6, 0.6)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Distance from edge\n(log-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# tiff(file = "~/HMSNO/Figure5.tiff", res = 600, width = 8, height = 5, units = "in")
# Figure5
# dev.off()
# 
# # beta4$species <- factor(beta4$species, levels = c("callipygus", "dorsalis", "harveyi", "leucogaster", 
# #                                                   "monticola", "moschatus",  "ogilbyi", "spadix", 
# #                                                   "nigrifrons", "scriptus", "spekii", "silvicultor",
# #                                                   "community"))
# 
# Figure6 <- ggplot(beta3, aes(x = species)) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
#                stat = "identity", size = 0.75, width = 0.5, alpha = 0) + 
#   coord_cartesian(ylim = c(-1, 1)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Effect of elevation\n(log-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# tiff(file = "~/HMSNO/Figure6.tiff", res = 600, width = 8, height = 5, units = "in")
# Figure6
# dev.off()

# omega4$species <- factor(omega4$species, levels = c("callipygus", "dorsalis", "harveyi", "leucogaster", 
#                                                     "monticola", "moschatus",  "ogilbyi", "spadix", 
#                                                     "nigrifrons", "scriptus", "spekii", "silvicultor",
#                                                     "community"))
# 
# Figure7 <- ggplot(omega4) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   coord_cartesian(ylim = c(-2, 2)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Effect of elevation\n(logit-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# tiff(file = "~/HMSNO/Figure7.tiff", res = 600, width = 8, height = 5, units = "in")
# Figure7
# dev.off()


# FigO0 <- ggplot(omega0) + 
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Apparent survival", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# 
# FigO1 <- ggplot(omega1 %>% filter(species %in% c("nigrifrons", "scriptus", "spekii", "silvicultor", "community"))) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Density effect (logit-scale)", x = expression()) +
#   geom_vline(xintercept=(5 - 0.5), linetype = "dotdash")


# 
# FigO2 <- ggplot(omega2) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   coord_cartesian(ylim = c(-2, 2)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Rainfall", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# #-Gains-#
# 
# FigG0 <- ggplot(gamma0) + 
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Gains per site", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# 
# FigG1 <- ggplot(gamma1 %>% filter(species %in% c("nigrifrons", "scriptus", "spekii", "silvicultor", "community"))) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Effect of density (log-scale)", x = expression()) +
#   geom_vline(xintercept=(5 - 0.5), linetype = "dotdash")
# 
# Fig.Gains <- ggplot(gains, aes(x = species, color = disturbance)) + 
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3, position = pos) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25, position = pos) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3, position = pos) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5)) +
#   labs(y ="Apparent survival", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# hist(gamma0.log[,5])
# hist(do.call(rbind, out.retro[c(1:nc)][,params[grep("mu.g1", params)]]))
# hist(gamma0.log[,5] + do.call(rbind, out.retro[c(1:nc)][,params[grep("mu.g1", params)]]))
# hist(exp(gamma0.log[,5] + do.call(rbind, out.retro[c(1:nc)][,params[grep("mu.g1", params)]])))
# hist(exp(gamma0.log[,5]))
# 
# #-Initial abundance-#
# 
# FigB1 <- ggplot(beta1 %>% filter(species %in% c("nigrifrons", "scriptus", "spekii", "silvicultor", "community"))) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   coord_cartesian(ylim = c(-5, 5)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Density effect (log-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# 
# FigB2 <- ggplot(beta2) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   coord_cartesian(ylim = c(-0.75, 0.75)) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Edge effect (log-scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")
# 
# FigB3 <- ggplot(beta3 %>% filter(species %in% c("nigrifrons", "scriptus", "spekii", "silvicultor", "community"))) + 
#   geom_hline(yintercept = 0, col = "red") +
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.3) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#   labs(y ="Interaction: edge x density (log scale)", x = expression()) +
#   geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")

#------------------------------#
#-Quasi-extinction probability-#
#------------------------------#

# #-Rcp 4.5-#
# 
# N.rcp45 <- array(NA, dim = c(nspecs,max(parkE),22,ni))
# 
# for(s in 1:nspecs){
#   for(i in parkS[s]:parkE[s]){
#     for(t in nstart[i]:22){
#       N.rcp45[s,i,t,] <- unlist(out.rcp45[c(1:nc)][,paste("Npark[", t, ", ", i, ", ", s, "]", sep = "")])
#     }
#   }
# }
# 
# meanN.rcp45 <- apply(N.rcp45, MARGIN = c(1,2,3), mean, na.rm = TRUE)
# N.rcp45.97.5 <- apply(N.rcp45, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
# N.rcp45.2.5 <- apply(N.rcp45, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)
# 
# data.rcp45 <- reshape2::melt(meanN.rcp45, varnames = c("species", "parkID", "yearID"), value.name = "abundance")
# data.rcp45$species <- factor(data.rcp45$species, labels = spec.names)
# data.rcp45$park <- factor(data$park, labels = park.names)
# data.rcp45$year <- as.numeric(year.names(data.rcp45$yearID))
# 
# data.rcp45$`97.5%` <- reshape2::melt(N.rcp45.97.5, varnames = c("species", "parkID", "yearID"), value.name = "97.5%")$`97.5%`
# data.rcp45$`2.5%` <- reshape2::melt(N.rcp45.2.5, varnames = c("species", "parkID", "yearID"), value.name = "2.5%")$`2.5%`
# 
# 
# data.rcp45$effort <- nsite[data.rcp45$parkID]
# 
# data.rcp45 <- data.rcp45 %>%
#   mutate(abundance_site = abundance/effort,
#          `97.5%_site` = `97.5%`/effort,
#          `2.5%_site` = `2.5%`/effort)
# 
# 
# 
# #Summary statistics
# measurements <- c("Initial", "Final", "Projected", "Threshold", "QSE")
# 
# #Extract summary statistics
# Data.rcp45 <- array(NA, dim = c(nspecs,max(parkE),length(measurements),ni))
# for(s in 1:nspecs){
#   for(i in parkS[s]:parkE[s]){
#     Data.rcp45[s,i,1,] <- N.rcp45[s,i,nstart[i],]
#     Data.rcp45[s,i,2,] <- N.rcp45[s,i,nend[i],]
#     Data.rcp45[s,i,3,] <- N.rcp45[s,i,22,]
#     Data.rcp45[s,i,4,] <- N.rcp45[s,i,11,] * 0.2
#     Data.rcp45[s,i,5,] <- N.rcp45[s,i,22,] < (N.rcp45[s,i,11,] * 0.2)
#   }
# }
# 
# #Calculate quantiles from MCMC
# Rcp45.mean <- apply(Data.rcp45, MARGIN = c(1,2,3), mean, na.rm = TRUE)
# Rcp45.97.5 <- apply(Data.rcp45, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
# Rcp45.2.5 <- apply(Data.rcp45, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)
# 
# #Format summary statistics into dataframe
# qse45 <- full_join(full_join(dcast(melt(Rcp45.mean, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements), 
#                   dcast(melt(Rcp45.2.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
#                  by = c("species", "parkID")),
#                  dcast(melt(Rcp45.97.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
#                  by = c("species", "parkID"))
# colnames(qse45) <- c("species", "park", paste0(measurements, rep(c("_mean", "_lower", "_upper"), each = length(measurements))))
# qse45$species <- factor(qse45$species, labels = spec.names)
# qse45$park <- factor(qse45$park, labels = park.names)
# qse45 <- qse45[complete.cases(qse45),]
# qse45 <- qse45[,c("species", "park", paste0(rep(measurements, each = 3), c("_mean", "_lower", "_upper")))]
# 
# #Remove columns
# qse45 <- qse45[,-c(16,17)]
# # qse$Percent_change <- (qse[,"Final_mean"] - qse[,"Initial_mean"])/qse[,"Initial_mean"]
# # qse$Projected_change <- (qse[,"Projected_mean"] - qse[,"Final_mean"])/qse[,"Final_mean"]
# 
# #-Rcp 8.5-#
# 
# N.rcp85 <- array(NA, dim = c(nspecs,max(parkE),22,ni))
# 
# for(s in 1:nspecs){
#   for(i in parkS[s]:parkE[s]){
#     for(t in nstart[i]:22){
#       N.rcp85[s,i,t,] <- unlist(out.rcp85[c(1:nc)][,paste("Npark[", t, ", ", i, ", ", s, "]", sep = "")])
#     }
#   }
# }
# 
# meanN.rcp85 <- apply(N.rcp85, MARGIN = c(1,2,3), mean, na.rm = TRUE)
# N.rcp85.97.5 <- apply(N.rcp85, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
# N.rcp85.2.5 <- apply(N.rcp85, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)
# 
# data.rcp85 <- reshape2::melt(meanN.rcp85, varnames = c("species", "parkID", "yearID"), value.name = "abundance")
# data.rcp85$species <- factor(data.rcp85$species, labels = spec.names)
# data.rcp85$park <- factor(data$park, labels = park.names)
# data.rcp85$year <- as.numeric(year.names(data.rcp85$yearID))
# 
# high.rcp85 <- reshape2::melt(N.rcp85.97.5, varnames = c("species", "parkID", "yearID"), value.name = "97.5%")$`97.5%`
# low.rcp85 <- reshape2::melt(N.rcp85.2.5, varnames = c("species", "parkID", "yearID"), value.name = "2.5%")$`2.5%`
# data.rcp85$`97.5%` <- high.rcp85
# data.rcp85$`2.5%` <- low.rcp85
# 
# data.rcp85$effort <- nsite[data.rcp85$parkID]
# 
# data.rcp85 <- data.rcp85 %>%
#   mutate(abundance_site = abundance/effort,
#          `97.5%_site` = `97.5%`/effort,
#          `2.5%_site` = `2.5%`/effort)
# 
# 
# #Summary statistics
# measurements <- c("Initial", "Final", "Projected", "Threshold", "QSE")
# 
# #Extract summary statistics
# Data.rcp85 <- array(NA, dim = c(nspecs,max(parkE),length(measurements),ni))
# for(s in 1:nspecs){
#   for(i in parkS[s]:parkE[s]){
#     Data.rcp85[s,i,1,] <- N.rcp85[s,i,nstart[i],]
#     Data.rcp85[s,i,2,] <- N.rcp85[s,i,nend[i],]
#     Data.rcp85[s,i,3,] <- N.rcp85[s,i,22,]
#     Data.rcp85[s,i,4,] <- N.rcp85[s,i,11,] * 0.2
#     Data.rcp85[s,i,5,] <- N.rcp85[s,i,22,] < (N.rcp85[s,i,11,] * 0.2)
#   }
# }
# 
# #Calculate quantiles from MCMC
# Rcp85.mean <- apply(Data.rcp85, MARGIN = c(1,2,3), mean, na.rm = TRUE)
# Rcp85.97.5 <- apply(Data.rcp85, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
# Rcp85.2.5 <- apply(Data.rcp85, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)
# 
# #Format summary statistics into dataframe
# qse85 <- full_join(full_join(dcast(melt(Rcp85.mean, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements), 
#                              dcast(melt(Rcp85.2.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
#                              by = c("species", "parkID")),
#                    dcast(melt(Rcp85.97.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
#                    by = c("species", "parkID"))
# colnames(qse85) <- c("species", "park", paste0(measurements, rep(c("_mean", "_lower", "_upper"), each = length(measurements))))
# qse85$species <- factor(qse85$species, labels = spec.names)
# qse85$park <- factor(qse85$park, labels = park.names)
# qse85 <- qse85[complete.cases(qse85),]
# qse85 <- qse85[,c("species", "park", paste0(rep(measurements, each = 3), c("_mean", "_lower", "_upper")))]
# 
# #Remove columns
# qse85 <- qse85[,-c(16,17)]

#-BPVA-#

N.bpva <- array(NA, dim = c(nspecs,max(parkE),22,ni))

for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    for(t in nstart[i]:22){
      N.bpva[s,i,t,] <- unlist(out.bpva[c(1:nc)][,paste("Npark[", t, ", ", i, ", ", s, "]", sep = "")])
    }
  }
}

meanN.bpva <- apply(N.bpva, MARGIN = c(1,2,3), mean, na.rm = TRUE)
N.bpva.97.5 <- apply(N.bpva, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
N.bpva.2.5 <- apply(N.bpva, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)

data.bpva <- reshape2::melt(meanN.bpva, varnames = c("species", "parkID", "yearID"), value.name = "abundance")
data.bpva$species <- factor(data.bpva$species, labels = spec.names)
data.bpva$park <- factor(data.bpva$park, labels = park.names)
data.bpva$park <- factor(data.bpva$park, levels = c("KRP", "NNNP", "UDZ", "BIF", "NFNP", "VNP"))
data.bpva$year <- as.numeric(year.names(data.bpva$yearID))

data.bpva$`97.5%` <- reshape2::melt(N.bpva.97.5, varnames = c("species", "parkID", "yearID"), value.name = "97.5%")$`97.5%`
data.bpva$`2.5%` <- reshape2::melt(N.bpva.2.5, varnames = c("species", "parkID", "yearID"), value.name = "2.5%")$`2.5%`


data.bpva$effort <- nsite[data.bpva$parkID]

data.bpva <- data.bpva %>%
  mutate(abundance_site = abundance/effort,
         `97.5%_site` = `97.5%`/effort,
         `2.5%_site` = `2.5%`/effort)

data.bpva <- data.bpva %>% mutate(endyear = as.numeric(as.character(factor(parkID, labels = nend))))

#Summary statistics
measurements <- c("Initial", "Final", "Projected", "Threshold", "QSE")

#Extract summary statistics
Data.bpva <- array(NA, dim = c(nspecs,max(parkE),length(measurements),ni))
for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    Data.bpva[s,i,1,] <- N.bpva[s,i,nstart[i],]
    Data.bpva[s,i,2,] <- N.bpva[s,i,nend[i],]
    Data.bpva[s,i,3,] <- N.bpva[s,i,22,]
    Data.bpva[s,i,4,] <- N.bpva[s,i,11,] * 0.5
    Data.bpva[s,i,5,] <- N.bpva[s,i,22,] < (mean(N.bpva[s,i,11,]) * 0.5)
  }
}

#Calculate quantiles from MCMC
bpva.mean <- apply(Data.bpva, MARGIN = c(1,2,3), mean, na.rm = TRUE)
bpva.97.5 <- apply(Data.bpva, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
bpva.2.5 <- apply(Data.bpva, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)

#Format summary statistics into dataframe
qse <- full_join(full_join(dcast(melt(bpva.mean, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements), 
                             dcast(melt(bpva.2.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
                             by = c("species", "parkID")),
                   dcast(melt(bpva.97.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
                   by = c("species", "parkID"))
colnames(qse) <- c("species", "park", paste0(measurements, rep(c("_mean", "_lower", "_upper"), each = length(measurements))))
qse$species <- factor(qse$species, labels = spec.names)
qse$park <- factor(qse$park, labels = park.names)
qse$park <- factor(qse$park, levels = c("KRP", "NNNP", "UDZ", "BIF", "NFNP", "VNP"))
qse <- qse[complete.cases(qse),]
qse <- qse[,c("species", "park", paste0(rep(measurements, each = 3), c("_mean", "_lower", "_upper")))]

#Remove columns
qse <- qse[,-c(16,17)]
# qse$Percent_change <- (qse[,"Final_mean"] - qse[,"Initial_mean"])/qse[,"Initial_mean"]
# qse$Projected_change <- (qse[,"Projected_mean"] - qse[,"Final_mean"])/qse[,"Final_mean"]

qse %>%
  mutate(`Quasi-extinction probability` = format(round(QSE_mean, digits=2), nsmall = 2)) %>%
  arrange(park, species) %>%
  select(park, species, `Quasi-extinction probability`) %>%
  filter(!(park == "VNP"|park == "NFNP")) %>%
  write.csv(., file = "./SupportingInformation/AppS3TableS5.csv")

#--------------------#
#-Summary statistics-#
#--------------------#

# #Inital abundances by species
# print(qse45 %>% group_by(species) %>% select(Initial_mean, Initial_lower, Initial_upper),
#       n = 23)
# 
# print(qse45 %>% group_by(species) %>% summarize(sd(Initial_mean)),
#       n = 23)
# 
# #Inital abundances by park
# print(qse45 %>% group_by(park) %>% arrange(park) %>% select(Initial_mean, Initial_lower, Initial_upper),
#       n = 23)
# 
# print(qse45 %>% group_by(park) %>% arrange(park) %>% summarize(sd(Initial_mean)),
#       n = 23)
# 
# #Final abundances by species
# print(qse45 %>% group_by(species) %>% select(Final_mean, Final_lower, Final_upper),
#       n = 23)
# 
# #Final abundances by park
# print(qse45 %>% group_by(park) %>% arrange(park) %>% select(Final_mean, Final_lower, Final_upper),
#       n = 23)

# #Quasi-extinction by species
# print(qse45 %>% group_by(species) %>% select(QSE_mean),
#       n = 23)
# qse45 %>% group_by(species) %>% summarize(median(QSE_mean), max(QSE_mean))
# #Quasi-extinction by park
# qse45 %>% group_by(park) %>% summarize(median(QSE_mean), max(QSE_mean))
# 
# #Quasi-extinction by species
# print(qse85 %>% group_by(species) %>% select(QSE_mean),
#       n = 23)
# qse85 %>% group_by(species) %>% summarize(median(QSE_mean), max(QSE_mean))
# #Quasi-extinction by park
# qse85 %>% group_by(park) %>% summarize(median(QSE_mean), max(QSE_mean))
# 

#Quasi-extinction by species
print(qse %>% group_by(species) %>% select(park, QSE_mean),
      n = 23)
qse %>% group_by(species) %>% summarize(median(QSE_mean), max(QSE_mean))
#Quasi-extinction by park
qse %>% group_by(park) %>% summarize(median(QSE_mean), max(QSE_mean))


data.rcp45 %>% 
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(data = data.rcp45 %>% filter(yearID < 12), aes(col = park)) + 
  geom_ribbon(data = data.rcp45 %>% filter(yearID < 12), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.5) +
  geom_line(data = data.rcp45 %>% filter(yearID >= 11), aes(col = park), linetype = "dashed") + 
  geom_ribbon(data = data.rcp45 %>% filter(yearID >= 11), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.25) +
  facet_wrap(. ~ species) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  labs(y = "Density (abundance/site")

data.rcp85 %>% 
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(data = data.rcp85 %>% filter(yearID < 12), aes(col = park)) + 
  geom_ribbon(data = data.rcp85 %>% filter(yearID < 12), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.5) +
  geom_line(data = data.rcp85 %>% filter(yearID >= 11), aes(col = park), linetype = "dashed") + 
  geom_ribbon(data = data.rcp85 %>% filter(yearID >= 11), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.25) +
  facet_wrap(. ~ species) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  labs(y = "Density (abundance/site")

data.rcp45 %>% 
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(data = data.rcp45 %>% filter(yearID < 12), aes(col = species)) + 
  geom_ribbon(data = data.rcp45 %>% filter(yearID < 12), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  geom_line(data = data.rcp45 %>% filter(yearID >= 11), aes(col = species), linetype = "dashed") + 
  geom_ribbon(data = data.rcp45 %>% filter(yearID >= 11), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.25) +
  facet_wrap(. ~ park, scales = 'free_x') +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  labs(y = "Density (abundance/site")

data.rcp85 %>% 
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(data = data.rcp85 %>% filter(yearID < 12), aes(col = species)) + 
  geom_ribbon(data = data.rcp85 %>% filter(yearID < 12), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  geom_line(data = data.rcp85 %>% filter(yearID >= 11), aes(col = species), linetype = "dashed") + 
  geom_ribbon(data = data.rcp85 %>% filter(yearID >= 11), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.25) +
  facet_wrap(. ~ park, scales = 'free_x') +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  labs(y = "Density (abundance/site")



Figure3 <-  ggplot(data.bpva %>% filter(!(park == "VNP"|park == "NFNP")), aes(x = year, y = abundance_site)) +
  # geom_point(col = "grey") +
  geom_line(data = data.bpva %>% filter(!(park == "VNP"|park == "NFNP") & yearID <= endyear), aes(col = species)) + 
  geom_ribbon(data = data.bpva %>% filter(!(park == "VNP"|park == "NFNP") & yearID <= endyear), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  geom_line(data = data.bpva %>% filter(!(park == "VNP"|park == "NFNP") & yearID >= endyear), aes(col = species), linetype = "dashed") + 
  geom_ribbon(data = data.bpva %>% filter(!(park == "VNP"|park == "NFNP") & yearID >= endyear), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.25) +
  facet_wrap(. ~ park, scales = 'free', ncol = 2, dir = "v") +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic"),
    # legend.key.width = unit(0.5, "lines"),
    # legend.key.height = unit(5, "pt")
  ) +
  labs(y = "Density (abundance/site)\n", x = "Year")


ggplot(data.bpva, aes(x = year, y = abundance_site)) +
  # geom_point(col = "grey") +
  geom_line(data = data.bpva %>% filter(yearID <= endyear), aes(col = species)) + 
  geom_ribbon(data = data.bpva %>% filter(yearID <= endyear), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  # geom_line(data = data.bpva %>% filter(yearID > endyear), aes(col = species), linetype = "dashed") + 
  # geom_ribbon(data = data.bpva %>% filter(yearID > endyear), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.25) +
  facet_wrap(. ~ park, scales = 'free', ncol = 2, dir = "v") +
  scale_x_continuous(breaks = c(2010, 2015, 2020)) +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic"),
    # legend.key.width = unit(0.5, "lines"),
    # legend.key.height = unit(5, "pt")
  ) +
  labs(y = "Density (abundance/site)\n", x = "Year")

ggplot(data.bpva %>% drop_na(.) %>% filter(yearID <= endyear), aes(x = year, y = abundance_site)) +
  geom_point(col = "grey") +
  geom_line(aes(col = species)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ park, scales = 'free', ncol = 2, dir = "v") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
  # guides(colour = guide_legend(nrow = 1)) +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic"),
    # legend.key.width = unit(0.5, "lines"),
    # legend.key.height = unit(5, "pt")
  ) +
  labs(y = "Density (abundance/site)\n", x = "Year")


tiff(file = "~/HMSNO/Figure3_draft15.tiff", res = 600, width = 6, height = 6, units = "in")
Figure3
dev.off()

Figure1.rcp85 <-  ggplot(data.rcp85, aes(x = year, y = abundance_site)) +
  geom_line(data = data.rcp85 %>% filter(yearID < 12), aes(col = species)) + 
  geom_ribbon(data = data.rcp85 %>% filter(yearID < 12), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  geom_line(data = data.rcp85 %>% filter(yearID >= 11), aes(col = species), linetype = "dashed") + 
  geom_ribbon(data = data.rcp85 %>% filter(yearID >= 11), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.25) +
  facet_wrap(. ~ park, scales = 'free', ncol = 2, dir = "v") +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  labs(y = "Density (abundance/site)\n", x = "Year")

tiff(file = "~/HMSNO/Figure1_rcp85.tiff", res = 600, width = 6, height = 6, units = "in")
Figure1.rcp85
dev.off()


Figure2.rcp45 <- ggplot(data.rcp45 %>% 
                    filter(species %in% c("nigrifrons", "scriptus", "spekii", "silvicultor")) %>%
                    mutate(y_max = ifelse(species == "nigrifrons", 3.75,
                                          ifelse(species == "scriptus", 1.75,
                                                 ifelse(species == "spekii", 0.3, 2.5)))),
                  aes(x = year, y = abundance_site)) +
  geom_line(aes(col = park, linetype = park), size = 1) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.1) +
  facet_wrap(. ~ species, scales = 'free', 
             labeller = labeller(species = c("nigrifrons" = "C. nigrifrons",
                                             "scriptus" = "T. scriptus", 
                                             "spekii" = "T. spekii", 
                                             "silvicultor" = "C. silvicultor"))) +
  geom_blank(aes(y = y_max)) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  scale_fill_manual(values = c("KRP" = "#2166ac", "NNNP" = "#2166ac", "UDZ" = "#2166ac",
                               "BIF" = "#b2182b", "NFNP" = "#b2182b", "VNP" = "#b2182b")) +
  scale_color_manual(values = c("KRP" = "#2166ac", "NNNP" = "#2166ac", "UDZ" = "#2166ac",
                                "BIF" = "#b2182b", "NFNP" = "#b2182b", "VNP" = "#b2182b")) +
  scale_linetype_manual(values = c("KRP" = "solid", "NNNP" = "dashed", "UDZ" = "dotted",
                                   "BIF" = "solid", "NFNP" = "dashed", "VNP" = "dotted")) +
  theme_few() +
  theme(
    strip.text = element_text(face = "italic")) +
  labs(y = "Density (abundance/site)\n", x = "Year", fill = "Park", linetype = "Park", color = "Park")

tiff(file = "~/HMSNO/Figure2_rcp45.tiff", res = 600, width = 6, height = 5, units = "in")
Figure2.rcp45
dev.off()

Figure2.rcp85 <- ggplot(data.rcp85 %>% 
                          filter(species %in% c("nigrifrons", "scriptus", "spekii", "silvicultor")) %>%
                          mutate(y_max = ifelse(species == "nigrifrons", 3.75,
                                                ifelse(species == "scriptus", 1.75,
                                                       ifelse(species == "spekii", 0.3, 2.5)))),
                        aes(x = year, y = abundance_site)) +
  geom_line(aes(col = park, linetype = park), size = 1) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.1) +
  facet_wrap(. ~ species, scales = 'free', 
             labeller = labeller(species = c("nigrifrons" = "C. nigrifrons",
                                             "scriptus" = "T. scriptus", 
                                             "spekii" = "T. spekii", 
                                             "silvicultor" = "C. silvicultor"))) +
  geom_blank(aes(y = y_max)) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  scale_fill_manual(values = c("KRP" = "#2166ac", "NNNP" = "#2166ac", "UDZ" = "#2166ac",
                               "BIF" = "#b2182b", "NFNP" = "#b2182b", "VNP" = "#b2182b")) +
  scale_color_manual(values = c("KRP" = "#2166ac", "NNNP" = "#2166ac", "UDZ" = "#2166ac",
                                "BIF" = "#b2182b", "NFNP" = "#b2182b", "VNP" = "#b2182b")) +
  scale_linetype_manual(values = c("KRP" = "solid", "NNNP" = "dashed", "UDZ" = "dotted",
                                   "BIF" = "solid", "NFNP" = "dashed", "VNP" = "dotted")) +
  theme_few() +
  theme(
    strip.text = element_text(face = "italic")) +
  labs(y = "Density (abundance/site)\n", x = "Year", fill = "Park", linetype = "Park", color = "Park")

tiff(file = "~/HMSNO/Figure2_rcp85.tiff", res = 600, width = 6, height = 5, units = "in")
Figure2.rcp85
dev.off()


#------------------------#
#-Rainfall visualization-#
#------------------------#



#-----------------------------#
#-Bayesian variable selection-#
#-----------------------------#

z1 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output$statistics[grep("z1", params),"Mean"], output$statistics["psi1","Mean"]),
                     l2.5 = c(output$quantiles[grep("z1", params),"2.5%"], output$quantiles["psi1","2.5%"]),
                     l25 = c(output$quantiles[grep("z1", params),"25%"], output$quantiles["psi1","25%"]),
                     u75 = c(output$quantiles[grep("z1", params),"75%"], output$quantiles["psi1","75%"]),
                     u97.5 = c(output$quantiles[grep("z1", params),"97.5%"], output$quantiles["psi1","97.5%"]))

FigZ1 <- ggplot() + 
  geom_hline(yintercept = 0.5, col = "red") +
  geom_point(data = z1, aes(x = species, y = mean),
                size = 3, shape = 1) +
  geom_errorbar(data = subset(z1, species == "community"), aes(x = species, ymin = mean, ymax = mean),
                width = 0.275) +
  geom_errorbar(data = subset(z1, species == "community"), aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(data = subset(z1, species == "community"), aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Apparent survival probability", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")


#----------#
#-Figure 1-#
#----------#

#Load lat long
sites <- NULL
sheets <- excel_sheets("~/HMSNO/DataFormat/RawData/trap locations and dates.xlsx")
for(i in 1:length(sheets)){
  Temp <- read_excel(path = "~/HMSNO/DataFormat/RawData/trap locations and dates.xlsx", sheet = sheets[i])
  Temp <- Temp %>% mutate(park = as.factor(sheets[i]),
                          latitude = as.numeric(as.character(latitude)),
                          longitude = as.numeric(as.character(longitude)),
                          start_date = as.Date(start_date),
                          end_date = as.Date(end_date)) 
  sites <- bind_rows(sites, Temp)
}

#Set coordinate system
sites <- sites %>% 
  drop_na("longitude"|"latitude") %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326)

#Domain
africa <- rnaturalearth::ne_countries(continent = "africa", returnclass = "sf")
park <- st_read(dsn = "~/HMSNO/DataFormat/Parks", layer = "d4be616e-6b41-48f9-a554-0744b7474e53202049-1-1aajr88.7n8i")

UDZ <- park %>% filter(name == "Udzungwa Mountains")
VIR <- park %>% filter(name == "Parc national des Volcans")
NFNP <- park %>% filter(name == "Nyungwe")
BIF <- park %>% filter(name == "Bwindi Impenetrable National Park")
NNNP <- park %>% filter(name == "Nouabalé-Ndoki")
KRP <- park %>% filter(name == "Korup")

PARK <- park %>% filter(name == "Udzungwa Mountains" | 
                        name == "Parc national des Volcans" |
                        name == "Nyungwe" |
                        name == "Bwindi Impenetrable National Park" |
                        name == "Nouabalé-Ndoki" |
                        name == "Korup")

Fig1A <- ggplotGrob(ggplot() +
  geom_sf(data = africa) +
  geom_sf(data = sites, aes(color = park), size = 3) +
  #geom_sf(data = PARK, aes(fill = name)) +
  scale_color_brewer(palette = "Accent") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "bottom",
        legend.title = element_blank()
  ))

Fig1B <- ggplotGrob(ggplot() +
  geom_sf(data = UDZ, fill = RColorBrewer::brewer.pal(6, name = "Accent")[1]) +
  geom_sf(data = sites %>% filter(park == "UDZ"), col = "black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
  ))

Fig1C <- ggplotGrob(ggplot() +
  geom_sf(data = VIR, fill =  RColorBrewer::brewer.pal(6, name = "Accent")[2]) +
  geom_sf(data = sites %>% filter(park == "VIR"), col = "black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
  ))

Fig1D <- ggplotGrob(ggplot() +
  geom_sf(data = NFNP, fill = RColorBrewer::brewer.pal(6, name = "Accent")[3]) +
  geom_sf(data = sites %>% filter(park == "NFNP"), col = "black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
  ))

Fig1E <- ggplotGrob(ggplot() +
  geom_sf(data = BIF, fill = RColorBrewer::brewer.pal(6, name = "Accent")[4]) +
  geom_sf(data = sites %>% filter(park == "BIF"), col = "black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
  ))

Fig1F <- ggplotGrob(ggplot() +
  geom_sf(data = NNNP, fill = RColorBrewer::brewer.pal(6, name = "Accent")[5]) +
  geom_sf(data = sites %>% filter(park == "NNNP"), col = "black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
  ))

Fig1G <- ggplotGrob(ggplot() +
  geom_sf(data = KRP, fill = RColorBrewer::brewer.pal(6, name = "Accent")[6]) +
  geom_sf(data = sites %>% filter(park == "KRP"), col = "Black", size = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"
  ))


png(file = "~/HMSNO/PostAnalysis/Figure1.png", res = 600, width = 5, height = 8, units = "in", bg = "transparent")
gridExtra::grid.arrange(gridExtra::arrangeGrob(Fig1A, Fig1B, Fig1C, Fig1D, Fig1E, Fig1F, Fig1G, layout_matrix = matrix(c(1,2,4,6,1,3,5,7), nrow = 2, byrow = TRUE)))
dev.off()

#----------#
#-Figure 2-#
#----------#

Fig2 <- ggplotGrob(data %>% 
  filter(yearID <= nend[park]) %>%
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(aes(col = species)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ park, scales = 'free_x', ncol = 3) +
  scale_x_continuous(breaks = c(2010, 2014, 2018), limits = c(2009, 2019)) +
  theme_few() +
  labs(y = "density (abundance/site)\n", x = expression()))


loc <- Fig2$layout[grep("strip", Fig2$layout$name),1:4]
loc <- loc %>% arrange(t,l)

for(i in 1:max(parkE)){
  Fig2 <- gtable::gtable_add_grob(Fig2, textGrob(LETTERS[i], x = 0.05, y = 0.5, gp = gpar(fontsize = 14, fontface = 2)),
                                  t = loc$t[i], l = loc$l[i])
}

tiff(file = "~/HMSNO/PostAnalysis/Figure2_tws.tiff", res = 600, width = 8, height = 5, units = "in")
grid.arrange(Fig2)
dev.off()

#----------#
#-Figure 3-#
#----------#

Fig3A <- ggplotGrob(ggplot(omega0) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.3) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(y = "Annual apparent\nsurvival probability", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash"))

Fig3B <- ggplotGrob(ggplot(gamma0) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.3) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Annual number of individuals\ngained per site", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash"))

Fig3B$widths <- Fig3A$widths

tiff(file = "~/HMSNO/PostAnalysis/Figure3.tiff", res = 600, width = 6.5, height = 5, units = "in")
grid.arrange(arrangeGrob(Fig3A, Fig3B, ncol = 1, nrow = 2))
grid.text("A", x = 0.15, y = 0.96)
grid.text("B", x = 0.15, y = 0.46)
dev.off()

#----------#
#-Figure 4-#
#----------#

Fig4A <- ggplotGrob(FigO1)
Fig4B <- ggplotGrob(FigZ1)

#----------#
#-Figure 5-#
#----------#

Fig5 <- ggplotGrob(ggplot() +
                   geom_line(data = data %>% filter(yearID <= nend[park]), aes(x = year, y = abundance_site, col = species)) + 
                   geom_ribbon(data = data %>% filter(yearID <= nend[park]), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
                   geom_line(data = data %>% filter(yearID >= nend[park]), aes(x = year, y = abundance_site, col = species), linetype = "dotted") + 
                   geom_ribbon(data = data %>% filter(yearID >= nend[park]), aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.25) +
                   facet_wrap(. ~ park, scales = 'free_x', ncol = 3) +
                   scale_x_continuous(breaks = c(2010, 2020, 2030), limits = c(2009, 2031)) +
                   theme_few() +
                   labs(y = "density (abundance/site)\n", x = expression()))


loc <- Fig5$layout[grep("strip", Fig5$layout$name),1:4]
loc <- loc %>% arrange(t,l)

for(i in 1:max(parkE)){
  Fig5 <- gtable::gtable_add_grob(Fig5, textGrob(LETTERS[i], x = 0.05, y = 0.5, gp = gpar(fontsize = 14, fontface = 2)),
                                  t = loc$t[i], l = loc$l[i])
}

tiff(file = "~/HMSNO/PostAnalysis/Figure5_tws.tiff", res = 600, width = 8, height = 5, units = "in")
grid.arrange(Fig5)
dev.off()
