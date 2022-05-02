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
library(abind)

#-----------#
#-Load data-#
#-----------#

load(file = "~/HMSNO/DataFormat/HMSNO.data.Rdata")
load(file = "~/HMSNO/DataFormat/HMSNO.con.Rdata")

attach(HMSNO.data)
attach(HMSNO.con)

#-------------------#
#-Load model output-#
#-------------------#

pattern <- "chain"

#Retrospective analysis
files <- list.files(path = "~/HMSNO/DataAnalysis/CaseStudy", pattern = pattern, full.names = TRUE)

nc <- 5

for(i in 1:nc){
  load(files[i])
}

out <- mcmc.list(mget(ls()[grep(pattern, ls())]))

rm(list = ls()[grep(pattern, ls())])

#------------#
#-Parameters-#
#------------#

params <- attr(out[[1]], "dimnames")[[2]][!grepl("Npark|park.surv|pop.surv|park.gain|pop.gain", attr(out[[1]], "dimnames")[[2]])]

#-------------#
#-Convergence-#
#-------------#

Rhat <- gelman.diag(out[c(1:nc)][,params])
if(all(Rhat[[1]][,1] < 1.1)){
  print("Converged")
}else{
  tmp <- as.numeric(which(Rhat[[1]][,1] > 1.1))
  print("Not converged")
  print(params[tmp])
  traceplot(out[c(1:nc)][,params[tmp]])
}

#------------#
#-Extraction-#
#------------#

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
       
spec.names <- c("C. callipygus", "C. dorsalis", "C. harveyi", "C. leucogaster", 
                  "P. monticola", "N. moschatus", "C. nigrifrons", "C. ogilbyi",
                  "T. scriptus", "C. spadix", "T. spekii", "C. silvicultor")

park.names <- c("UDZ", "VNP", "NFNP", "BIF", "NNNP", "KRP")

year.names <- function(x){
  return(format(as.Date("2008", format = "%Y") + lubridate::years(x), "%Y"))
}

ni <- nc * length(out[[1]][,1])

#---------#
#-Results-#
#---------#

N <- array(NA, dim = c(nspecs,max(parkE),11,ni))
park.surv <- array(NA, dim = c(max(parkE),ni))
pop.surv <- array(NA, dim = c(nspecs,max(parkE),ni))

for(i in 1:6){
  park.surv[i,] <- unlist(out[c(1:nc)][,paste("park.surv[", i, "]", sep = "")])
}

for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    pop.surv[s,i,] <- unlist(out[c(1:nc)][,paste("pop.surv[", i, ", ", s, "]", sep = "")])
    for(t in nstart[i]:nend[i]){
      N[s,i,t,] <- unlist(out[c(1:nc)][,paste("Npark[", t, ", ", i, ", ", s, "]", sep = "")])
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

#-------------------#
#-Covariate effects-#
#-------------------#

output <- summary(out[c(1:nc)][,params])

omega0 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(plogis(output$statistics[grep("omega0", params),"Mean"]), output$statistics["mu.o0","Mean"]),
                     l2.5 = c(plogis(output$quantiles[grep("omega0", params),"2.5%"]), output$quantiles["mu.o0","2.5%"]),
                     l25 = c(plogis(output$quantiles[grep("omega0", params),"25%"]), output$quantiles["mu.o0","25%"]),
                     u75 = c(plogis(output$quantiles[grep("omega0", params),"75%"]), output$quantiles["mu.o0","75%"]),
                     u97.5 = c(plogis(output$quantiles[grep("omega0", params),"97.5%"]), output$quantiles["mu.o0","97.5%"]))

popdata <- full_join(popdata, omega0 %>% mutate(park = "TEAM"))
popdata <- full_join(popdata, parkdata %>% mutate(species = "community"))
popdata$park <- factor(popdata$park, levels = c("KRP", "NNNP", "UDZ", "BIF", "NFNP", "VNP", "TEAM"))
popdata$species <- factor(popdata$species, levels = c(spec.names, "community"))

#-Gains-#

gamma0 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(exp(output$statistics[grep("gamma0", params),"Mean"]), exp(output$statistics["mu.g0","Mean"])),
                     l2.5 = c(exp(output$quantiles[grep("gamma0", params),"2.5%"]), exp(output$quantiles["mu.g0","2.5%"])),
                     l25 = c(exp(output$quantiles[grep("gamma0", params),"25%"]), exp(output$quantiles["mu.g0","25%"])),
                     u75 = c(exp(output$quantiles[grep("gamma0", params),"75%"]), exp(output$quantiles["mu.g0","75%"])),
                     u97.5 = c(exp(output$quantiles[grep("gamma0", params),"97.5%"]), exp(output$quantiles["mu.g0","97.5%"])))

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

#TEAM community mean
popdata[popdata$species == "community" & popdata$park == "TEAM",]

#Range of species survival
sd.o0 <- mean(1/sqrt(unlist(out[c(1:nc)][,"tau.o0"])))
sd.o0; quantile(1/sqrt(unlist(out[c(1:nc)][,"tau.o0"])), c(0.025, 0.975))

plogis(logit(output$statistics["mu.o0","Mean"]) + sd.o0) -
plogis(logit(output$statistics["mu.o0","Mean"]) - sd.o0)

#Range of park survival
sd.eps.o <- mean(1/sqrt(unlist(out[c(1:nc)][,"tau.eps.o"])))
sd.eps.o; quantile(1/sqrt(unlist(out[c(1:nc)][,"tau.eps.o"])), c(0.025, 0.975))

plogis(logit(output$statistics["mu.o0","Mean"]) + sd.eps.o) -
  plogis(logit(output$statistics["mu.o0","Mean"]) - sd.eps.o)


#Species vs park magnitude
data.frame(species = factor(spec.names), 
           mean = output$statistics[grep("omega0", params),"Mean"] - logit(output$statistics["mu.o0","Mean"]),
           l2.5 = output$quantiles[grep("omega0", params),"2.5%"] - logit(output$quantiles["mu.o0","2.5%"]),
           l25 = output$quantiles[grep("omega0", params),"25%"] - logit(output$quantiles["mu.o0","25%"]),
           u75 = output$quantiles[grep("omega0", params),"75%"] - logit(output$quantiles["mu.o0","75%"]),
           u97.5 = output$quantiles[grep("omega0", params),"97.5%"] - logit(output$quantiles["mu.o0","97.5%"]))

data.frame(park = factor(park.names),
           mean = output$statistics[grep("eps.o", params),"Mean"][-7],
           l2.5 = output$quantiles[grep("eps.o", params),"2.5%"][-7],
           l25 = output$quantiles[grep("eps.o", params),"25%"][-7],
           u75 = output$quantiles[grep("eps.o", params),"75%"][-7],
           u97.5 = output$quantiles[grep("eps.o", params),"97.5%"][-7])

#Maximum species survival
max(popdata[popdata$species != "community",'mean'])

popdata[popdata$species != "community",]

popdata[popdata$species == "community" & popdata$park != "TEAM",]

#TEAM community gains
gamma0[gamma0$species == 'community',]

#----------#
#-Figure 2-#
#----------#

pos <- position_dodge(0.75, preserve = "total")

data$park <- factor(data$park, labels = c("KRP" = "Korup National Park",
                                          "NNNP" = "Nouabale-Ndoki National Park",
                                          "UDZ" = "Udzungwa National Park",
                                          "BIF" = "Bwindi National Park",
                                          "NFNP" = "Nyungwe Forest National Park",
                                          "VNP" = "Volcanoes National Park"))

Figure2 <-  ggplot(data %>% drop_na(.), aes(x = year, y = abundance_site)) +
  geom_point(col = "grey") +
  geom_line(aes(col = species)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ park, scales = 'free', ncol = 2, dir = "v") +
  scale_x_continuous(breaks = c(2010, 2012, 2014, 2016, 2018, 2020)) +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "italic")) +
  labs(y = "Relative abundance/site\n", x = "Year")


tiff(file = "~/HMSNO/PostAnalysis/Figure2.tiff", res = 600, width = 6, height = 6, units = "in")
Figure2
dev.off()

#----------#
#-Figure 3-#
#----------#

Fig3A <- popdata %>% filter(park == "TEAM" | species == "community") %>%
  mutate(xname = ifelse(park == "TEAM", as.character(species), as.character(park))) %>%
  mutate(xname = factor(xname, levels = c(spec.names, "community", park.names))) %>%
  ggplot(., aes(x = xname, color = park)) +
  geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
               stat = "identity", size = 0.75, width = 0.5) +
  scale_color_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "black")) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  labs(y = "Apparent survival", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")

Fig3B <- gamma0 %>%
  ggplot(., aes(x = species)) +
  geom_boxplot(aes(ymin = l2.5, lower = l25, middle = mean, upper = u75, ymax = u97.5),
               stat = "identity", size = 0.75, width = 0.5, col = "black") +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5, face = "italic"),
        legend.title = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  labs(y = "Individuals gained \nper site", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")

Fig3A <- ggplotGrob(Fig3A)
Fig3B <- ggplotGrob(Fig3B)

Fig3B$heights <- Fig3A$heights

Legend <- Fig3A$grobs[[15]]

Fig3A$grobs[[15]] <- nullGrob()

Fig3B$widths <- Fig3A$widths

Fig3A$widths[9] <- unit(0, "cm")

Fig3B$widths[9] <- sum(unit(5.4, "cm"), unit(2.5, "point"))

tiff(file = "~/HMSNO/PostAnalysis/Figure3.tiff", res = 600, width = 8, height = 5, units = "in")
grid.arrange(arrangeGrob(Fig3A, Fig3B), bottom = "\n \n", left = "")
grid.draw(textGrob("A", x = unit(0.5, "cm"), y = unit(12.25, "cm"),
                   gp=grid::gpar(fontsize=14, fontface = 2)))
grid.draw(textGrob("B", x = unit(0.5, "cm"), y = unit(6.75, "cm"),
                   gp=grid::gpar(fontsize=14, fontface = 2)))
pushViewport(viewport(x = unit(17, "cm"), y = unit(4.75, "cm")))
grid.draw(Legend)
dev.off()

#------------------#
#-Simulation study-#
#------------------#

#List of all filenames of 75 site output
filenames_75 <- list.files(path = "~/HMSNO/DataAnalysis/SimulationStudy/SimulationOutput/output_75", full.names = TRUE)
filenames_25 <- list.files(path = "~/HMSNO/DataAnalysis/SimulationStudy/SimulationOutput/output_25", full.names = TRUE)

filenames <- c(filenames_75, filenames_25)

#File length of each scenario
n75 <- length(filenames_75)
n25 <- length(filenames_25)

#Load first file
load(filenames[1])

#Initialize vector for all output
out <- output[,2:7]

#Unconverged index
ii <- NULL

#Harvout parameters from files
for(i in 2:(n75 + n25)){
  load(filenames[i])
  if(max(output[,c(5,8)], na.rm = TRUE) < 1.1){
    out <- abind(out, output[,2:7], along = 3)
  }else{
    ii <- c(ii, i)
  }
}

jj <- seq(1,c(n75+n25))[-ii]
nn75 <- max(which(jj<500))
nn25 <- length(jj)

#-Summary dataframe-#

ms75_75 <- apply((out[1:12,2,1:nn75] - out[1:12,1,1:nn75])/out[1:12,1,1:nn75], MARGIN = 1, FUN = quantile, probs = 0.75, na.rm = TRUE)
ms50_75 <- apply((out[1:12,2,1:nn75] - out[1:12,1,1:nn75])/out[1:12,1,1:nn75], MARGIN = 1, FUN = quantile, probs = 0.5, na.rm = TRUE)
ms25_75 <- apply((out[1:12,2,1:nn75] - out[1:12,1,1:nn75])/out[1:12,1,1:nn75], MARGIN = 1, FUN = quantile, probs = 0.25, na.rm = TRUE)
msmax_75 <-  ((ms75_75 - ms25_75) * 1.5) + ms75_75
msmin_75 <- ms25_75 - ((ms75_75 - ms25_75) * 1.5)

ss75_75 <- apply((out[1:12,5,1:nn75] - out[1:12,1,1:nn75])/out[1:12,1,1:nn75], MARGIN = 1, FUN = quantile, probs = 0.75, na.rm = TRUE)
ss50_75 <- apply((out[1:12,5,1:nn75] - out[1:12,1,1:nn75])/out[1:12,1,1:nn75], MARGIN = 1, FUN = quantile, probs = 0.5, na.rm = TRUE)
ss25_75 <- apply((out[1:12,5,1:nn75] - out[1:12,1,1:nn75])/out[1:12,1,1:nn75], MARGIN = 1, FUN = quantile, probs = 0.25, na.rm = TRUE)
ssmax_75 <-  ((ss75_75 - ss25_75) * 1.5) + ss75_75
ssmin_75 <- ss25_75 - ((ss75_75 - ss25_75) * 1.5)

ms75_25 <- apply((out[1:12,2,(nn75+1):nn25] - out[1:12,1,(nn75+1):nn25])/out[1:12,1,(nn75+1):nn25], MARGIN = 1, FUN = quantile, probs = 0.75, na.rm = TRUE)
ms50_25 <- apply((out[1:12,2,(nn75+1):nn25] - out[1:12,1,(nn75+1):nn25])/out[1:12,1,(nn75+1):nn25], MARGIN = 1, FUN = quantile, probs = 0.5, na.rm = TRUE)
ms25_25 <- apply((out[1:12,2,(nn75+1):nn25] - out[1:12,1,(nn75+1):nn25])/out[1:12,1,(nn75+1):nn25], MARGIN = 1, FUN = quantile, probs = 0.25, na.rm = TRUE)
msmax_25 <-  ((ms75_25 - ms25_25) * 1.5) + ms75_25
msmin_25 <- ms25_25 - ((ms75_25 - ms25_25) * 1.5)

ss75_25 <- apply((out[1:12,5,(nn75+1):nn25] - out[1:12,1,(nn75+1):nn25])/out[1:12,1,(nn75+1):nn25], MARGIN = 1, FUN = quantile, probs = 0.75, na.rm = TRUE)
ss50_25 <- apply((out[1:12,5,(nn75+1):nn25] - out[1:12,1,(nn75+1):nn25])/out[1:12,1,(nn75+1):nn25], MARGIN = 1, FUN = quantile, probs = 0.5, na.rm = TRUE)
ss25_25 <- apply((out[1:12,5,(nn75+1):nn25] - out[1:12,1,(nn75+1):nn25])/out[1:12,1,(nn75+1):nn25], MARGIN = 1, FUN = quantile, probs = 0.25, na.rm = TRUE)
ssmax_25 <-  ((ss75_25 - ss25_25) * 1.5) + ss75_25
ssmin_25 <- ss25_25 - ((ss75_25 - ss25_25) * 1.5)

data <- data.frame(rbind(cbind(msmin_75, ms25_75, ms50_75, ms75_75, msmax_75), 
                         cbind(ssmin_75, ss25_75, ss50_75, ss75_75, ssmax_75),
                         cbind(msmin_25, ms25_25, ms50_25, ms75_25, msmax_25), 
                         cbind(ssmin_25, ss25_25, ss50_25, ss75_25, ssmax_25)))

params <- as.factor(output[,1])
data$params <- rep(params[1:12], 4)
data$framework <- rep(rep(c("multi-species", "single-species"), each = 12), 2)
data$sites <- rep(c("75", "25"), each = 24)
colnames(data)[1:5] <- c("qmin","q25","q50","q75","qmax")
data <- data %>% separate(params, c("params", "species.type"))
data$framework <- factor(data$framework, labels = c("multi-species" = "Multi-species",
                                                    "single-species" = "Single-species"))
data$species.type <- factor(data$species.type, labels = c("common" = "Common", 
                      "elusive" = "Elusive", 
                      "rare" = "Rare"))
data$params <- factor(data$params,
                      labels = c("gamma" = "Gains", 
                                 "lambda" = "Initial Abundance", 
                                "omega" = "Apparent Survival", 
                                "r" = "Detection Probability"))
data$params <- factor(data$params, levels = c("Detection Probability", "Initial Abundance", "Apparent Survival", "Gains"))
data$sites <- factor(data$sites, labels = c("25" = "25 sites",
                                            "75" = "75 sites"))
data$species.type <- factor(data$species.type, levels = c("Common", "Rare", "Elusive"))

#----------#
#-Figure 1-#
#----------#

Figure1 <- ggplotGrob(ggplot(data, aes(x = sites, fill = framework)) +
  geom_boxplot(aes(ymin = qmin, lower = q25, middle = q50, upper = q75, ymax = qmax), 
               stat = "identity", size = 0.75) +
  scale_fill_manual(values = c("white", "grey")) +
  facet_grid(params ~ species.type, scale = "free", 
             labeller = label_wrap_gen(width = 2, multi_line = TRUE)) +
  geom_hline(yintercept = 0, col = "black", size = 1) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  labs(y = "Relative bias", x = expression()))

loc <- Figure1$layout[grep("panel", Figure1$layout$name),1:4]
loc <- loc %>% arrange(t,l)

for(i in 1:dim(loc)[1]){
  Figure1 <- gtable::gtable_add_grob(Figure1, grid::textGrob(LETTERS[i], x = unit(0.1, "in"), y = unit(1.5, "in"),
                                                             gp = grid::gpar(fontsize = 14,
                                                                             fontface = 2)),
                                     t = loc$t[i], l = loc$l[i])
}


tiff(file = "~/HMSNO/PostAnalysis/Figure1.tiff", res = 600, width = 6, height = 8, units = "in")
grid::grid.draw(Figure1)
dev.off()