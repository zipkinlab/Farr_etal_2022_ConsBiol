#TO DO: effort by year

#-----------#
#-Libraries-#
#-----------#

library(coda)
library(tidyverse)
library(ggthemes)
library(reshape2)
library(grid)
library(gridExtra)
library(rnaturalearth)

#-----------#
#-Load data-#
#-----------#

files <- list.files(path = "Z:/HMSNO/DataAnalysis", pattern = "chain", full.names = TRUE)

nc <- length(files)

for(i in 1:nc){
  load(files[i])
}

out <- mcmc.list(mget(paste0("chain", 1:nc)))

load(file = "Z:/HMSNO/DataFormat/HMSNO.Adata.Rdata")
load(file = "Z:/HMSNO/DataFormat/HMSNO.Acon.Rdata")
# load(file = "Z:/HMSNO/DataFormat/Effort.Rdata")

#-------------#
#-Attach data-#
#-------------#

attach(HMSNO.data)
attach(HMSNO.con)

#------------#
#-Parameters-#
#------------#

projection <- 10

params <- attr(out[[1]], "dimnames")[[2]][!grepl("beta|Npark", attr(out[[1]], "dimnames")[[2]])]
for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    params <- c(params, paste0("beta0[", i, ", ", s, "]"))
  }
}

for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    for(t in nstart[i]:(nend[i]+projection)){
      params <- c(params, paste0("Npark[", t, ", ", i, ", ", s, "]"))
    }
  }
}

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

#-----------#
#-Extaction-#
#-----------#
       
output <- summary(out[c(1:nc)][,params])

# names <- factor(c("Peters's", "Bay", "Harvey's", "White-bellied", "Blue",
#                   "Black-fronted", "Ogilby's",  "Abbott's", "Yellow-backed", "Community"),
#           levels =  c("Abbott's", "Bay", "Black-fronted", "Blue",
#           "Harvey's", "Ogilby's", "Peters's","White-bellied", 
#           "Yellow-backed", "Community"))

spec.names <- c("callipygus", "dorsalis", "harveyi", "leucogaster", 
                  "monticola", "moschatus", "nigrifrons", "olgilvy",
                  "scriptus", "spadix", "spekii", "sylvilocutor")

park.names <- c("UDZ", "VIR", "NFNP", "BIF", "NNNP", "KRP")

year.names <- function(x){
  return(format(as.Date("2008", format = "%Y") + lubridate::years(x), "%Y"))
}

ni <- nc * length(out[[1]][,1])

N <- array(NA, dim = c(nspecs,max(parkE),(max(nend)+projection),ni))

for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    for(t in nstart[i]:(nend[i]+projection)){
      N[s,i,t,] <- unlist(out[c(1:nc)][,paste("Npark[", t, ", ", i, ", ", s, "]", sep = "")])
    }
  }
}

meanN <- apply(N, MARGIN = c(1,2,3), mean, na.rm = TRUE)
N97.5 <- apply(N, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
N2.5 <- apply(N, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)

data <- reshape2::melt(meanN, varnames = c("species", "parkID", "yearID"), value.name = "abundance")
data$species <- factor(data$species, labels = spec.names)
data$park <- factor(data$park, labels = park.names)
data$year <- as.numeric(year.names(data$yearID))

high <- reshape2::melt(N97.5, varnames = c("species", "parkID", "yearID"), value.name = "97.5%")$`97.5%`
low <- reshape2::melt(N2.5, varnames = c("species", "parkID", "yearID"), value.name = "2.5%")$`2.5%`
data$`97.5%` <- high
data$`2.5%` <- low

# data <- data %>% 
#   left_join(., Effort, by = c("parkID", "yearID"))

data$effort <- nsite[data$parkID]

# for(i in 1:dim(data)[1]){
#   if(is.na(data$nsite[i])){
#     data$nsite[i] <- print(max(data[data$park == data$park[i],"nsite"], na.rm = TRUE))
#   }
# }

#FIX Effort to be by year
data <- data %>%
  mutate(abundance_site = abundance/effort,
         `97.5%_site` = `97.5%`/effort,
         `2.5%_site` = `2.5%`/effort)

# data2 <- data.frame(t(rep(NA, dim(data)[2])))
# 
# for(s in 1:nspec){
#   for(i in parkspec[,s]){
#     data2 <- rbind(data2, data %>% filter(species == names[s] & parkID == i))
#   }
# }

data %>% 
  ggplot(., aes(x = year, y = abundance)) +
  geom_line(aes(col = park)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = park), alpha = 0.5) +
  facet_wrap(. ~ species) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few()

data %>% 
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(aes(col = park)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = park), alpha = 0.5) +
  facet_wrap(. ~ species) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  labs(y = "Abundance")

data %>% 
  ggplot(., aes(x = year, y = abundance)) +
  geom_line(aes(col = species)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ park) +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few()
  # theme(panel.background = element_rect(fill = "transparent", color = NA),
  #       axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
  #       legend.position = "none") +
  # labs(y ="Detection probability", x = expression())

data %>% 
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(aes(col = species)) + 
  geom_ribbon(aes(x = year, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ park, scales = 'free_x') +
  scale_x_continuous(breaks = c(2010, 2015, 2020, 2025, 2030)) +
  theme_few() +
  labs(y = "density (abundance/site)")



#------------------------------#
#-Quasi-extinction probability-#
#------------------------------#

#Summary statistics
measurements <- c("Initial", "Final", "Projected", "Threshold", "QSE", "PopGrowth")

#Extract summary statistics
Data.pop <- array(NA, dim = c(nspecs,max(parkE),length(measurements),ni))
for(s in 1:nspecs){
  for(i in parkS[s]:parkE[s]){
    Data.pop[s,i,1,] <- N[s,i,nstart[i],]
    Data.pop[s,i,2,] <- N[s,i,nend[i],]
    Data.pop[s,i,3,] <- N[s,i,nend[i]+projection,]
    Data.pop[s,i,4,] <- N[s,i,nstart[i],] * 0.2
    Data.pop[s,i,5,] <- N[s,i,nend[i]+projection,] < (N[s,i,nstart[i],] * 0.2)
    tmp2 <- 1
    for(t in (nstart[i]+1):nend[i]){
      tmp1 <- N[s,i,t,]/N[s,i,t-1,]
      tmp2 <- tmp1 * tmp2
    }
    Data.pop[s,i,6,] <- tmp2^(1/(nend[i]-nstart[i]))
  }
}

#Calculate quantiles from MCMC
Data.mean <- apply(Data.pop, MARGIN = c(1,2,3), mean, na.rm = TRUE)
Data.97.5 <- apply(Data.pop, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
Data.2.5 <- apply(Data.pop, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)

#Format summary statistics into dataframe
qse <- full_join(full_join(dcast(melt(Data.mean, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements), 
                  dcast(melt(Data.2.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
                 by = c("species", "parkID")),
                 dcast(melt(Data.97.5, varnames = c("species", "parkID", "measurements")), formula = species + parkID ~ measurements),
                 by = c("species", "parkID"))
colnames(qse) <- c("species", "park", paste0(measurements, rep(c("_mean", "_lower", "_upper"), each = length(measurements))))
qse$species <- factor(qse$species, labels = spec.names)
qse$park <- factor(qse$park, labels = park.names)
qse <- qse[complete.cases(qse),]
qse <- qse[,c("species", "park", paste0(rep(measurements, each = 3), c("_mean", "_lower", "_upper")))]

#Remove columns
qse <- qse[,-c(16,17)]
qse$Percent_change <- (qse[,"Final_mean"] - qse[,"Initial_mean"])/qse[,"Initial_mean"]
# qse$Projected_change <- (qse[,"Projected_mean"] - qse[,"Final_mean"])/qse[,"Final_mean"]

#--------------------#
#-Summary statistics-#
#--------------------#

#Inital abundances by species
print(qse %>% group_by(species) %>% select(Initial_mean, Initial_lower, Initial_upper),
      n = 26)

print(qse %>% group_by(species) %>% summarize(sd(Initial_mean)),
      n = 26)

#Inital abundances by park
print(qse %>% group_by(park) %>% arrange(park) %>% select(Initial_mean, Initial_lower, Initial_upper),
      n = 26)

print(qse %>% group_by(park) %>% arrange(park) %>% summarize(sd(Initial_mean)),
      n = 26)

#Final abundances by species
print(qse %>% group_by(species) %>% select(Final_mean, Final_lower, Final_upper),
      n = 26)

#Final abundances by park
print(qse %>% group_by(park) %>% arrange(park) %>% select(Final_mean, Final_lower, Final_upper),
      n = 26)

#Population growth by species
print(qse %>% group_by(species) %>% select(PopGrowth_mean, PopGrowth_lower, PopGrowth_upper),
      n = 26)
#Population growth by park
print(qse %>% group_by(park) %>% arrange(park) %>% select(PopGrowth_mean, PopGrowth_lower, PopGrowth_upper),
      n = 26)

#Percent change by species
print(qse %>% group_by(species) %>% select(Percent_change),
      n = 26)
#Percent change by park
print(qse %>% group_by(park) %>% arrange(park) %>% select(Percent_change),
      n = 26)


#Quasi-extinction by species
qse %>% group_by(species) %>% summarize(median(QSE_mean), max(QSE_mean))
#Quasi-extinction by park
qse %>% group_by(park) %>% summarize(median(QSE_mean), max(QSE_mean))


#-----------------------------#
#-Species-specific parameters-#
#-----------------------------#

omega0 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(plogis(output$statistics[grep("omega0", params),"Mean"]), output$statistics["mu.o0","Mean"]),
                     l2.5 = c(plogis(output$quantiles[grep("omega0", params),"2.5%"]), output$quantiles["mu.o0","2.5%"]),
                     l25 = c(plogis(output$quantiles[grep("omega0", params),"25%"]), output$quantiles["mu.o0","25%"]),
                     u75 = c(plogis(output$quantiles[grep("omega0", params),"75%"]), output$quantiles["mu.o0","75%"]),
                     u97.5 = c(plogis(output$quantiles[grep("omega0", params),"97.5%"]), output$quantiles["mu.o0","97.5%"]))

omega1 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(output$statistics[grep("omega1", params),"Mean"], output$statistics["mu.o1","Mean"]),
                     l2.5 = c(output$quantiles[grep("omega1", params),"2.5%"], output$quantiles["mu.o1","2.5%"]),
                     l25 = c(output$quantiles[grep("omega1", params),"25%"], output$quantiles["mu.o1","25%"]),
                     u75 = c(output$quantiles[grep("omega1", params),"75%"], output$quantiles["mu.o1","75%"]),
                     u97.5 = c(output$quantiles[grep("omega1", params),"97.5%"], output$quantiles["mu.o1","97.5%"]))

gamma0 <- data.frame(species = factor(c(spec.names, "community"), levels = c(spec.names, "community")), mean = c(exp(output$statistics[grep("gamma0", params),"Mean"]), exp(output$statistics["mu.g0","Mean"])),
                     l2.5 = c(exp(output$quantiles[grep("gamma0", params),"2.5%"]), exp(output$quantiles["mu.g0","2.5%"])),
                     l25 = c(exp(output$quantiles[grep("gamma0", params),"25%"]), exp(output$quantiles["mu.g0","25%"])),
                     u75 = c(exp(output$quantiles[grep("gamma0", params),"75%"]), exp(output$quantiles["mu.g0","75%"])),
                     u97.5 = c(exp(output$quantiles[grep("gamma0", params),"97.5%"]), exp(output$quantiles["mu.g0","97.5%"])))



FigO1 <- ggplot(omega1) + 
  geom_hline(yintercept = 0, col = "red") +
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.3) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(-5, 5)) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Edge effect (log-scale)", x = expression()) +
  geom_vline(xintercept=(13 - 0.5), linetype = "dotdash")


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
#-Figure 2-#
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

sites <- as.factor()

#Set coordinate system
sites <- sites %>% 
  drop_na("longitude"|"latitude") %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326)

#Domain
africa <- ne_countries(continent = "africa", returnclass = "sf")

ggplot() +
  geom_sf(data = africa) +
  geom_sf(data = sites, aes(color = park), size = 3) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())


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
