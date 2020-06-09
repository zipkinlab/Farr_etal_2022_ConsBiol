#-Libraries-#

library(coda)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(reshape2)

#-Load data-#

load(file = "~/HMSNO/MSdynout.Rdata")
load(file = "~/HMSNO/DataFormat/HMSNO.data.Rdata")
load(file = "~/HMSNO/DataFormat/HMSNO.con.Rdata")
load(file = "~/HMSNO/DataFormat/Effort.Rdata")


# Effort <- Data %>% 
#   mutate(siteID = as.numeric(as.factor(`deployment ID`)),
#          year = as.integer(as.factor(Year))) %>% 
#   group_by(park, year) %>% 
#   summarize(nsite = n_distinct(siteID), days = sum(days))
# Effort$park <- as.factor(Effort$park)

# Effort <- data.frame(park = HMSNO.data$park, t(HMSNO.data$days)) %>% 
#   melt(id = "park") %>%
#   mutate(year = as.numeric(variable)) %>%
#   group_by(park, year) %>%
#   summarize(nsite = n(), days = sum(value, na.rm = TRUE))
#   
#   group_by(park)

nspec <- HMSNO.con$nspec
nparks <- max(HMSNO.data$park)
nyrs <- max(HMSNO.con$yr)
nyrstr <- HMSNO.con$nyrstr
nyrend <- HMSNO.con$nyrend
parkspec <- HMSNO.con$parkspec
#-Summary-#

names <- factor(c("Peters's", "Bay", "Harvey's", "White-bellied", "Blue",
                  "Black-fronted", "Ogilby's",  "Abbott's", "Yellow-backed", "Community"),
          levels =  c("Abbott's", "Bay", "Black-fronted", "Blue",
          "Harvey's", "Ogilby's", "Peters's","White-bellied", 
          "Yellow-backed", "Community"))

output <- MSdyn.o
nchains <- length(output)
out <- summary(output[c(1:nchains)][,543:578])
ni <- length(unlist(output[c(1)][,1])) * nchains

N <- array(NA, dim = c(nspec,nparks,nyrs,ni))

for(s in 1:nspec){
  for(i in 1:nparks){
    for(t in nyrstr[i]:nyrend[i]){
      N[s,i,t,] <- unlist(output[c(1:nchains)][,paste("Npark[", t, ", ", i, ", ", s, "]", sep = "")])
    }
  }
}


meanN <- apply(N, MARGIN = c(1,2,3), mean, na.rm = TRUE)
N97.5 <- apply(N, MARGIN = c(1,2,3), quantile, probs = 0.975, na.rm = TRUE)
N2.5 <- apply(N, MARGIN = c(1,2,3), quantile, probs = 0.025, na.rm = TRUE)

data <- reshape2::melt(meanN, varnames = c("species", "parkID", "yearID"), value.name = "abundance")
data$species <- factor(data$species,
                       labels = c("Peters's", "Bay", "Harvey's", "White-bellied", "Blue",
                                      "Black-fronted", "Ogilby's",  "Abbott's", "Yellow-backed"))
data$park <- as.factor(data$park)

high <- reshape2::melt(N97.5, varnames = c("species", "parkID", "yearID"), value.name = "97.5%")$`97.5%`
low <- reshape2::melt(N2.5, varnames = c("species", "parkID", "yearID"), value.name = "2.5%")$`2.5%`
data$`97.5%` <- high
data$`2.5%` <- low

data <- data %>% 
  left_join(., Effort, by = c("parkID", "yearID")) %>%
  mutate(abundance_site = abundance/nsite,
         `97.5%_site` = `97.5%`/nsite,
         `2.5%_site` = `2.5%`/nsite,
         abundance_days = abundance/days,
         `97.5%_days` = `97.5%`/days,
         `2.5%_days` = `2.5%`/days,
         parkID = as.factor(parkID))

# data2 <- data.frame(t(rep(NA, dim(data)[2])))
# 
# for(s in 1:nspec){
#   for(i in parkspec[,s]){
#     data2 <- rbind(data2, data %>% filter(species == names[s] & parkID == i))
#   }
# }

data %>% 
  #filter(park != "4" & park != "6") %>%
  ggplot(., aes(x = yearID, y = abundance)) +
  geom_line(aes(col = parkID)) + 
  geom_ribbon(aes(x = yearID, ymin = `2.5%`, ymax = `97.5%`, fill = parkID), alpha = 0.5) +
  facet_wrap(. ~ species) +
  theme_few()

data %>% 
  #filter(park != "4" & park != "6") %>%
  ggplot(., aes(x = yearID, y = abundance_site)) +
  geom_line(aes(col = parkID)) + 
  geom_ribbon(aes(x = yearID, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = parkID), alpha = 0.5) +
  facet_wrap(. ~ species) +
  theme_few() +
  labs(y = "Abundance")

data %>% 
  #filter(park != "4" & park != "6") %>%
  ggplot(., aes(x = yearID, y = abundance)) +
  geom_line(aes(col = species)) + 
  geom_ribbon(aes(x = yearID, ymin = `2.5%`, ymax = `97.5%`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ parkID) +
  theme_few()
  # theme(panel.background = element_rect(fill = "transparent", color = NA),
  #       axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
  #       legend.position = "none") +
  # labs(y ="Detection probability", x = expression())

Fig1 <- data %>% 
  mutate(year = factor(as.character(yearID))) %>%
  #filter(park != "4" & park != "6") %>%
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(aes(col = species)) + 
  geom_line(aes(x = yearID, col = species)) +
  geom_ribbon(aes(x = yearID, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ parkID, scales='free_x', labeller = labeller(parkID = c("1" = "Korup",
                                                        "2" = "Nouabale-Ndoki", 
                                                        "3" = "Nyungwe", 
                                                        "4" = "Bwindi", 
                                                        "5" = "Virunga", 
                                                        "6" = "Udzungwa"))) +
  labs(y = "Abundance per site", x = "Year") +
  theme_few()

Fig2 <- data %>% 
  mutate(year = factor(as.character(yearID + 2008))) %>%
  #filter(park != "4" & park != "6") %>%
  ggplot(., aes(x = year, y = abundance_site)) +
  geom_line(aes(col = species)) + 
  geom_line(aes(x = yearID, col = species)) +
  geom_ribbon(aes(x = yearID, ymin = `2.5%_site`, ymax = `97.5%_site`, fill = species), alpha = 0.5) +
  facet_wrap(. ~ parkID, scales='free', labeller = labeller(parkID = c("1" = "Korup",
                                                                       "2" = "Nouabale-Ndoki", 
                                                                       "3" = "Nyungwe", 
                                                                       "4" = "Bwindi", 
                                                                       "5" = "Virunga", 
                                                                       "6" = "Udzungwa"))) +
  labs(y = "Abundance per site", x = "Year") +
  theme_few()

# alpha0 <- data.frame(species = names, mean = c(plogis(out$statistics[1:9,"Mean"]), out$statistics[13,"Mean"]),
#                      l2.5 = c(plogis(out$quantiles[1:9,"2.5%"]), out$quantiles[13,"2.5%"]),
#                      l25 = c(plogis(out$quantiles[1:9,"25%"]), out$quantiles[13,"25%"]),
#                      u75 = c(plogis(out$quantiles[1:9,"75%"]), out$quantiles[13,"75%"]),
#                      u97.5 = c(plogis(out$quantiles[1:9,"97.5%"]), out$quantiles[13,"97.5%"]))

omega0 <- data.frame(species = names, mean = c(plogis(out$statistics[6:14,"Mean"]), out$statistics[3,"Mean"]),
                     l2.5 = c(plogis(out$quantiles[6:14,"2.5%"]), out$quantiles[3,"2.5%"]),
                     l25 = c(plogis(out$quantiles[6:14,"25%"]), out$quantiles[3,"25%"]),
                     u75 = c(plogis(out$quantiles[6:14,"75%"]), out$quantiles[3,"75%"]),
                     u97.5 = c(plogis(out$quantiles[6:14,"97.5%"]), out$quantiles[3,"97.5%"]))

omega1 <- data.frame(species = names, mean = c(out$statistics[15:23,"Mean"], out$statistics[4,"Mean"]),
                     l2.5 = c(out$quantiles[15:23,"2.5%"], out$quantiles[4,"2.5%"]),
                     l25 = c(out$quantiles[15:23,"25%"], out$quantiles[4,"25%"]),
                     u75 = c(out$quantiles[15:23,"75%"], out$quantiles[4,"75%"]),
                     u97.5 = c(out$quantiles[15:23,"97.5%"], out$quantiles[4,"97.5%"]))

omega2 <- data.frame(species = names, mean = c(out$statistics[24:32,"Mean"], out$statistics[5,"Mean"]),
                     l2.5 = c(out$quantiles[24:32,"2.5%"], out$quantiles[5,"2.5%"]),
                     l25 = c(out$quantiles[24:32,"25%"], out$quantiles[5,"25%"]),
                     u75 = c(out$quantiles[24:32,"75%"], out$quantiles[5,"75%"]),
                     u97.5 = c(out$quantiles[24:32,"97.5%"], out$quantiles[5,"97.5%"]))

# FigA0 <- ggplot(alpha0) + 
#   geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                 width = 0.275) +
#   geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
#                 width = 0, size = 1.25) +
#   geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
#                 width = 0, size = 3) +
#   coord_cartesian(ylim = c(0, 0.65)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   geom_hline(yintercept = 0, alpha = 0.75) +
#   theme_few() +
#   theme(panel.background = element_rect(fill = "transparent", color = NA),
#         axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#         legend.position = "none") +
#     labs(y ="Detection probability", x = expression())

FigO0 <- ggplot(omega0) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.275) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Apparent survival probability", x = expression())


FigO1 <- ggplot(omega1) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.275) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(-7, 7)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Effect of density", x = expression())

FigO2 <- ggplot(omega2) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.275) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Effect of distance", x = expression())

#-Single-species-#

# species <- c("callipygus", "dorsalis", "harveyi", "leucogaster", "monticola", "nigrifrons")
# names <- c("Peters's", "Bay", "Harvey's", "White-bellied", "Blue", "Black-fronted")
# 
# alpha0 <- beta0 <- beta1 <- beta2 <- beta3 <- 
#   beta4 <- gamma0 <- omega0 <- omega1 <- omega2 <- omega3 <- rep(NA, 6)
# 
# for(i in 1:length(species)){
#   load(file = paste("Z:/HMSNO/", species[i], ".Rdata", sep=""))
#   output <- summary(as.mcmc.list(out))
#   alpha0 <- rbind(alpha0, c("mean" = output$statistics[1,1], output$quantiles[1,]))
#   beta0 <- rbind(beta0, c("mean" = output$statistics[3,1], output$quantiles[3,]))
#   beta1 <- rbind(beta1, c("mean" = output$statistics[4,1], output$quantiles[4,]))
#   beta2 <- rbind(beta2, c("mean" = output$statistics[5,1], output$quantiles[5,]))
#   beta3 <- rbind(beta3, c("mean" = output$statistics[6,1], output$quantiles[6,]))
#   beta4 <- rbind(beta4, c("mean" = output$statistics[7,1], output$quantiles[7,]))
#   gamma0 <- rbind(gamma0, c("mean" = output$statistics[8,1], output$quantiles[8,]))
#   omega0 <- rbind(omega0, c("mean" = output$statistics[9,1], output$quantiles[9,]))
#   omega1 <- rbind(omega1, c("mean" = output$statistics[10,1], output$quantiles[10,]))
#   omega2 <- rbind(omega2, c("mean" = output$statistics[11,1], output$quantiles[11,]))
#   omega3 <- rbind(omega3, c("mean" = output$statistics[12,1], output$quantiles[12,]))
# }
# 
# alpha0 <- alpha0[-1,]; beta0 <- beta0[-1,]; beta1 <- beta1[-1,]; beta2 <- beta2[-1,]; beta3 <- beta3[-1,]; beta4 <- beta4[-1,];
# gamma0 <- gamma0[-1,]; omega0 <- omega0[-1,]; omega1 <- omega1[-1,]; omega2 <- omega2[-1,]; omega3 <- omega3[-1,]
# 
# out <- list(alpha0, beta0, beta1, beta2, beta3, beta4, gamma0, omega0, omega1, omega2, omega3)
# 
# for(i in 1:11){
#   data.frame(species = names, out[[i]]) %>%
#     ggplot(.) + 
#     geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
#                   width = 0.275) +
#     geom_errorbar(aes(x = species, ymin = X2.5., ymax = X97.5.),
#                   width = 0, size = 1.25) +
#     geom_errorbar(aes(x = species, ymin = X25., ymax = X75.),
#                   width = 0, size = 3) +
#     geom_hline(yintercept = 0, alpha = 0.75) +
#     theme_few() +
#     theme(panel.background = element_rect(fill = "transparent", color = NA),
#           axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
#           legend.position = "none") +
#     labs(y = expression(), x = expression())
#   
# }
