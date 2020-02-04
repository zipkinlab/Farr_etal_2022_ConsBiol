#-Libraries-#

library(coda)
library(ggplot2)
library(ggthemes)

#-Load data-#

load(file  = "~/HMSNO/DataAnalysis/MSoutput.Rdata")


#-Summary-#

names <- as.factor(c("Peters's", "Bay", "Harvey's", "White-bellied", "Blue",
          "Ogilby's", "Black-fronted", "Abbott's", "Yellow-backed", "Community"))

out <- summary(output)

alpha0 <- data.frame(species = names, mean = plogis(out$statistics[1:9,"Mean"]),
                     l2.5 = plogis(out$quantiles[1:9,"2.5%"]),
                     l25 = plogis(out$quantiles[1:9,"25%"]),
                     u75 = plogis(out$quantiles[1:9,"75%"]),
                     u97.5 = plogis(out$quantiles[1:9,"97.5%"]))

omega0 <- data.frame(species = names, mean = plogis(out$statistics[15:23,"Mean"]),
                     l2.5 = plogis(out$quantiles[15:23,"2.5%"]),
                     l25 = plogis(out$quantiles[15:23,"25%"]),
                     u75 = plogis(out$quantiles[15:23,"75%"]),
                     u97.5 = plogis(out$quantiles[15:23,"97.5%"]))

FigA0 <- ggplot(alpha0) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.275) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  scale_y_continuous(expand = c(0, 0)) +
  #coord_cartesian(ylim = c(0, 1)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
    labs(y ="Detection probability", x = expression())

FigO0 <- ggplot(omega0) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.275) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(0.25, 1)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Apparent survival probability", x = expression())
