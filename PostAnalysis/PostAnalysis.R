#-Libraries-#

library(coda)
library(ggplot2)
library(ggthemes)

#-Load data-#

load(file  = "~/HMSNO/DataAnalysis/MSoutput.Rdata")


#-Summary-#

names <- factor(c("Peters's", "Bay", "Harvey's", "White-bellied", "Blue",
          "Ogilby's", "Black-fronted", "Abbott's", "Yellow-backed", "Community"),
          levels =  c("Abbott's", "Bay", "Black-fronted", "Blue",
          "Harvey's", "Ogilby's", "Peters's","White-bellied", 
          "Yellow-backed", "Community"))

out <- summary(output)

alpha0 <- data.frame(species = names, mean = c(plogis(out$statistics[1:9,"Mean"]), out$statistics[13,"Mean"]),
                     l2.5 = c(plogis(out$quantiles[1:9,"2.5%"]), out$quantiles[13,"2.5%"]),
                     l25 = c(plogis(out$quantiles[1:9,"25%"]), out$quantiles[13,"25%"]),
                     u75 = c(plogis(out$quantiles[1:9,"75%"]), out$quantiles[13,"75%"]),
                     u97.5 = c(plogis(out$quantiles[1:9,"97.5%"]), out$quantiles[13,"97.5%"]))

omega0 <- data.frame(species = names, mean = c(plogis(out$statistics[15:23,"Mean"]), out$statistics[14,"Mean"]),
                     l2.5 = c(plogis(out$quantiles[15:23,"2.5%"]), out$quantiles[14,"2.5%"]),
                     l25 = c(plogis(out$quantiles[15:23,"25%"]), out$quantiles[14,"25%"]),
                     u75 = c(plogis(out$quantiles[15:23,"75%"]), out$quantiles[14,"75%"]),
                     u97.5 = c(plogis(out$quantiles[15:23,"97.5%"]), out$quantiles[14,"97.5%"]))

FigA0 <- ggplot(alpha0) + 
  geom_errorbar(aes(x = species, ymin = mean, ymax = mean),
                width = 0.275) +
  geom_errorbar(aes(x = species, ymin = l2.5, ymax = u97.5),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = species, ymin = l25, ymax = u75),
                width = 0, size = 3) +
  coord_cartesian(ylim = c(0, 0.65)) +
  scale_y_continuous(expand = c(0, 0)) +
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
  coord_cartesian(ylim = c(0.3, 1)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  vline
  theme_few() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  labs(y ="Apparent survival probability", x = expression())
