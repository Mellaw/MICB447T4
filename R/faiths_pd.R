#Calculating and plotting Faith's phylogentic distance
#Date last modified: Dec 19 2019

#libraries
library(tidyverse)
library(ggplot2)
library(ggsignif)

# Load data
feed_faith_pd <- read_tsv("feedfaithspd.tsv")
ad_faith_pd <- read_tsv("adfaithspd.tsv")

# Create Faith's PD box plot for feeding method
ggplot(feed_faith_pd, aes( x = feed, y = faith_pd)) +
  geom_boxplot() +
  xlab("Feeding Method") + 
  ylab("Faith's PD") +
  geom_signif(comparisons = list(c("breast", "formula")), 
              map_signif_level=TRUE) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

# Create Faith's PD box plot for AD status
ggplot(ad_faith_pd, aes( x = eczema_inf, y = faith_pd)) +
  geom_boxplot() +
  xlab("AD Status") + 
  ylab("Faith's PD") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
