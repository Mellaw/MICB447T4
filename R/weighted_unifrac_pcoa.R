#Plotting weighted UniFrac PCoA
#Date last modified: Dec 19 2021

# libraries ---------------------------------------------------------------
library(phyloseq)
library(DESeq2)
library(tidyverse)
library(data.table)
library(cluster)
library(dplyr)
library(ape)
library(vegan)

library(factoextra)
library(dendextend)
library(svglite)

#not all library are not needed

# Data Prep ---------------------------------------------------------------
set.seed(800)

metadata <- import_qiime_sample_data("filtered_infant_metadata.tsv") 
tree <- read_tree("tree.nwk")
tree <- multi2di(tree)

#Samples were alreadt subset for life_stage = infant in qiime
biom <- import_biom("table-with-taxonomy.biom")


#Create Physeq object
physeq = merge_phyloseq(biom, metadata, tree)

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(physeq) 


#Check sequencing depth
depth <- sample_sums(physeq)
unique(depth)
count(sample_sums(physeq) >= 10000)

#Subset to only include samples with a min read depth of 10000
min_10k <- prune_samples(sample_sums(physeq) >= 10000, physeq)

feed_type <- c("breast", "formula")
ad <-  c("no", "yes")

#Subset samples to metadata categories of interest
#remove samples with no data for feed and eczema_inf
inf_sub <- subset_samples(min_10k,
                          life_stage == "Infant" & 
                            feed %in% feed_type & 
                            eczema_inf %in% ad) 


# Rarefy ------------------------------------------------------------------
#Rarefy at 26750 reads #  TO CHANGE
physeq_rar <- rarefy_even_depth(physeq, sample.size = 26750, rngseed = TRUE)

# ß-diversity PCoA plots for feed types------------------------------------------------------------------
physeq_rar <- rarefy_even_depth(physeq, sample.size = 26750)
ord <- ordinate(physeq_rar, method = "PCoA", distance = "wunifrac")
plot_ordination(physeq_rar,
                ord,
                color = "feed") +
  # Define title of plot
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 15))

# ß-diversity PCoA plots for AD------------------------------------------------------------------
physeq_rar <- rarefy_even_depth(physeq, sample.size = 26750)
ord <- ordinate(physeq_rar, method = "PCoA", distance = "wunifrac")
plot_ordination(physeq_rar,
                ord,
                color = "eczema_inf") +
  # Define title of plot
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 15))

# ß-diversity PCoA plots for time points/life stage------------------------------------------------------------------ 
physeq_rar <- rarefy_even_depth(physeq, sample.size = 26750)
ord <- ordinate(physeq_rar, method = "PCoA", distance = "wunifrac")
plot_ordination(physeq_rar,
                ord,
                color = "age_category") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 15)) +
  scale_colour_manual(values = c("lightsteelblue1", "lightskyblue3", "lightskyblue4", "dodgerblue4", "black"),
                      labels = c("0.5 months", "2 months", "4 months", "6 months", "9 months"))



# alpha-diversity boxplots for feedtype------------------------------------------------------------------ 
feed_comp <- list(c("breast", "formula"))

faith <- estimate_pd(physeq_rar) %>%
  setDT(keep.rownames="#SampleID") %>%
  inner_join(meta)

faith_plot <- ggplot(faith, aes(x = as.factor(feed), y = PD)) +
  geom_boxplot() +
  xlab("Feeding type") +
  ylab("Faith's Phylogenetic Diversity") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  stat_compare_means(comparisons = feed_comp, method = "kruskal.test",
                     label = "p.format", label.y = 19, label.x = 1.3, size = 6)
