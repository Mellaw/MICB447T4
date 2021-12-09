# libraries ---------------------------------------------------------------
library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(ggsignif)

library(phyloseq)
library(dplyr)
library(DESeq2)
library(ape)
library(vegan)
library(svglite)

# helper functions --------------------------------------------------------

# Calculate relative abundance
calculate_relative_abundance <- function(x) x / sum(x)

# Calculate geometric mean 
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Data Prep ---------------------------------------------------------------
metadata <- import_qiime_sample_data("QIIME imports/inf_full/infant_metadata.txt") 

tree <- read_tree("QIIME imports/inf_full/tree.nwk")
tree <- multi2di(tree)

#Samples were already subset for life_stage = infant in qiime
biom <- import_biom("QIIME imports/inf_full/table-with-taxonomy.biom")

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


#Subset samples to metadata categories of interest
#remove samples with no data for feed and eczema_inf
ad <-  c("no", "yes")

breast_sub <- subset_samples(min_10k,life_stage == "Infant" & 
                              feed == "breast" & 
                              eczema_inf %in% ad)

formula_sub <- subset_samples(min_10k,life_stage == "Infant" & 
                               feed == "formula" & 
                               eczema_inf %in% ad)


# Common variables for differential & relative abundance ------------------
ab_tax <- function (physeq_ob, rank) {
  total_counts <- taxa_sums(physeq_ob)
  relative_abundance <- calculate_relative_abundance(total_counts)
  abundant <- relative_abundance > 0.0005 
  
  abundant_taxa_diff <- prune_taxa(abundant, physeq_ob) 
  abundant_genera <- tax_glom(abundant_taxa_diff, taxrank = rank, NArm = TRUE)
  ab_genera <- prune_samples(sample_sums(abundant_genera) > 0, abundant_genera)
}

breast_phy <- ab_tax(breast_sub, "Genus")
formula_phy <- ab_tax(formula_sub, "Genus")

# Differential Abundance AD --------------------------------------------------
sample_data(breast_phy)$eczema_inf  <-
  factor(sample_data(breast_phy)$eczema_inf ,
         levels = c("no", "yes"))

sample_data(formula_phy)$eczema_inf  <-
  factor(sample_data(formula_phy)$eczema_inf ,
         levels = c("no", "yes"))

deseq_ob <- function (physeq_ob) {
  deseq_inf <- phyloseq_to_deseq2(physeq_ob, ~ eczema_inf)
  geo_means <- apply(counts(deseq_inf), 1, calculate_gm_mean)
  deseq_inf <- estimateSizeFactors(deseq_inf, geoMeans = geo_means) 
  deseq_inf <- DESeq(deseq_inf, fitType = "local")
  inf_diff_abund <- results(deseq_inf)
}

sig_change <- function(diff_abund, alpha) {
  significant_inf <- as.data.frame(diff_abund)
  significant_inf <- filter(significant_inf, padj < alpha)
}

sig_tax <- function (ab_genera, significant_inf) {
  genera_df <- as.data.frame(tax_table(ab_genera))
  significant_inf <- merge(significant_inf, genera_df, by = "row.names")
  significant_inf <- arrange(significant_inf, log2FoldChange)
}


breast_do <- deseq_ob(breast_phy)
sig_breast <- sig_change(breast_do, 0.05)
sig_taxa_br <- sig_tax(breast_phy,sig_breast)

formula_do <- deseq_ob(formula_phy)
sig_formula <- sig_change(formula_do, 0.05)
sig_taxa_fm <- sig_tax(formula_phy,sig_formula)

dim(sig_taxa_br)
dim(sig_taxa_fm)



# Differential Abundance plots ------------------------------------------------
sig_genera <- function (significant_inf) {
  mutate(significant_inf,
         Genus = factor(Genus, levels = Genus))
}

format_genera_lab <- function (significant_inf) {
  significant_inf <- significant_inf %>% 
    mutate(Genus = str_remove(Genus, "g__")) %>% 
    mutate(Genus = str_remove(Genus, "_1"))
}

diff_genera_br <- sig_genera(sig_taxa_br)
diff_genera_br_ft <- format_genera_lab(diff_genera_br)

diff_genera_fm <- sig_genera(sig_taxa_fm)
diff_genera_fm_ft <- format_genera_lab(diff_genera_fm)

diff_genera_br_ft <- diff_genera_br_ft %>% 
  mutate(ad_non = ifelse(log2FoldChange > 0, "AD", "Non AD"))

diff_genera_fm_ft <- diff_genera_fm_ft %>% 
  mutate(ad_non = ifelse(log2FoldChange > 0, "AD", "Non AD"))

  
  
ggplot(diff_genera_br_ft, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential Abunance",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_classic()+
  theme(legend.title = element_blank())

ggplot(diff_genera_fm_ft, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential Abunance",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_classic()+
  theme(legend.title = element_blank())


# Relative Abundance breastfed  ------------------------------------------------

#step 4
br_RA <- transform_sample_counts(breast_sub, calculate_relative_abundance)

#step 5
br_counts <- taxa_sums(breast_sub)
relative_abundance_br <- calculate_relative_abundance(br_counts)
abundant_br <- relative_abundance_br > 0.0005
abundant_br_RA_taxa <- prune_taxa(abundant_br, br_RA)


#step 6
br_abundant_genera_rel <- tax_glom(abundant_br_RA_taxa, taxrank = "Genus")

#Parabacteroides
para <- subset_taxa(br_abundant_genera_rel, Genus == "g__Parabacteroides")
otu_table(para)
para_long <- psmelt(para)

#Erysipelatoclostridium
ery <- subset_taxa(br_abundant_genera_rel, Genus == "g__Erysipelatoclostridium")
otu_table(ery)
ery_long <- psmelt(ery)

#Enterococcus
enter <- subset_taxa(br_abundant_genera_rel, Genus == "g__Enterococcus") 
otu_table(enter)
enter_long <- psmelt(enter)

#Clostridium 
clost <- subset_taxa(br_abundant_genera_rel, Genus == "g__Clostridium_sensu_stricto_1")
otu_table(clost)
clost_long <- psmelt(clost)


# Relative abundance formula-fed -------------------------------------------

#step 4
fm_RA <- transform_sample_counts(formula_sub, calculate_relative_abundance)

#step 5
fm_counts <- taxa_sums(formula_sub)
relative_abundance_fm <- calculate_relative_abundance(fm_counts)
abundant_fm <- relative_abundance_fm > 0.0005
abundant_fm_RA_taxa <- prune_taxa(abundant_fm, fm_RA)

#step 6
fm_abundant_genera_rel <- tax_glom(abundant_fm_RA_taxa, taxrank = "Genus")

#Faecalitalea
fae <- subset_taxa(fm_abundant_genera_rel, Genus == "g__Faecalitalea")
otu_table(fae)
fae_long <- psmelt(fae)

#Candidatus_Stoquefichus
can <- subset_taxa(fm_abundant_genera_rel, Genus == "g__Candidatus_Stoquefichus")
otu_table(can)
can_long <- psmelt(can)

#Acinetobacter
aci <- subset_taxa(fm_abundant_genera_rel, Genus == "g__Acinetobacter") 
otu_table(aci)
aci_long <- psmelt(aci)


# Plotting Relative abundance Breastfed ---------------------------------------------
para_gg <- ggplot(para_long, aes(x = eczema_inf, y = Abundance*100, fill = eczema_inf)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme(legend.position = "none") +
  labs(title = "Parabacteroides",
       x     = "AD Status",
       y     = "Relative Abundance")+
  theme_classic() +
  theme(legend.position = "none")

ery_gg <- ggplot(ery_long, aes(x = eczema_inf, y = Abundance*100, fill = eczema_inf)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme(legend.position = "none") +
  labs(title = "Erysipelatoclostridium",
       x     = "AD Status",
       y     = "Relative Abundance [%]")+
  theme_classic() +
  theme(legend.position = "none")

enter_gg <- ggplot(enter_long, aes(x = eczema_inf, y = Abundance*100, fill = eczema_inf)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 2.75)) +
  theme(legend.position = "none") +
  labs(title = "Enterococcus",
       x     = "AD Status",
       y     = "Relative Abundance [%]")+
  theme_classic() +
  theme(legend.position = "none")

clost_gg <- ggplot(clost_long, aes(x = eczema_inf, y = Abundance*100, fill = eczema_inf)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme(legend.position = "none") +
  labs(title = "Clostridium sensu stricto",
       x     = "AD Status",
       y     = "Relative Abundance [%]")+
  theme_classic() +
  theme(legend.position = "none")

ggsave(file="./images/breastfed/clost.png", plot=clost_gg,
       width=3.5, height=2.5, units ="in")
ggsave(file="./images/breastfed/enter.png", plot=enter_gg,
       width=3.5, height=2.5, units ="in")
ggsave(file="./images/breastfed/ery.png", plot=ery_gg,
       width=3.5, height=2.5, units ="in")
ggsave(file="./images/breastfed/para.png", plot=para_gg,
       width=3.5, height=2.5, units ="in")

# Plotting relative abundance formula-fed ---------------------------------

fae_gg <- ggplot(fae_long, aes(x = eczema_inf, y = Abundance*100, fill = eczema_inf)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 0.025)) +
  theme(legend.position = "none") +
  labs(title = "Faecalitalea",
       x     = "AD Status",
       y     = "Relative Abundance [%]")+
  theme_classic()+
  theme(legend.position = "none")

can_gg <- ggplot(can_long, aes(x = eczema_inf, y = Abundance*100, fill = eczema_inf)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 0.12)) +
  theme(legend.position = "none") +
  labs(title = "Candidatus Stoquefichus",
       x     = "AD Status",
       y     = "Relative Abundance [%]")+
  theme_classic() +
  theme(legend.position = "none")

aci_gg <- ggplot(aci_long, aes(x = eczema_inf, y = Abundance*100, fill = eczema_inf)) +
  geom_boxplot(color = "black", fill = "grey") +
  #scale_fill_manual(values = c("#93cde6", "#b0d7f7"),
                    #labels = c("No AD", "AD"))+
  coord_cartesian(ylim = c(0, 0.12)) +
  theme(legend.position = "none") +
  labs(title = "Acinetobacter",
       x     = "AD Status",
       y     = "Relative Abundance[%]") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(file="./images/formula/fae.png", plot=fae_gg,
       width=3.5, height=2.5, units ="in")
ggsave(file="./images/formula/can.png", plot=can_gg,
       width=3.5, height=2.5, units ="in")
ggsave(file="./images/formula/aci.png", plot=aci_gg,
       width=3.5, height=2.5, units ="in")


# Abundance for C. difficile ----------------------------------------------

tax_br_cdiff <- tax_glom(abundant_br_RA_taxa, taxrank = "Species")
tax_fm_cdiff <- tax_glom(abundant_fm_RA_taxa, taxrank = "Species")


c_diff_br <- subset_taxa(tax_br_cdiff, Species == "s__Clostridioides_difficile")
otu_table(c_diff_br )
c_diff_br_long <- psmelt(c_diff_br)

c_diff_fm <- subset_taxa(tax_fm_cdiff, Species == "s__Clostridioides_difficile")
otu_table(c_diff_fm )
c_diff_fm_long <- psmelt(c_diff_fm)

#Plot relative abundance

c_diff_br <- ggplot(c_diff_br_long, aes(x = eczema_inf, y = Abundance)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 0.1)) +
  theme(legend.position = "none") +
  labs(title = "Clostridioides difficile",
       x     = "AD Status",
       y     = "Relative Abundance")+
  theme_classic()

c_diff_fm <- ggplot(c_diff_fm_long, aes(x = eczema_inf, y = Abundance)) +
  geom_boxplot(color = "black", fill = "grey") +
  coord_cartesian(ylim = c(0, 0.1)) +
  theme(legend.position = "none") +
  labs(title = "Clostridioides difficile",
       x     = "AD Status",
       y     = "Relative Abundance")+
  theme_classic()

ggsave(file="./images/breastfed/c_diff_br.png", plot=c_diff_br,
       width=3.5, height=3.5, units ="in")
ggsave(file="./images/formula/c_diff_fm.png", plot=c_diff_fm,
       width=3.5, height=3.5, units ="in")





