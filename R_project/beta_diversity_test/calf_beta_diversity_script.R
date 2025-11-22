library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(vegan)

#### Load in RData ####
load("calf_phyloseq_rare.RData")
load("calf_phyloseq_final.RData")


#subset out samples with diet & milk
calf_phyloseq_no_diet <- subset_samples(calf_phyloseq_rare, host_age != "not applicable")
nsamples(calf_phyloseq_no_diet) #check sample number
samp_dat_wdiv_no_diet <- data.frame(sample_data(calf_phyloseq_no_diet), estimate_richness(calf_phyloseq_no_diet))

set.seed(42)
#### Beta diversity #####
calf_dm_braycurtis <- vegdist(t(otu_table(calf_phyloseq_no_diet)), method="bray") # Bray-curtis
calf_pcoa_bc <- ordinate(calf_phyloseq_no_diet, method="PCoA", distance=calf_dm_braycurtis)

adonis2(calf_dm_braycurtis ~ host_age, data=samp_dat_wdiv_no_diet)


#Time Points & Sex
calf_gg_pcoa <- plot_ordination(calf_phyloseq_no_diet, calf_pcoa_bc, color = "host_age", shape="host_sex") +
  labs(pch="Sex", col = "Timepoints")
calf_gg_pcoa

ggsave("calf_gg_pcoa.png"
       , calf_gg_pcoa
       , height=4, width=5)

#Only Timepoints
calf_gg_pcoa_timepoints <- plot_ordination(calf_phyloseq_no_diet, calf_pcoa_bc, color = "host_age") +
  stat_ellipse(aes(color = host_age), type = "t", level = 0.95) +
  labs(col = "Timepoints")
calf_gg_pcoa_timepoints

ggsave("calf_gg_pcoa_timepoints.png"
       , calf_gg_pcoa_timepoints
       , height=4, width=5)

#Only Sex
calf_gg_pcoa_sex <- plot_ordination(calf_phyloseq_no_diet, calf_pcoa_bc, color = "host_sex") +
  stat_ellipse(aes(color = host_sex), type = "t", level = 0.95) +
  labs(col = "Sex")
calf_gg_pcoa_sex

ggsave("calf_gg_pcoa_sex.png"
       , calf_gg_pcoa_sex
       , height=4, width=5)

