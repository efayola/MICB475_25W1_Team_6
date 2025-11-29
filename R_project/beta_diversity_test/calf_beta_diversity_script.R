library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(vegan)

#### Load in RData ####
load("all_phyloseq_objects/calf_phyloseq_rare.RData")
load("all_phyloseq_objects/calf_phyloseq_no_diet.RData")
load("all_phyloseq_objects/calf_phyloseq_male.RData")
load("all_phyloseq_objects/calf_phyloseq_female.RData")

#calf_phyloseq_no_diet is rarefied (419 + 55 + 8 - 2)
calf_phyloseq_no_diet_and_T9 <- subset_samples(calf_phyloseq_no_diet, host_age != "T9")
nsamples(calf_phyloseq_no_diet_and_T9)


set.seed(42)
#### Beta diversity #####
calf_dm_braycurtis <- vegdist(t(otu_table(calf_phyloseq_no_diet_and_T9)), method="bray") # Bray-curtis
calf_pcoa_bc <- ordinate(calf_phyloseq_no_diet_and_T9, method="PCoA", distance=calf_dm_braycurtis)

#PERMANOVA test
samp_dat_wdiv <- data.frame(sample_data(calf_phyloseq_no_diet_and_T9), estimate_richness(calf_phyloseq_no_diet_and_T9))
adonis2(calf_dm_braycurtis ~ host_age, data=samp_dat_wdiv)
#           Df SumOfSqs     R2      F Pr(>F)    
# Model      7   17.421 0.1033 6.7636  0.001 ***
# Residual 411  151.232 0.8967                  
# Total    418  168.654 1.0000    

adonis2(calf_dm_braycurtis ~ host_sex, data=samp_dat_wdiv)
#           Df SumOfSqs      R2      F Pr(>F)  
# Model      1    0.654 0.00388 1.6241  0.019 *
# Residual 417  167.999 0.99612                
# Total    418  168.654 1.00000                


#Time Points & Sex
calf_gg_pcoa <- plot_ordination(calf_phyloseq_no_diet_and_T9, calf_pcoa_bc, color = "host_age", shape="host_sex") +
  labs(pch="Sex", col = "Timepoints")
calf_gg_pcoa

ggsave("beta_diversity_test/calf_gg_pcoa.png"
       , calf_gg_pcoa
       , height=4, width=5)

#Only Timepoints
calf_gg_pcoa_timepoints <- plot_ordination(calf_phyloseq_no_diet_and_T9, calf_pcoa_bc, color = "host_age") +
  stat_ellipse(aes(color = host_age), type = "t", level = 0.95) +
  labs(col = "Timepoints")
calf_gg_pcoa_timepoints

ggsave("beta_diversity_test/calf_gg_pcoa_timepoints.png"
       , calf_gg_pcoa_timepoints
       , height=4, width=5)

#Only Sex
calf_gg_pcoa_sex <- plot_ordination(calf_phyloseq_no_diet_and_T9, calf_pcoa_bc, color = "host_sex") +
  stat_ellipse(aes(color = host_sex), type = "t", level = 0.95) +
  labs(col = "Sex")
calf_gg_pcoa_sex

ggsave("beta_diversity_test/calf_gg_pcoa_sex.png"
       , calf_gg_pcoa_sex
       , height=4, width=5)

# Only Male Beta Diversity
set.seed(42)
calf_dm_braycurtis_male <- vegdist(t(otu_table(calf_phyloseq_male)), method="bray") # Bray-curtis
calf_pcoa_bc_male <- ordinate(calf_phyloseq_male, method="PCoA", distance=calf_dm_braycurtis_male)

calf_gg_pcoa_only_male_timepoints <- plot_ordination(calf_phyloseq_male, calf_pcoa_bc_male, color = "host_age") +
  stat_ellipse(aes(color = host_age), type = "t", level = 0.95) +
  labs(col = "Timepoints")
calf_gg_pcoa_only_male_timepoints

ggsave("beta_diversity_test/calf_gg_pcoa_only_male_timepoints.png"
       , calf_gg_pcoa_only_male_timepoints
       , height=4, width=5)


