library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(vegan)
library(ggplot2)

load("all_phyloseq_objects/calf_phyloseq_male.RData")

set.seed(42)
#Male Beta Diversity Heat Map
calf_dm_braycurtis_male <- vegdist(t(otu_table(calf_phyloseq_male)), method="bray") # Bray-curtis
#all possibility for male
samp_dat_wdiv_male <- data.frame(sample_data(calf_phyloseq_male), estimate_richness(calf_phyloseq_male))
adonis2(calf_dm_braycurtis_male ~ host_age, data=samp_dat_wdiv_male)

set.seed(42)
# idea
# Subset the phyloseq object different timepoints
calf_phyloseq_tp1_2 <- subset_samples(calf_phyloseq_male, host_age %in% c("T1", "T2"))
#Calculate Bray-Curtis distance for the subsetted data
calf_dm_braycurtis_tp1_2 <- vegdist(t(otu_table(calf_phyloseq_tp1_2)), method="bray")
# Create sample data frame for the subsetted object
samp_dat_wdiv_tp1_2 <- data.frame(sample_data(calf_phyloseq_tp1_2), estimate_richness(calf_phyloseq_tp1_2))
# Run PERMANOVA (adonis2) on subsetted timepoints
adonis2(calf_dm_braycurtis_tp1_2 ~ host_age, data=samp_dat_wdiv_tp1_2)
#save PERMANOVA results
permanova_result <- adonis2(calf_dm_braycurtis_tp1_2 ~ host_age, data=samp_dat_wdiv_tp1_2)
# only extract the first row of PERMANOVA results
model_stats <- permanova_result[1, ]
# Convert first row of PERMANOVA results into a data frame
permanova_df <- data.frame(
  Variable = "host_age",
  R2 = model_stats$R2,
  F_value = model_stats$F,
  P_value = model_stats$`Pr(>F)`
)


# data_matrix <- matrix(rnorm(100, mean = 0, sd = 1), nrow = 10, ncol = 10)
# ggplot(data_matrix) + 
#   geom_tile()

