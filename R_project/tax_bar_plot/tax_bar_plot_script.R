library(DESeq2)
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(vegan)
library(ggplot2)
library(patchwork)
library(ggrepel)

#### Load Data ####
# load("all_phyloseq_objects/calf_phyloseq_rare.RData")
# rank_names(calf_phyloseq_rare)
# 
# #merges species that have the same taxonomy at genus level
# calf_phyloseq_rare_phylum <- tax_glom(calf_phyloseq_rare, taxrank = "Phylum")
# save(calf_phyloseq_rare_phylum, file="all_phyloseq_objects/calf_phyloseq_rare_phylum.RData")

load("all_phyloseq_objects/calf_phyloseq_rare_phylum.RData")
nsamples(calf_phyloseq_rare_phylum) #480 removed 3 by rarefractino

calf_phyloseq_no_diet_and_T9 <- subset_samples(calf_phyloseq_rare_phylum, host_age != "not applicable" & host_age !="T9")
nsamples(calf_phyloseq_no_diet_and_T9)

# Convert to relative abundance
calf_RA <- transform_sample_counts(calf_phyloseq_no_diet_and_T9, function(x) x/sum(x))

# Code to glom later by families
# # To remove black bars, "glom" by phylum first
# calf_phylum <- tax_glom(calf_RA, taxrank = "Phylum", NArm=FALSE)

#Aggregate relative abundances by Phylum
calf_RA_top10 <- prune_taxa(taxa_sums(calf_RA), calf_RA)

gg_taxa <- plot_bar(top10_phyla, fill="Phylum") + 
  facet_wrap(.~host_age, scales = "free_x")
gg_taxa

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)