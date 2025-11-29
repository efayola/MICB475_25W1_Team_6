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

load("all_phyloseq_objects/calf_phyloseq_rare.RData")

# Plot bar plot of taxonomy
plot_bar(calf_phyloseq_rare, fill="Phylum") 

# Convert to relative abundance
calf_RA <- transform_sample_counts(calf_phyloseq_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
calf_phylum <- tax_glom(calf_RA, taxrank = "Phylum", NArm=FALSE)

gg_taxa <- plot_bar(calf_phylum, fill="Phylum") + 
  facet_wrap(.~subject, scales = "free_x")
gg_taxa

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)