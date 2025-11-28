library(DESeq2)
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(vegan)
library(ggplot2)

set.seed(42)
#### Load data ####
load("all_phyloseq_objects/calf_phyloseq_final.RData")
rank_names(calf_phyloseq_final)
calf_phyloseq_final_genus <- tax_glom(calf_phyloseq_final, taxrank = "Genus")
#save(calf_phyloseq_final_genus, file="all_phyloseq_objects/calf_phyloseq_final_genus.RData")
load("all_phyloseq_objects/calf_phyloseq_final_genus.RData")

#### DESeq sex ####
calf_plus1 <- transform_sample_counts(calf_phyloseq_final, function(x) x+1)
calf_sex_deseq <- phyloseq_to_deseq2(calf_plus1, ~`host_sex`)
DESEQ_calf_sex <- DESeq(calf_sex_deseq)
res_calf_sex <- results(DESEQ_calf_sex, name = "host_sex_male_vs_female", tidy = TRUE)
View(res_calf_sex)