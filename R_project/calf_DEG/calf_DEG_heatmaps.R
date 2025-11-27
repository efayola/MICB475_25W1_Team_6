library(tidyverse)
library(phyloseq)
library(DESeq2)

#load data
load("all_phyloseq_objects/calf_phyloseq_rare.RData")

# #### DESeq ####
calf_plus1 <- transform_sample_counts(calf_phyloseq_rare, function(x) x+1)
calf_timepoints_deseq <- phyloseq_to_deseq2(calf_plus1, ~`host_age`)
DESEQ_calf_timepoints <- DESeq(calf_timepoints_deseq)
res_calf_timepoints <- results(DESEQ_calf_timepoints,
                               name = "host_age_vs_age",
                               tidy = TRUE)
View(res_calf_sex)