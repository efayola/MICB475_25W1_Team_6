#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("all_phyloseq_objects/calf_phyloseq_final.RData")

#### "core" microbiome ####

# Convert to relative abundance
calf_RA <- transform_sample_counts(calf_phyloseq_final, fun=function(x) x/sum(x))

# subset by T7 and T8
calf_T7 <- subset_samples(calf_RA, `host_age` == "T7")
calf_T8 <- subset_samples(calf_RA, `host_age` == "T8")

# further subset by sex
calf_T7_female <- subset_samples(calf_T7, host_sex == "female")
calf_T7_male <- subset_samples(calf_T7, host_sex == "male")

calf_T8_female <- subset_samples(calf_T8, host_sex == "female")
calf_T8_male <- subset_samples(calf_T8, host_sex == "male")

# extract core taxa
core_T7_female <- core_members(calf_T7_female, detection = 0.001, prevalence = 0.3)
core_T7_male <- core_members(calf_T7_male, detection = 0.001, prevalence = 0.3)

core_T8_female <- core_members(calf_T8_female, detection = 0.001, prevalence = 0.3)
core_T8_male <- core_members(calf_T8_male, detection = 0.001, prevalence = 0.3)

# make Venn diagram
# 4-way Venn diagram of core taxa overlap across sex and T7/T8 
list_all <- list("T7 Female" = core_T7_female, "T7 Male" = core_T7_male,
                "T8 Female" = core_T8_female, "T8 Male" = core_T8_male)

all_venn <- ggVennDiagram(x = list_all) +
            scale_fill_gradient(low = "#E7E1F7", high = "#AD92F0") +
            labs(title = "Core Taxa Overlap Across Sex and Timepoints 7 & 8") +
            theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# save Venn diagram
ggsave("venn_sex_T7_T8.png", all_venn, width = 10, height = 10, dpi = 300)
