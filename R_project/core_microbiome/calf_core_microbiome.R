#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("R_project/calf_phyloseq_final.RData")

#### "core" microbiome ####

# Convert to relative abundance
calf_RA <- transform_sample_counts(calf_phyloseq_final, fun=function(x) x/sum(x))

# subset by pre-weaning and weaning stages
calf_pre <- subset_samples(calf_RA, `host_age` %in% c("T1","T2","T3","T4","T5"))
calf_weaning <- subset_samples(calf_RA, `host_age` %in% c("T6","T7","T8"))

# further subset by sex
calf_pre_female <- subset_samples(calf_pre, host_sex == "female")
calf_pre_male <- subset_samples(calf_pre, host_sex == "male")

calf_weaning_female <- subset_samples(calf_weaning, host_sex == "female")
calf_weaning_male <- subset_samples(calf_weaning, host_sex == "male")

# extract core taxa
core_pre_female <- core_members(calf_pre_female, detection = 0.001, prevalence = 0.3)
core_pre_male <- core_members(calf_pre_male, detection = 0.001, prevalence = 0.3)

core_weaning_female <- core_members(calf_weaning_female, detection = 0.001, prevalence = 0.3)
core_weaning_male <- core_members(calf_weaning_male,   detection = 0.001, prevalence = 0.3)


# make Venn diagram
# Venn diagram of female vs male pre-weaning calves 
pre_venn <- ggVennDiagram(x = list(Female = core_pre_female, Male = core_pre_male)) +
            scale_fill_gradient(low = "grey95", high ="steelblue") +
            ggtitle("Core Taxa: Pre-weaning (Females vs Males)") +
            theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave("venn_pre.png", pre_venn)

# Venn diagram of female vs male weaning calves                            
weaning_venn <- ggVennDiagram(x=list(Female = core_weaning_female, Male = core_weaning_male)) +
                scale_fill_gradient(low = "lavender", high ="purple") +
                ggtitle("Core Taxa: Weaning (Females vs Males)") +
                theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave("venn_weaning.png", weaning_venn)

# 4-way Venn diagram of core taxa overlap across sex and developmental stage
list_all <- list("Pre-Weaning Female" = core_pre_female, "Pre-Weaning Male" = core_pre_male,
                "Weaning Female" = core_weaning_female, "Weaning Male" = core_weaning_male)

all_venn <- ggVennDiagram(x = list_all) +
            scale_fill_gradient(low = "white", high = "skyblue") +
            labs(title = "Core Taxa Overlap Across Sex and Developmental Stage") +
            theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave("venn_all.png", all_venn)

