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


### T1 and T5 ###
# subset by t1 and t5
calf_T1 <- subset_samples(calf_RA, `host_age` == "T1")
calf_T5 <- subset_samples(calf_RA, `host_age` == "T5")

# further subset by sex
calf_T1_female <- subset_samples(calf_T1, host_sex == "female")
calf_T1_male <- subset_samples(calf_T1, host_sex == "male")

calf_T5_female <- subset_samples(calf_T5, host_sex == "female")
calf_T5_male <- subset_samples(calf_T5, host_sex == "male")

# extract core taxa
core_T1_female <- core_members(calf_T1_female, detection = 0.001, prevalence = 0.3)
core_T1_male <- core_members(calf_T1_male, detection = 0.001, prevalence = 0.3)

core_T5_female <- core_members(calf_T5_female, detection = 0.001, prevalence = 0.3)
core_T5_male <- core_members(calf_T5_male, detection = 0.001, prevalence = 0.3)

# make 4-way Venn diagram of core taxa overlap across sex and T1/T5 
list_T1_T5_sex <- list("T1 Female" = core_T1_female, "T1 Male" = core_T1_male,
                 "T5 Female" = core_T5_female, "T5 Male" = core_T5_male)

T1_T5_sex_venn <- ggVennDiagram(x = list_T1_T5_sex, 
                                label_alpha = 0, 
                                edge_size = 0.8,     
                                set_size = 5,         
                                label_size = 4) +
                  scale_fill_gradient(low = "#E3F9FC", high = "#82C1FF") 

# save Venn diagram
ggsave("venn_sex_T1_T5.png", T1_T5_sex_venn, width = 10, height = 10, dpi = 300)


### T5 and T8 ###
# subset by T5 and T8
calf_T5 <- subset_samples(calf_RA, `host_age` == "T5")
calf_T8 <- subset_samples(calf_RA, `host_age` == "T8")

# further subset by sex
calf_T5_female <- subset_samples(calf_T5, host_sex == "female")
calf_T5_male <- subset_samples(calf_T5, host_sex == "male")

calf_T8_female <- subset_samples(calf_T8, host_sex == "female")
calf_T8_male <- subset_samples(calf_T8, host_sex == "male")

# extract core taxa
core_T5_female <- core_members(calf_T5_female, detection = 0.001, prevalence = 0.3)
core_T5_male <- core_members(calf_T5_male, detection = 0.001, prevalence = 0.3)

core_T8_female <- core_members(calf_T8_female, detection = 0.001, prevalence = 0.3)
core_T8_male <- core_members(calf_T8_male, detection = 0.001, prevalence = 0.3)

# make 4-way Venn diagram of core taxa overlap across sex and T5/T8 
list_T5_T8_sex <- list("T5 Female" = core_T5_female, "T5 Male" = core_T5_male,
                       "T8 Female" = core_T8_female, "T8 Male" = core_T8_male)

T5_T8_sex_venn <- ggVennDiagram(x = list_T5_T8_sex, 
                                label_alpha = 0, 
                                edge_size = 0.8,     
                                set_size = 5,         
                                label_size = 4) +
                  scale_fill_gradient(low = "#E5FFDB", high = "#77D184")

# save Venn diagram
ggsave("venn_sex_T5_T8.png", T5_T8_sex_venn, width = 10, height = 10, dpi = 300)


### T7 and T8 ###
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

# make 4-way Venn diagram of core taxa overlap across sex and T7/T8 
list_T7_T8_sex <- list("T7 Female" = core_T7_female, "T7 Male" = core_T7_male,
                "T8 Female" = core_T8_female, "T8 Male" = core_T8_male)

T7_T8_sex_venn <- ggVennDiagram(x = list_T7_T8_sex, 
                                label_alpha = 0, 
                                edge_size = 0.8,     
                                set_size = 5,         
                                label_size = 4) +
                  scale_fill_gradient(low = "#EEEBF5", high = "#B082FF")

# save Venn diagram
ggsave("venn_sex_T7_T8.png", T7_T8_sex_venn, width = 10, height = 10, dpi = 300)
