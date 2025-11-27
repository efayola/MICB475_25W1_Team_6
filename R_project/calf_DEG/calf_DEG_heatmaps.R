library(DESeq2)
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(vegan)
library(ggplot2)

#load data
load("all_phyloseq_objects/calf_phyloseq_rare.RData")

# #### DESeq ####
# calf_plus1 <- transform_sample_counts(calf_phyloseq_rare, function(x) x+1)
# calf_timepoints_deseq <- phyloseq_to_deseq2(calf_plus1, ~`host_age`)
# calf_timepoints_deseq$host_age <- relevel(calf_timepoints_deseq$host_age, ref = "T1")
# DESEQ_calf_timepoints <- DESeq(calf_timepoints_deseq)

#save(DESEQ_calf_timepoints, file="DESEQ_calf_timepoints.RData")
load("calf_DEG/DESEQ_calf_timepoints.RData")

resultsNames(DESEQ_calf_timepoints)
# [1] "Intercept"                     "host_age_not.applicable_vs_T1" "host_age_T2_vs_T1"            
# [4] "host_age_T3_vs_T1"             "host_age_T4_vs_T1"             "host_age_T5_vs_T1"            
# [7] "host_age_T6_vs_T1"             "host_age_T7_vs_T1"             "host_age_T8_vs_T1"            
# [10] "host_age_T9_vs_T1"   


res_calf_timepoints_T2_T1 <- results(DESEQ_calf_timepoints, 
                               name = "host_age_T2_vs_T1",
                               tidy = TRUE)
View(res_calf_timepoints_T2_T1)
# To get table of results
sigASVs_T2_T1 <- res_calf_timepoints_T2_T1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_T2_T1)
#To get ASV
sigASVs_T2_T1_vec <- sigASVs_T2_T1 %>%
  pull(ASV)

T1_T2_pair <- subset_samples(calf_phyloseq_rare, host_age %in% c("T1", "T2")) %>%
  transform_sample_counts(function(x) x + 1)

# Prune phyloseq file
T1_T2_DESeq <- prune_taxa(sigASVs_T2_T1_vec,T1_T2_pair)
sigASVs_T2_T1 <- tax_table(T1_T2_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_T2_T1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


# 1. Define the contrast names to iterate over
# We start from the T2 contrast up to the T9 contrast
contrast_names <- c(
  "host_age_T2_vs_T1",
  "host_age_T3_vs_T1",
  "host_age_T4_vs_T1",
  "host_age_T5_vs_T1",
  "host_age_T6_vs_T1",
  "host_age_T7_vs_T1",
  "host_age_T8_vs_T1"
)

# 2. Initialize an empty list to store the results from each loop
all_timepoints_results <- list()

# 3. Loop through each contrast name
for (contrast in contrast_names) {
  
  # a. Extract the time point (e.g., "T2" from "host_age_T2_vs_T1")
  # We use a regular expression to extract the second time point
  timepoint_label <- sub("host_age_(T\\d+)_vs_T1", "\\1", contrast)
  
  # b. Get the DESeq2 results for the current contrast
  res_current <- results(DESEQ_calf_timepoints,
                         name = contrast,
                         tidy = TRUE)
  
  # c. Filter for significant ASVs (padj < 0.01 & |log2FoldChange| > 2)
  # and rename the row column to ASV
  sigASVs_current <- res_current %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
    dplyr::rename(ASV = row)
  
  # d. Get a vector of significant ASV names
  sigASVs_vec_current <- sigASVs_current %>% pull(ASV)
  
  # f. Subset the phyloseq object for T1 and the current time point (e.g., T2)
  # This part is needed to get the taxonomic information associated with the ASVs.
  # We extract the target timepoint from the contrast name
  
  # Get the two time points for subsetting (T1 and the current one)
  timepoints_to_subset <- c("T1", timepoint_label)
  
  T_pair <- subset_samples(calf_phyloseq_rare, host_age %in% timepoints_to_subset) %>%
    transform_sample_counts(function(x) x + 1)
  
  # g. Prune the phyloseq file to keep only the significant ASVs
  T_pair_DESeq <- prune_taxa(sigASVs_vec_current, T_pair)
  
  # h. Extract taxonomy, join with differential expression results, and prepare for plotting
  sigASVs_final <- tax_table(T_pair_DESeq) %>%
    as.data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    right_join(sigASVs_current, by = "ASV") %>%
    arrange(log2FoldChange) %>%
    # Use make.unique on Genus to ensure a unique label for plotting if needed
    mutate(Genus = make.unique(Genus)) %>%
    # Add the timepoint label column!
    mutate(host_age = timepoint_label)
  
  # i. Store the final data frame in the list
  all_timepoints_results[[timepoint_label]] <- sigASVs_final
  
}

# 4. Combine all the results into a single data frame
final_heatmap_data <- bind_rows(all_timepoints_results)

Tx_count <- final_heatmap_data %>%
  filter(host_age == "T8") %>%
  nrow()
Tx_count

# You may want to re-factor the Genus and host_age for better plotting order
final_heatmap_data <- final_heatmap_data %>%
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
  mutate(host_age = factor(host_age, levels = c("T2", "T3", "T4", "T5", "T6", "T7", "T8")))

ggplot(final_heatmap_data,aes(x = host_age, y = Family, fill = log2FoldChange)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "viridis",  # A good high-contrast sequential palette
    direction = -1,    # Reverse direction for higher FoldChange = darker/more intense color
    name = expression(bold(log2FoldChange)) # Label the legend
  )
