library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(vegan)

set.seed(42)
##### Male Beta Diversity Heat Map
#all possibility for male
samp_dat_wdiv_male <- data.frame(sample_data(calf_phyloseq_male), estimate_richness(calf_phyloseq_male))
#adonis2(calf_dm_braycurtis_male ~ host_age, data=samp_dat_wdiv_male)
#order data
samp_dat_wdiv_male$host_age <- factor(samp_dat_wdiv_male$host_age, 
                                      levels = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9"))
#get unique timepoints
timepoints <- levels(samp_dat_wdiv_male$host_age)
# Initialize a list to store the PERMANOVA results
pairwise_permanova_results <- list()

#Use combn to generate all unique pairs and loop through them
combn(timepoints, 2, simplify = FALSE, FUN = function(pair) {
  # a. Filter the metadata (samp_dat_wdiv_male) for the current pair
  subset_data <- subset(samp_dat_wdiv_male, host_age %in% pair)
  # b. Get the sample names for filtering the distance matrix
  sample_names <- rownames(subset_data)
  # c. Subset the Bray-Curtis distance matrix
  # You need to subset the distance matrix to include only the samples in the current pair.
  # The distance matrix must be subsetted on both rows and columns.
  subset_dm <- as.matrix(calf_dm_braycurtis_male)[sample_names, sample_names]
  # d. Run the PERMANOVA (adonis2) on the subset data
  permanova_result <- adonis2(subset_dm ~ host_age, data = subset_data)
  # e. Store the results in the list
  pair_name <- paste(pair[1], "vs", pair[2])
  pairwise_permanova_results[[pair_name]] <<- permanova_result
})
print(pairwise_permanova_results[["T8 vs T9"]])