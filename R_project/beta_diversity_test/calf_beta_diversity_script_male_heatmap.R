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
#Get the unique time points
timepoints <- unique(sample_data(calf_phyloseq_male)$host_age)
timepoints <- as.character(timepoints)
# Generate all unique combinations of two time points
pairwise_combinations <- combn(timepoints, 2, simplify = FALSE)

# Initialize an empty data frame to store all results with two comparison columns
all_permanova_results_split <- data.frame(
  TP1 = character(),     
  TP2 = character(),     
  Df = numeric(),
  R2 = numeric(),
  F_value = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# 2. Loop through every pairwise combination
for (pair in pairwise_combinations) {
  tp1 <- pair[1]
  tp2 <- pair[2]
  
  # --- Subsetting and Calculation ---
  
  # Subset the phyloseq object for the current pair of timepoints
  calf_phyloseq_subset <- subset_samples(calf_phyloseq_male, host_age %in% c(tp1, tp2))
  
  # Calculate Bray-Curtis distance for the subsetted data
  calf_dm_braycurtis_subset <- vegdist(t(otu_table(calf_phyloseq_subset)), method="bray")
  
  # Create sample data frame for the subsetted object (only host_age needed for adonis2)
  samp_dat_subset <- data.frame(sample_data(calf_phyloseq_subset))
  
  # --- Run PERMANOVA and Save Results ---
  
  # Run PERMANOVA (adonis2)
  permanova_result <- adonis2(calf_dm_braycurtis_subset ~ host_age, data=samp_dat_subset)
  
  # Extract the results for the model (the first row)
  model_stats <- permanova_result[1, ]
  
  # Create a data frame for the current comparison with split columns
  current_df <- data.frame(
    TP1 = tp1,                   # Time Point 1
    TP2 = tp2,                   # Time Point 2
    Df = model_stats$Df,
    R2 = model_stats$R2,
    F_value = model_stats$F,
    P_value = model_stats$`Pr(>F)`,
    stringsAsFactors = FALSE
  )
  
  # Append the results to the main results data frame
  all_permanova_results_split <- rbind(all_permanova_results_split, current_df)
}

#make duplicate data
df_all_permanova_results_split <- all_permanova_results_split %>% rename(TP1 = TP2,
                                                                         TP2 = TP1)
all_permanova_results_merged<- rbind(all_permanova_results_split, df_all_permanova_results_split) %>% 
  arrange(TP1, TP2) %>%
  filter(TP1 < TP2)

#Remove Timeline 9
all_permanova_results_merged_no_T9<- rbind(all_permanova_results_split, df_all_permanova_results_split) %>% 
  arrange(TP1, TP2) %>%
  filter(TP1 < TP2) %>%
  filter(TP2 != "T9")

#generate lower triangle heatmap
calf_beta_diversity_heatmap_male <- ggplot(all_permanova_results_merged_no_T9, aes(x = TP1, y = TP2, fill = R2)) +
  
  # Add tiles (the squares) and map the R2 value to the fill color
  geom_tile() +
  # You can adjust the number of digits shown
  geom_text(aes(label = round(R2, 3)), color = "black", size = 3) +
  
  # Use a color palette that clearly shows the magnitude of R2 (e.g., Viridis or sequential)
  scale_fill_viridis_c(
    option = "viridis",  # A good high-contrast sequential palette
    direction = -1,    # Reverse direction for higher R2 = darker/more intense color
    name = expression(bold(R^2~Value)) # Label the legend with R^2
  ) + 
# Ensure the axes cover all the data and don't stretch the plot unnaturally
  scale_x_discrete(limits = c("T1", "T2", "T3", "T4", "T5", "T6", "T7","T8"),
                   expand = c(0, 0.5)) +
  #reverse y-axis label
  scale_y_discrete(limits = c("T8", "T7", "T6", "T5", "T4", "T3", "T2", "T1"),
                   expand = c(0, 0.5)) +
  theme(
    # Rotate axis labels for better readability if time points are long
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    #remove Axis Title
    axis.title = element_blank(),
    # Make the plot area square and centered
    aspect.ratio = 1,
    #remove background grid and color
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  )

calf_beta_diversity_heatmap_male

ggsave("beta_diversity_test/calf_beta_diversity_heatmap_male.png"
       , calf_beta_diversity_heatmap_male
       , height=4, width=5)
