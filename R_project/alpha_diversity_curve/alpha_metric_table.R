library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)
library(writexl)

#### Load in RData ####
load("all_phyloseq_objects/calf_phyloseq_rare.RData")
load("all_phyloseq_objects/calf_phyloseq_final.RData")

nsamples(calf_phyloseq_rare)

#### Alpha diversity ######
plot_richness(calf_phyloseq_rare) 
plot_richness(calf_phyloseq_rare, measures = c("Shannon","Chao1")) 
estimate_richness(calf_phyloseq_rare)


#subset out samples with diet & milk
calf_phyloseq_no_diet <- subset_samples(calf_phyloseq_rare, host_age != "not applicable")
nsamples(calf_phyloseq_no_diet) #check sample number
calf_phyloseq_no_diet_and_T9 <- subset_samples(calf_phyloseq_no_diet, host_age != "T9")
nsamples(calf_phyloseq_no_diet_and_T9)

set.seed(42)
alpha_diversity_data <- estimate_richness(calf_phyloseq_no_diet_and_T9, measures = c("Shannon", "Chao1"))
sample_metadata <- sample_data(calf_phyloseq_no_diet_and_T9)
alpha_df <- data.frame(alpha_diversity_data, sample_metadata)

my_comparisons <- list(c("T1", "T2"), c("T1", "T3"), c("T1", "T4"), c("T1", "T5"), c("T1", "T6"), c("T1", "T7"), c("T1", "T8"),
                       c("T2", "T3"), c("T2", "T4"), c("T2", "T5"), c("T2", "T6"), c("T2", "T7"), c("T2", "T8"),
                       c("T3", "T4"), c("T3", "T5"), c("T3", "T6"), c("T3", "T7"), c("T3", "T8"),
                       c("T4", "T5"), c("T4", "T6"), c("T4", "T7"), c("T4", "T8"),
                       c("T5", "T6"), c("T5", "T7"), c("T5", "T8"),
                       c("T6", "T7"), c("T6", "T8"),
                       c("T7", "T8"))

#### Shannon DATA Statistics ####
shannon_stats_table <- compare_means(
  Shannon ~ host_age,          # Formula: Shannon diversity vs. host_age
  data = alpha_df,
  method = "wilcox.test",     # Wilcoxon rank-sum test
  group.by = NULL,            # No grouping needed
  comparisons = my_comparisons # Use the specific comparisons list
)

#duplicate data to have nice arranged comparison
shannon_stats_table <- shannon_stats_table %>% rename(group1 = 'group1',
                                                      group2 = 'group2')

shannon_stats_table_d1 <- shannon_stats_table %>% rename(group1 = 'group2',
                                                         group2 = 'group1')
shannon_stats_table_all<- rbind(shannon_stats_table, shannon_stats_table_d1) %>% 
  arrange(group1, group2) %>%
  filter(group1 < group2) %>%
  select(.y.,group1, group2, p, p.adj, p.signif)

write_xlsx(shannon_stats_table_all, "alpha_diversity_curve/shannon_timepoints.xlsx")





#### CHAO DATA Statistics ####
chao1_stats_table <- compare_means(
  Chao1 ~ host_age,            # Formula: Chao1 richness vs. host_age
  data = alpha_df,
  method = "wilcox.test",       # Wilcoxon rank-sum test
  group.by = NULL,
  comparisons = my_comparisons
)

chao1_stats_table <- chao1_stats_table %>% rename(group1 = 'group1',
                                                      group2 = 'group2')

chao1_stats_table_d1 <- chao1_stats_table %>% rename(group1 = 'group2',
                                                         group2 = 'group1')
chao1_stats_table_all<- rbind(chao1_stats_table, chao1_stats_table_d1) %>% 
  arrange(group1, group2) %>%
  filter(group1 < group2) %>%
  select(.y.,group1, group2, p, p.adj, p.signif)

write_xlsx(chao1_stats_table_all, "alpha_diversity_curve/chao1_timepoints.xlsx")

#### Shannon between sex at different timepoints
timepoints_to_include <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8")
alpha_df_filtered <- alpha_df %>%
  dplyr::filter(host_age %in% timepoints_to_include) %>%
  dplyr::mutate(host_sex = as.factor(host_sex),
                host_age = factor(host_age, levels = timepoints_to_include))

shannon_stats_table_sex <- compare_means(
  Shannon ~ host_sex,           # Compare Shannon by host_sex
  data = alpha_df_filtered,
  method = "wilcox.test",      # Use Wilcoxon test
  group.by = "host_age"       # Group by time point
) %>%
  select(.y.,host_age,group1, group2, p, p.adj, p.signif)

write_xlsx(shannon_stats_table_sex, "alpha_diversity_curve/shannon_sex_at_timepoints.xlsx")

chao1_stats_table_sex <- compare_means(
  Chao1 ~ host_sex,           # Compare chao1 by host_sex
  data = alpha_df_filtered,
  method = "wilcox.test",      # Use Wilcoxon test
  group.by = "host_age"       # Group by time point
) %>%
  select(.y.,host_age,group1, group2, p, p.adj, p.signif)

write_xlsx(chao1_stats_table_sex, "alpha_diversity_curve/chao1_sex_at_timepoints.xlsx")
