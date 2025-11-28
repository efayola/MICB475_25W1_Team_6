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

set.seed(42)
#### Load data ####
# load("all_phyloseq_objects/calf_phyloseq_final.RData")
# rank_names(calf_phyloseq_final)
# 
# #merges species that have the same taxonomy at genus level
# calf_phyloseq_final_genus <- tax_glom(calf_phyloseq_final, taxrank = "Genus")
# #save(calf_phyloseq_final_genus, file="all_phyloseq_objects/calf_phyloseq_final_genus.RData")

load("all_phyloseq_objects/calf_phyloseq_final_genus.RData")
nsamples(calf_phyloseq_final_genus)

#subset out samples with diet & milk & T9
calf_phyloseq_no_diet_and_T9 <- subset_samples(calf_phyloseq_final_genus, host_age != "not applicable" & host_age !="T9")
nsamples(calf_phyloseq_no_diet_and_T9) #check sample number


#### DESeq sex ####
calf_plus1 <- transform_sample_counts(calf_phyloseq_no_diet_and_T9, function(x) x+1)
calf_sex_deseq <- phyloseq_to_deseq2(calf_plus1, ~`host_sex`)
calf_sex_deseq$host_sex <- relevel(calf_sex_deseq$host_sex, ref = "female")
DESEQ_calf_sex <- DESeq(calf_sex_deseq)

resultsNames(DESEQ_calf_sex)
#[1] "Intercept"               "host_sex_male_vs_female"

res_calf_sex <- results(DESEQ_calf_sex, name = "host_sex_male_vs_female", tidy = TRUE)
View(res_calf_sex)

# To get table of results
sigASVs_sex <- res_calf_sex %>% 
  filter(padj<0.01 & abs(log2FoldChange)>1.9) %>%
  dplyr::rename(ASV=row)
View(sigASVs_sex)
#To get ASV
sigASVs_sex_vec <- sigASVs_sex %>%
  pull(ASV)

# Prune phyloseq file
male_female_DESeq <- prune_taxa(sigASVs_sex_vec,calf_phyloseq_no_diet_and_T9)
sigASVs_sex <- tax_table(male_female_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_sex) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
View(sigASVs_sex)

#graph volcano plot
vol_plot_sex_upd <- res_calf_sex %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>1.9) %>%
  ggplot() +
    geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant),size = 2, alpha = 0.8) + 
    scale_color_manual(
      values = c("#A9A9A9", "#F8766D"),     # blue = not sig, red = sig
      labels = c("Not significant", "Significant"),
      name   = "Significance",
    ) + 
    geom_label_repel(
      data = sigASVs_sex, # Use the filtered data for only the points you want to label
      aes(x = log2FoldChange, y = -log10(padj), label = Genus),
      size = 3.5,
      box.padding = 0.9,
      point.padding = 0.5,
    ) +
    geom_vline(xintercept = c(-1.9, 1.9), linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgrey") +
    labs(title = "M vs F",
         x = expression(Log[2]* "Fold Change"),
         y = expression(-Log[10]*"(adj. p-value)")) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
vol_plot_sex_upd
ggsave("calf_DEG/labelled_graph/vol_plot_sex_upd.png", width = 8, height = 6.5, dpi = 300, bg = "white")





#### DESeq T6 vs T1 ####
calf_plus1_T6_T1 <- transform_sample_counts(calf_phyloseq_no_diet_and_T9, function(x) x+1)
calf_sex_deseq_T6_T1 <- phyloseq_to_deseq2(calf_plus1_T6_T1, ~`host_age`)
calf_sex_deseq_T6_T1$host_age <- relevel(calf_sex_deseq_T6_T1$host_age, ref = "T1")
DESEQ_calf_T6_T1<- DESeq(calf_sex_deseq_T6_T1)

resultsNames(DESEQ_calf_T6_T1)
# [1] "Intercept"         "host_age_T2_vs_T1" "host_age_T3_vs_T1" "host_age_T4_vs_T1" "host_age_T5_vs_T1"
# [6] "host_age_T6_vs_T1" "host_age_T7_vs_T1" "host_age_T8_vs_T1"

res_calf_T6_T1 <- results(DESEQ_calf_T6_T1, name = "host_age_T6_vs_T1", tidy = TRUE)
View(res_calf_T6_T1)

#graph volcano plot
vol_plot_T6_T1_upd <- res_calf_T6_T1 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant),size = 2, alpha = 0.8) + 
  scale_color_manual(
    values = c("#A9A9A9", "#F8766D"),     # blue = not sig, red = sig
    labels = c("Not significant", "Significant"),
    name   = "Significance",
  ) + 
  # geom_label_repel(
  #   data = sigASVs_sex, # Use the filtered data for only the points you want to label
  #   aes(x = log2FoldChange, y = -log10(padj), label = Genus),
  #   size = 3.5,
  #   box.padding = 0.9,
  #   point.padding = 0.5,
  # ) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgrey") +
  labs(title = "T6 vs T1",
       x = expression(Log[2]* "Fold Change"),
       y = expression(-Log[10]*"(adj. p-value)")) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
vol_plot_T6_T1_upd
