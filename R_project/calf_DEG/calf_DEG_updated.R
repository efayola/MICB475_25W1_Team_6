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
calf_phyloseq_final_genus_no_diet <- subset_samples(calf_phyloseq_final_genus, host_age != "not applicable")
nsamples(calf_phyloseq_final_genus_no_diet) #check sample number
calf_phyloseq_no_diet_and_T9 <- subset_samples(calf_phyloseq_final_genus_no_diet, host_age != "T9")
nsamples(calf_phyloseq_no_diet_and_T9)


#### DESeq sex ####
calf_plus1 <- transform_sample_counts(calf_phyloseq_no_diet_and_T9, function(x) x+1)
calf_sex_deseq <- phyloseq_to_deseq2(calf_plus1, ~`host_sex`)
calf_sex_deseq$host_sex <- relevel(calf_sex_deseq$host_sex, ref = "female")
DESEQ_calf_sex <- DESeq(calf_sex_deseq)

resultsNames(DESEQ_calf_sex)
#[1] "Intercept"               "host_sex_male_vs_female"

res_calf_sex <- results(DESEQ_calf_sex, name = "host_sex_male_vs_female", tidy = TRUE)
View(res_calf_sex)


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
    geom_vline(xintercept = c(-1.9, 1.9), linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkgrey") +
    # labs(title = title,
    #      x = expression(Log[2]* "Fold Change"),
    #      y = expression(-Log[10]*"(adj. p-value)")) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) 
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank())
vol_plot_sex_upd
