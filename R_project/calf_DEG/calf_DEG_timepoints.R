
BiocManager::install("DESeq2")
library(tidyverse)
library(phyloseq)
library(DESeq2)


#### Load data ####
load("calf_phyloseq_final.RData")
tax_glom(calf_phyloseq_final, taxrank = "Family")

#### DESeq sex ####
calf_plus1 <- transform_sample_counts(calf_phyloseq_final, function(x) x+1)
calf_sex_deseq <- phyloseq_to_deseq2(calf_plus1, ~`host_sex`)
DESEQ_calf_sex <- DESeq(calf_sex_deseq)
res_calf_sex <- results(DESEQ_calf_sex, name = "host_sex_male_vs_female", tidy = TRUE)
View(res_calf_sex)

# 2. See exactly what result names DESeq2 created
resultsNames(DESEQ_calf_sex)

# Look at results 

## Volcano plot: effect size VS significance
ggplot(res_calf_sex) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot_sex <- res_calf_sex %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="calf_sex_vol_plot.png",vol_plot_sex)
# To get table of results
sigASVs_calf_sex <- res_calf_sex %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_calf_sex)
# Get only asv names
sigASVs_calf_sex_vec <- sigASVs_calf_sex %>%
  pull(ASV)

# Prune phyloseq file
calf_sex_DESeq <- prune_taxa(sigASVs_calf_sex_vec,calf_phyloseq_rare)
sigASVs_calf_sex <- tax_table(calf_sex_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_calf_sex) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggsave("DEG_sex_bar.png",
       ggplot(sigASVs_calf_sex) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)))

#### T1 and T5 ####
T1_ref  <- "T1"
T5_test <- "T5"

T1_T5_pair <- subset_samples(calf_phyloseq_rare, host_age %in% c(T1_ref, T5_test)) %>%
  transform_sample_counts(function(x) x + 1)

calf_deseq_T1_T5 <- phyloseq_to_deseq2(T1_T5_pair, ~`host_age`)
calf_deseq_T1_T5$host_age <- relevel(factor(calf_deseq_T1_T5$host_age), ref = T1_ref)
calf_deseq_T1_T5 <- DESeq(calf_deseq_T1_T5)

res_T1_T5 <- results(calf_deseq_T1_T5, contrast = c("host_age", T5_test, T1_ref), tidy = TRUE)

vol_plot_T1_T5 <- res_T1_T5 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))


ggsave(filename="DEG_T1_T5_vol_plot.png",vol_plot_T1_T5)


# To get table of results
sig_T1_T5_ASVs <- res_T1_T5 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sig_T1_T5_ASVs)

# Get only asv names
sigASVs_T1_T5_vec <- sig_T1_T5_ASVs %>%
  pull(ASV)

# Prune phyloseq file
calf_T1_T5_DESeq <- prune_taxa(sigASVs_T1_T5_vec,T1_T5_pair)
sig_T1_T5_ASVs <- tax_table(calf_T1_T5_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sig_T1_T5_ASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggsave("DEG_T1_T5_bar.png", 
       ggplot(sig_T1_T5_ASVs) +
         geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
         geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
         theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
)

# Family
sig_T1_T5_family_ASVs <- tax_table(calf_T1_T5_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sig_T1_T5_ASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Family = make.unique(Family)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

ggsave("DEG_T1_T5_family_bar.png", 
       ggplot(sig_T1_T5_ASVs) +
         geom_bar(aes(x = Family, y = log2FoldChange), stat = "identity") +
         geom_errorbar(aes(x = Family,
                           ymin = log2FoldChange - lfcSE,
                           ymax = log2FoldChange + lfcSE)) +
         theme(axis.text.x = element_text(angle = 90,
                                          hjust = 1,
                                          vjust = 0.5))
)

#### T5 and T8 ####

T5_ref  <- "T5"
T8_test <- "T8"

T5_T8_pair <- subset_samples(calf_phyloseq_rare, host_age %in% c(T5_ref, T8_test)) %>%
  transform_sample_counts(function(x) x + 1)

calf_deseq_T5_T8 <- phyloseq_to_deseq2(T5_T8_pair, ~`host_age`)
calf_deseq_T5_T8$host_age <- relevel(factor(calf_deseq_T5_T8$host_age), ref = T5_ref)
calf_deseq_T5_T8 <- DESeq(calf_deseq_T5_T8)

res_T5_T8 <- results(calf_deseq_T5_T8, contrast = c("host_age", T8_test, T5_ref), tidy = TRUE)

vol_plot_T5_T8 <- res_T5_T8 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="DEG_T5_T8_vol_plot.png",vol_plot_T5_T8)

# To get table of results
sig_T5_T8_ASVs <- res_T5_T8 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sig_T5_T8_ASVs)

# Get only asv names
sigASVs_T5_T8_vec <- sig_T5_T8_ASVs %>%
  pull(ASV)

# Prune phyloseq file
calf_T5_T8_DESeq <- prune_taxa(sigASVs_T5_T8_vec,T5_T8_pair)
sig_T5_T8_ASVs <- tax_table(calf_T5_T8_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sig_T5_T8_ASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggsave("DEG_T5_T8_bar.png", 
       ggplot(sig_T5_T8_ASVs) +
         geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
         geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
         theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
)

# Family
sig_T5_T8_family_ASVs <- tax_table(calf_T1_T5_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sig_T5_T8_ASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Family = make.unique(Family)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

ggsave("DEG_T5_T8_family_bar.png", 
       ggplot(sig_T5_T8_ASVs) +
         geom_bar(aes(x = Family, y = log2FoldChange), stat = "identity") +
         geom_errorbar(aes(x = Family,
                           ymin = log2FoldChange - lfcSE,
                           ymax = log2FoldChange + lfcSE)) +
         theme(axis.text.x = element_text(angle = 90,
                                          hjust = 1,
                                          vjust = 0.5))
)

#### T1 and T8 ####

T1_ref  <- "T1"
T8_test <- "T8"

T1_T8_pair <- subset_samples(calf_phyloseq_rare, host_age %in% c(T1_ref, T8_test)) %>%
  transform_sample_counts(function(x) x + 1)

calf_deseq_T1_T8 <- phyloseq_to_deseq2(T1_T8_pair, ~`host_age`)
calf_deseq_T1_T8$host_age <- relevel(factor(calf_deseq_T1_T8$host_age), ref = T1_ref)
calf_deseq_T1_T8 <- DESeq(calf_deseq_T1_T8)

res_T1_T8 <- results(calf_deseq_T1_T8, contrast = c("host_age", T8_test, T1_ref), tidy = TRUE)

vol_plot_T1_T8 <- res_T1_T8 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="DEG_T1_T8_vol_plot.png",vol_plot_T5_T8)

#### Putting volcano plots together
install.packages("patchwork")
library(patchwork)
library(dplyr)
library(ggplot2)

make_volcano <- function(res, title) {
  res %>%
    mutate(sig = case_when(
      is.na(padj) | is.na(log2FoldChange)                    ~ "NA",
      padj < 0.05 & abs(log2FoldChange) > 1                  ~ "Significant",
      TRUE                                                   ~ "Not significant"
    ),
    sig = factor(sig, levels = c("Not significant", "Significant", "NA"))   # ← critical line
    ) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c("#4A7BB7", "#DD3D2D", "black"),     # blue = not sig, red = sig, black = NA
      labels = c("Not significant", "Significant", "NA"),
      name   = "padj < 0.05 & |log2FC| > 1",
      drop   = FALSE                                 # ← forces all 3 levels even if missing
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
    labs(title = title,
         x = expression(Log[2]*" Fold Change"),
         y = expression(-Log[10]*"(adj. p-value)")) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# Re-run
p1 <- make_volcano(res_T1_T5, "T5 vs T1")
p2 <- make_volcano(res_T5_T8, "T8 vs T5")
p3 <- make_volcano(res_calf_sex, "M vs F")

# Now this will give ONE clean shared legend
p1 + p2 + p3 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("vol_T5vsT1_T8vsT5_MvsF_final.png", width = 18, height = 7, dpi = 300, bg = "white")

p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("vol_T5vsT1_T8vsT5_final.png", width = 18, height = 7, dpi = 300, bg = "white")
