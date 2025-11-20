
BiocManager::install("DESeq2")
library(tidyverse)
library(phyloseq)
library(DESeq2)


#### Load data ####
load("calf_phyloseq_rare.RData")

#### DESeq ####
calf_deseq <- phyloseq_to_deseq2(calf_phyloseq_rare, ~`host_sex`)
DESEQ_calf <- DESeq(calf_deseq)

calf_plus1 <- transform_sample_counts(calf_phyloseq_rare, function(x) x+1)
calf_deseq <- phyloseq_to_deseq2(calf_plus1, ~`host_sex`)
DESEQ_calf <- DESeq(calf_deseq)
res <- results(DESEQ_calf, name = "host_sex_male_vs_female", tidy = TRUE)
View(res)

# 2. See exactly what result names DESeq2 created
resultsNames(DESEQ_calf)

# Look at results 

## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="vol_plot.png",vol_plot)
# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
calf_DESeq <- prune_taxa(sigASVs_vec,calf_phyloseq_rare)
sigASVs <- tax_table(calf_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
