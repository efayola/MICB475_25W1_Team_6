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
library(viridis)

#### Load Data ####
# load("all_phyloseq_objects/calf_phyloseq_rare.RData")
# rank_names(calf_phyloseq_rare)
# 
# #merges species that have the same taxonomy at genus level
# calf_phyloseq_rare_phylum <- tax_glom(calf_phyloseq_rare, taxrank = "Phylum")
# save(calf_phyloseq_rare_phylum, file="all_phyloseq_objects/calf_phyloseq_rare_phylum.RData")

load("all_phyloseq_objects/calf_phyloseq_rare_phylum.RData")
# nsamples(calf_phyloseq_rare_phylum) #480 removed 3 by rarefractino

calf_phyloseq_no_diet_and_T9 <- subset_samples(calf_phyloseq_rare_phylum, host_age != "not applicable" & host_age !="T9")
# nsamples(calf_phyloseq_no_diet_and_T9)

# Convert to relative abundance
calf_RA <- transform_sample_counts(calf_phyloseq_no_diet_and_T9, function(x) x/sum(x))

# Code to glom later by families
# # To remove black bars, "glom" by phylum first
# calf_phylum <- tax_glom(calf_RA, taxrank = "Phylum", NArm=FALSE)

#Aggregate relative abundances by Phylum for Top 10 Phylum
top10_taxa <- names(sort(taxa_sums(calf_RA), TRUE)[1:10]) 
calf_RA_top10 <- prune_taxa(top10_taxa, calf_RA)

# Other Phylum identification
other_taxa <- taxa_names(calf_RA)[!taxa_names(calf_RA) %in% top10_taxa]
calf_RA_grouped <- calf_RA
taxa_info <- as.data.frame(tax_table(calf_RA_grouped)) # Get the current tax_table
taxa_info[other_taxa, "Phylum"] <- "Other" # Rename 'other' phylum taxa to "Other"
tax_table(calf_RA_grouped) <- tax_table(as.matrix(taxa_info)) #update tax table
calf_RA_final <- tax_glom(calf_RA_grouped, taxrank = "Phylum") # reaggregate

# gg_taxa <- plot_bar(calf_RA_top10, fill="Phylum") + 
#   facet_wrap(.~host_age, scales = "free_x") + 
#   labs(y = "Relative Abundance") +
#   theme(
#     axis.text.x = element_blank(), # Removes the text labels
#     axis.ticks.x = element_blank() # Removes the tick marks/lines
#   )
# gg_taxa

gg_taxa <- ggplot(psmelt(calf_RA_final), aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    color = NA #remove bar black outline
  ) +
  facet_wrap(~host_age, scales = "free_x") +
  scale_y_continuous(
    # Specify the exact positions for the tick marks
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    # Specify the labels corresponding to those positions
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) + 
  scale_fill_viridis_d(option = "D") + 
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),# Removes the text labels
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_text(
      face = "bold", #bolding subheadings
      size = 10
    ),
    strip.background = element_rect(
      fill = "white", # Set the interior color to white
      colour = NA
    )
  ) +
  labs(y = "Relative Abundance (%)")
gg_taxa

ggsave("tax_bar_plot/plot_taxonomy.png"
       , gg_taxa
       , height=5.31, width =16)
