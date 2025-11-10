library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
install.packages("ggpubr")
#### Load in RData ####
load("calf_phyloseq_rare.RData")
load("calf_phyloseq_final.RData")

nsamples(calf_phyloseq_rare)

#### Alpha diversity ######
plot_richness(calf_phyloseq_rare) 
plot_richness(calf_phyloseq_rare, measures = c("Shannon","Chao1")) 
estimate_richness(calf_phyloseq_rare)


#subset out samples with diet & milk
calf_phyloseq_no_diet <- subset_samples(calf_phyloseq_rare, host_age != "not applicable")
nsamples(calf_phyloseq_no_diet) #check sample number


#plot alpha diversity metric for male & female in bulk
set.seed(42)
calf_gg_richness_sex <- plot_richness(calf_phyloseq_no_diet, x = "host_sex", measures = c("Shannon","Chao1")) +
  labs(x = expression(bold("Sample Category")), y = expression(bold("Alpha Diversity Measure")))+
  geom_point(aes(color = NULL), alpha = 0) + # Makes default points transparent
  geom_jitter(aes(color = host_sex), width = 0.2, size = 2, alpha = 0.8) +                  # disperses points
  geom_boxplot(alpha = 0) +                                                                 # transparent boxplot
  scale_color_manual(values = c("#1b9e77", "#d95f02")) +                                    # custom nice colors "#7570b3" "#66a61e"
  scale_x_discrete(limits = c("male", "female")) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5))

#remove overlapp dots from geom_points, keep disperses points
calf_gg_richness_sex_no_overlap <- calf_gg_richness_sex
calf_gg_richness_sex_no_overlap$layers <- calf_gg_richness_sex_no_overlap$layers[-(1:2)]
calf_gg_richness_sex_no_overlap

#save
ggsave(filename = "calf_gg_richness_sex_no_overlap.png"
       , calf_gg_richness_sex_no_overlap
       , height=4, width=6)

#added significance test
calf_gg_richness_sex_no_overlap_with_stats <- calf_gg_richness_sex_no_overlap + 
  stat_compare_means(
    comparisons = list(c("male", "female")),
    method = "wilcox.test",      # Use Wilcoxon test for 2 groups
    label = "p.signif"           # This will display asterisks (****)
  )
calf_gg_richness_sex_no_overlap_with_stats

#save
ggsave(filename = "calf_gg_richness_sex_no_overlap_with_stats.png"
       , calf_gg_richness_sex_no_overlap_with_stats
       , height=4, width=6)
