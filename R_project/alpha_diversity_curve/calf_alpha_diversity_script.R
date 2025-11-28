library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggpubr)
library(rstatix)

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

#plot alpha diversity metric for male & female in bulk w/o diet & T9
set.seed(42)
calf_gg_richness_sex <- plot_richness(calf_phyloseq_no_diet_and_T9, x = "host_sex", measures = c("Shannon","Chao1")) +
  labs(x = expression(bold("Sample Category")), y = expression(bold("Alpha Diversity Measure")))+
  geom_point(aes(color = NULL), alpha = 0) + # Makes default points transparent
  geom_jitter(aes(color = host_sex), width = 0.2, size = 2, alpha = 0.8) +                  # disperses points
  geom_boxplot(alpha = 0) +                                                                 # transparent boxplot
  # scale_color_manual(values = c("#1b9e77", "#d95f02")) +                                    # custom nice colors "#7570b3" "#66a61e"
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









#Compare between male and female at T6 ----------
#subset out samples with T6
calf_phyloseq_T6 <- subset_samples(calf_phyloseq_rare, host_age == "T6")
nsamples(calf_phyloseq_T6) #check sample number

#plot alpha diversity metric for male & female in bulk
set.seed(42)
calf_gg_richness_T6 <- plot_richness(calf_phyloseq_T6, x = "host_sex", measures = c("Shannon","Chao1")) +
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
calf_gg_richness_T6

#remove overlapp dots from geom_points, keep disperses points
calf_gg_richness_T6_no_overlap <- calf_gg_richness_T6
calf_gg_richness_T6_no_overlap$layers <- calf_gg_richness_T6_no_overlap$layers[-(1:2)]
calf_gg_richness_T6_no_overlap

#added significance test
calf_gg_richness_T6_no_overlap_with_stats <- calf_gg_richness_T6_no_overlap + 
  stat_compare_means(
    comparisons = list(c("male", "female")),
    method = "wilcox.test",      # Use Wilcoxon test for 2 groups
    label = "p.signif"           # This will display asterisks (****)
  )
calf_gg_richness_T6_no_overlap_with_stats









# compare between timeline in bulk--------------
set.seed(42)
calf_gg_richness_timepoints <- plot_richness(calf_phyloseq_no_diet, x = "host_age", measures = c("Shannon","Chao1")) +
  labs(x = expression(bold("Sample Category")), y = expression(bold("Alpha Diversity Measure"))) +
  geom_point(aes(color = NULL), alpha = 0) + # Makes default points transparent
  geom_jitter(aes(color = host_age),size = 2, alpha = 0.8) +                  # disperses points
  geom_boxplot(alpha = 0, position = position_dodge(width = 2)) +                 # transparent boxplot

  #scale_color_manual(name = "Host Sex",values = c("male" = "#1b9e77", "female" = "#d95f02"))+           # custom nice colors "#7570b3" "#66a61e"
  scale_x_discrete(limits = c("T1","T2","T3","T4","T5","T6","T7","T8","T9")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5))
calf_gg_richness_timepoints


#remove overlapp dots from geom_points, keep disperses points
calf_gg_richness_timepoints_no_overlap <- calf_gg_richness_timepoints
calf_gg_richness_timepoints_no_overlap$layers <- calf_gg_richness_timepoints_no_overlap$layers[-(1:2)]
calf_gg_richness_timepoints_no_overlap

ggsave(filename = "calf_gg_richness_timepoints_no_overlap.png"
       , calf_gg_richness_timepoints_no_overlap
       , height=4, width=6)

#added significance test
calf_gg_richness_timepoints_no_overlap_with_stats <- calf_gg_richness_timepoints_no_overlap + 
  stat_compare_means(
    comparisons = list(c("T1", "T2"),
                       c("T1", "T3"),
                       c("T1", "T4"),
                       c("T1", "T5"),
                       c("T1", "T6"),
                       c("T1", "T7"),
                       c("T1", "T8"),
                       c("T1", "T9"),
                       c("T2", "T3"),
                       c("T2", "T4"),
                       c("T2", "T5"),
                       c("T2", "T6"),
                       c("T2", "T7"),
                       c("T2", "T8"),
                       c("T2", "T9")
                       ),
    method = "wilcox.test",      # Use Wilcoxon test for 2 groups
    label = "p.signif"           # This will display asterisks (****)
  )
calf_gg_richness_timepoints_no_overlap_with_stats

ggsave(filename = "calf_gg_richness_timepoints_no_overlap_with_stats_T1_T2.png"
       , calf_gg_richness_timepoints_no_overlap_with_stats
       , height=4, width=6)










##### TimePoint with sex --------------
calf_gg_richness_timepoints_sex <- plot_richness(calf_phyloseq_no_diet, x = "host_sex", measures = c("Shannon","Chao1")) +
  labs(x = expression(bold("Sample Category")), y = expression(bold("Alpha Diversity Measure"))) +
  geom_point(aes(color = NULL), alpha = 0) + # Makes default points transparent
  geom_jitter(aes(color = host_sex),
              position = position_jitterdodge(jitter.width = 0.2,   # Your original jitter width
                                              dodge.width = 2),  # Default dodge width, adjust if needed
              size = 2, alpha = 0.8) +                  # disperses points
  geom_boxplot(aes(fill = host_sex), alpha = 0, position = position_dodge(width = 5)) +                 # transparent boxplot
  facet_grid(
    variable ~ host_age, 
    scales = "free_y"
  ) +
  #scale_color_manual(name = "Host Sex",values = c("male" = "#1b9e77", "female" = "#d95f02"))+           # custom nice colors "#7570b3" "#66a61e"
  scale_x_discrete(limits = c("male","female")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5))

#remove overlapp dots from geom_points, keep disperses points
calf_gg_richness_timepoints_sex_no_overlap <- calf_gg_richness_timepoints_sex
calf_gg_richness_timepoints_sex_no_overlap$layers <- calf_gg_richness_timepoints_sex_no_overlap$layers[-(1:2)]
calf_gg_richness_timepoints_sex_no_overlap
ggsave(filename = "calf_gg_richness_timepoints_sex_no_overlap.png"
       , calf_gg_richness_timepoints_sex_no_overlap
       , height=4, width=6)


#added significance test
calf_gg_richness_timepoints_sex_no_overlap_with_stats <- calf_gg_richness_timepoints_sex_no_overlap + 
  stat_compare_means(
    comparisons = list(c("male", "female")),
    method = "wilcox.test",      # Use Wilcoxon test for 2 groups
    label = "p.signif",           # This will display asterisks (****)
    method.args = list(exact = FALSE),
    vjust = 0.5,
    y.position = -1
  )
calf_gg_richness_timepoints_sex_no_overlap_with_stats

ggsave(filename = "calf_gg_richness_timepoints_sex_no_overlap_with_stats.png"
       , calf_gg_richness_timepoints_sex_no_overlap_with_stats
       , height=8, width=12)
