library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(readxl)

#### Load data ####
# Change file paths as necessary
metafp <- "cattle_metadata.xlsx"
meta <- read_excel(metafp)

otufp <- "single_table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "single_taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "single_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'Run'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
calf_phyloseq <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(calf_phyloseq)
sample_data(calf_phyloseq)
tax_table(calf_phyloseq)
phy_tree(calf_phyloseq)

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
calf_phyloseq_filt <- subset_taxa(calf_phyloseq,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
calf_phyloseq_filt_nolow <- filter_taxa(calf_phyloseq_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
calf_phyloseq_filt_nolow_samps <- prune_samples(sample_sums(calf_phyloseq_filt_nolow)>100, calf_phyloseq_filt_nolow)
# Remove samples where month is na
calf_phyloseq_final <- subset_samples(calf_phyloseq_filt_nolow_samps, !is.na(month) )

sample_variables(calf_phyloseq_final)
get_variable(calf_phyloseq_final, c("env_medium","host_age","host_sex"))
nsamples(calf_phyloseq_final)

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size

#Generate refraction curve
rarecurve(t(as.data.frame(otu_table(calf_phyloseq_final))), cex=0.1)

depths <- sample_sums(calf_phyloseq_final)
specific_depth <- depths["SRR29857246"]
specific_sampling_depth <- depths['SRR29857208']

calf_phyloseq_rare <- rarefy_even_depth(calf_phyloseq_final, rngseed = 1, sample.size = 32000)


##### Saving #####
save(calf_phyloseq_final, file="calf_phyloseq_final.RData")
#save(calf_phyloseq_rare, file="calf_phyloseq_rare.RData")