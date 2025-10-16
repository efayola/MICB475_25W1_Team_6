## Agenda
* HW for the week: proceed with processing data, **demultiplexing, denoising**, diversity metrics on QIIME2
 * ask Reetu/Evelyn after trimming and before proceeding with next steps
* Project: Tracking fecal microbiome over cow development time points.
* For proposal: explain why network analysis
* run it as a single read, instead of paired because data is too big
* **NOTE:** run the single read and paired read demux, and then try out denoising for paired only first

## Notes
* What to do with the data
  * filter out T9
  * filter out the NAs
* start demultiplexing ASAP, each set may be overnight
* separate by host sex - do the same processes in parallel - do so after creating phyloseq object (in R)

### Aim
1. boutique analysis - diversity analysis as line graph/box plot
   * trends over time
   * timepoint develop
   * diff btwn male and female
   * before vs after weening

choose top 3 interesting timepoints 

2. diversity metrics - box plot for the different time points 
3. core microbiome venn diagram 
4. differential abundance - are there shared microbes and are they fluctuating?
5. network analysis 

