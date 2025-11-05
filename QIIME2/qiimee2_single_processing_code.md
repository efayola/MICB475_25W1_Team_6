### Screen session
```screen -r team6_single```

### Import manifest and demultiplex
```
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/calf/manifest_fwd_only.txt \
  --output-path ./single_demux_seqs.qza
```

### Create visualization of demultiplexed samples
```
qiime demux summarize \
  --i-data single_demux_seqs.qza \
  --o-visualization single_demux.qzv
```
### Download
```
scp root@10.19.139.156:/data/team6_calf/single_demux.qzv .
```

### Determine ASVs with DADA2
```
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single_demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 228 \
  --o-representative-sequences single_rep-seqs.qza \
  --o-table single_table.qza \
  --o-denoising-stats single_stats.qza
```

### Visualize ASVs stats
```
qiime feature-table summarize \
  --i-table single_table.qza \
  --o-visualization single_table.qzv \
  --m-sample-metadata-file /datasets/project_2/calf/metadata.txt

qiime feature-table tabulate-seqs \
  --i-data single_rep-seqs.qza \
  --o-visualization single_rep-seqs.qzv
```

### Train classifier 
https://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf
```
qiime feature-classifier extract-reads \
  --i-sequences /datasets/classifiers/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 228 \
  --o-reads ref-seqs-trimmed.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /datasets/classifiers/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza
```

### Taxanomic analysis
```
qiime feature-classifier classify-sklearn \
--i-classifier classifier.qza \
--i-reads single_rep-seqs.qza \
--o-classification single_taxonomy.qza 

qiime metadata tabulate \
  --m-input-file single_taxonomy.qza \
  --o-visualization single_taxonomy.qzv

# Taxonomy barplots
qiime taxa barplot \
  --i-table single_table.qza  \
  --i-taxonomy single_taxonomy.qza \
  --m-metadata-file /datasets/project_2/calf/metadata.txt \
  --o-visualization taxa-bar-plots.qzv

qiime taxa filter-table \
  --i-table single_table.qza \
  --i-taxonomy single_taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /datasets/project_2/calf/metadata.txt

qiime feature-table filter-samples \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file /datasets/project_2/calf/metadata.txt \
  --p-where "[host_age] != 'T9'" \
  --o-filtered-table table-no-T9-no-mit-no-chlor.qza

qiime feature-table summarize \
  --i-table table-no-T9-no-mit-no-chlor.qza \
  --o-visualization table-no-T9-no-mit-no-chlor.qzv \
  --m-sample-metadata-file /datasets/project_2/calf/metadata.txt
```

### Alpha Rarefaction - not working for now
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences single_rep-seqs.qza \
  --p-parttree \
  --p-mask-max-gap-frequency 0.8 \
  --p-mask-min-conservation 0.6 \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity alpha-rarefaction \
  --i-table single_table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 470000 \
  --m-metadata-file /datasets/project_2/calf/metadata.txt \
  --o-visualization alpha-rarefaction.qzv
```

