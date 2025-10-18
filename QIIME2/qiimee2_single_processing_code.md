### Screen session
```screen -R team6_single```

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
