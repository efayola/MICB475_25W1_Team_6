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
