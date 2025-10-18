### Screen session
```screen -R team6_single```

### Importing manifest
```qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/calf/manifest_fwd_only.txt \
  --output-path ./single_demux_seqs.qza```
