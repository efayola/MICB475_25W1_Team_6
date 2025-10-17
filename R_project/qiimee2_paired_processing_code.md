#### Login Code Information

```ssh root@10.19.139.156 ```  
password: Biome2055

#### Session detached screen information
 ```screen -S team6_paired```

#### Importing Code into QIIME2
```
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \  
  --input-format PairedEndFastqManifestPhred33V2 \  
  --input-path /datasets/project_2/calf/manifest.txt \  
  --output-path ./demux_seqs.qza
```
