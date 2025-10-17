## Login Information

```ssh root@10.19.139.156 ```  
password: Biome2055

## Available Resources
- QIIME2 [plug-in help page](https://docs.qiime2.org/2024.10/plugins/available/)

## Code
**Session detached screen information**  
 ```screen -S team6_paired```

### Importing Code into QIIME2
```
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \  
  --input-format PairedEndFastqManifestPhred33V2 \  
  --input-path /datasets/project_2/calf/manifest.txt \  
  --output-path ./demux_seqs.qza
```
