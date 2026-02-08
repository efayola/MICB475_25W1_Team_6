## Login Information

```ssh root@10.19.139.182 ```  
password: Biome1228

## Available Resources

- Basic Markdown [Syntax](https://www.markdownguide.org/basic-syntax/)
- QIIME2 [plug-in help page](https://docs.qiime2.org/2024.10/plugins/available/)
- Assembly and maturation of calf gut microbiome from neonate to post-puberty ([Paper](https://www.nature.com/articles/s41597-025-04677-7#Sec2))

## Code
**Session detached screen information**  
 ```screen -r team6_paired```

### Importing Code into QIIME2
Code adapted from QIIME2 doc under “Fastq manifest” formats ([Here](https://docs.qiime2.org/2024.10/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq)).
```
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \  
  --input-format PairedEndFastqManifestPhred33V2 \  
  --input-path /datasets/project_2/calf/manifest.txt \  
  --output-path ./paired_demux_seqs.qza
```
### Create visualization of demultiplexed samples
```
qiime demux summarize \
  --i-data /datasets/project_2/calf/demux.qza \
  --o-visualization paired_demux.qzv
```
