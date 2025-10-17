## Login Information

```ssh root@10.19.139.156 ```  
password: Biome2055

## Available Resources

- Basic Markdown [Syntax](https://www.markdownguide.org/basic-syntax/)
- QIIME2 [plug-in help page](https://docs.qiime2.org/2024.10/plugins/available/)
- Assembly and maturation of calf gut microbiome from neonate to post-puberty ([Paper](https://www.nature.com/articles/s41597-025-04677-7#Sec2))

## Code
**Session detached screen information**  
 ```screen -S team6_paired```

### Importing Code into QIIME2
Code adapted from QIIME2 doc under “Fastq manifest” formats ([Here](https://docs.qiime2.org/2024.10/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq)).
```
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \  
  --input-format PairedEndFastqManifestPhred33V2 \  
  --input-path /datasets/project_2/calf/manifest.txt \  
  --output-path ./demux_seqs.qza
```
### Demultiplexing using `cutadapt`  

My reasoning for using `cutadapt` plugins instead of `demux` plugins for analysis:

- in the paper they used `cutadapt`
- `demux` is used for data processed from EMT protocol ([Here](https://docs.qiime2.org/2024.10/plugins/available/demux/emp-paired/)), but we used manifest.
- Code adapted from the `cutadapt` ([Here](https://docs.qiime2.org/2024.10/plugins/available/cutadapt/demux-paired/)).

```
qiime cutadapt demux-paired \
  --i-seqs demux_seqs.qza \
  --m-forward-barcodes-file /datasets/project_2/calf/metadata.txt \
  --m-forward-barcodes-column sample-id \
  --o-per-sample-sequences demux_paired.qza \
  --o-untrimmed-sequences untrimmed_seqs_paired.qza
```
