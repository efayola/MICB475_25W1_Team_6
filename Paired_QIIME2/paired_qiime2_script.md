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
```
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \  
  --input-format PairedEndFastqManifestPhred33V2 \  
  --input-path /datasets/project_2/calf/manifest.txt \  
  --output-path ./new_paired_demux_seqs.qza
```
### Create visualization of demultiplexed samples
```
qiime demux summarize \
  --i-data paired_demux_seqs.qza \
  --o-visualization paired_demux.qzv
```
```
qiime demux summarize \
  --i-data new_paired_demux_seqs.qza \
  --o-visualization new_paired_demux.qzv
```
**Error Log Fixed**
```
(qiime2-amplicon-2025.4) root@stu-13829:/data/team6_paired# cat /tmp/qiime2-q2cli-err-70n694_s.log
Traceback (most recent call last):
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/site-packages/q2cli/commands.py", line 529, in __call__
    results = self._execute_action(
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/site-packages/q2cli/commands.py", line 601, in _execute_action
    results = action(**arguments)
  File "<decorator-gen-888>", line 2, in summarize
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/site-packages/qiime2/sdk/action.py", line 221, in bound_callable
    outputs = self._callable_executor_(
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/site-packages/qiime2/sdk/action.py", line 408, in _callable_executor_
    ret_val = self._callable(output_dir=temp_dir, **view_args)
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/site-packages/q2_demux/_summarize/_visualizer.py", line 131, in summarize
    for seq in read_fastq_seqs(filename):
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/site-packages/q2_demux/_util.py", line 18, in read_fastq_seqs
    for seq_header, seq, qual_header, qual in itertools.zip_longest(*[fh] * 4):
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/gzip.py", line 314, in read1
    return self._buffer.read1(size)
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/_compression.py", line 68, in readinto
    data = self.read(len(byte_view))
  File "/opt/conda/envs/qiime2-amplicon-2025.4/lib/python3.10/gzip.py", line 507, in read
    raise EOFError("Compressed file ended before the "
EOFError: Compressed file ended before the end-of-stream marker was reached
```
### Download
```scp root@10.19.139.182:/data/team6_paired/paired_demux.qzv .```
### Denoising with DADA2
```
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs paired_demux_seqs.qza \
    --p-trunc-len-f 228 \
    --p-trunc-len-r 213 \
    --o-representative-sequences paired_representative-sequences.qza \
    --o-table paired_table.qza \
    --o-denoising-stats paired_denoising-stats.qza
```
**Checking Merging Quality**
```
qiime metadata tabulate \
  --m-input-file paired_denoising-stats.qza \
  --o-visualization paired_denoising-stats.qzv

scp root@10.19.139.182:/data/team6_paired/paired_denoising-stats.qzv .
```
**Forward 260 Reverse 200**
- Forward Reads: Above Q25 until 236 bp. Above Q20 until 266 bp.
- Reverse Reads: Above Q25 until 195 bp. Above Q20 until 214 bp.
```
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ../paired_demux_seqs.qza \
    --p-trunc-len-f 260 \
    --p-trunc-len-r 200 \
    --o-representative-sequences paired_representative-sequences_f260_r200.qza \
    --o-table paired_table_f260_r200.qza \
    --o-denoising-stats paired_denoising-stats_f260_r200.qza
```

```
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ../paired_demux_seqs.qza \
    --p-trunc-len-f 266 \
    --p-trunc-len-r 214 \
    --o-representative-sequences paired_representative-sequences_f266_r214.qza \
    --o-table paired_table_f266_r214.qza \
    --o-denoising-stats paired_denoising-stats_f266_r214.qza
```

```
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ../paired_demux_seqs.qza \
    --p-trunc-len-f 255 \
    --p-trunc-len-r 195 \
    --o-representative-sequences paired_representative-sequences_f255_r195.qza \
    --o-table paired_table_f255_r195.qza \
    --o-denoising-stats paired_denoising-stats_f255_r195.qza
```

### Training Classifiers
```
qiime feature-classifier extract-reads \
  --i-sequences /datasets/classifiers/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-min-length 400 \
  --p-max-length 500 \
  --o-reads ref-seqs-trimmed400-500.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed400-500.qza \
  --i-reference-taxonomy /datasets/classifiers/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier400-500.qza
```
