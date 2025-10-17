## Login Information

```ssh root@10.19.139.156 ```  
password: Biome2055

## Available Resources

- Basic Markdown [Syntax](https://www.markdownguide.org/basic-syntax/)
- QIIME2 [plug-in help page](https://docs.qiime2.org/2024.10/plugins/available/)
- Assembly and maturation of calf gut microbiome from neonate to post-puberty ([Paper](https://www.nature.com/articles/s41597-025-04677-7#Sec2))

## Code
**Session detached screen information**  
 ```screen -R team6_paired```

### Importing Code into QIIME2
Code adapted from QIIME2 doc under “Fastq manifest” formats ([Here](https://docs.qiime2.org/2024.10/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq)).
```
qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \  
  --input-format PairedEndFastqManifestPhred33V2 \  
  --input-path /datasets/project_2/calf/manifest.txt \  
  --output-path ./demux_seqs.qza
```
### Create visualization of demultiplexed samples
```
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux.qzv
```
#### error file
```/tmp/qiime2-q2cli-err-jngv22ej.log ```
'''
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
'''
