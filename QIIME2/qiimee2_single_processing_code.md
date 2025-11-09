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

### Export File from QIIME2

```
#export table.qza file
qiime tools export \
	--input-path ../team6_calf/single_table.qza \
	--output-path table_export

#To export a biom file into a text file:
biom convert \
-i feature-table.biom \
--to-tsv \
-o feature-table.txt

#To export single_taxonomy.qza file
qiime tools export \
	--input-path ../team6_calf/single_taxonomy.qza \
	--output-path taxonomy_export

#To export rooted-tree.qza file
qiime tools export \
	--input-path ../team6_calf/rooted-tree.qza \
	--output-path tree_export

scp -r root@10.19.139.156:/data/rstudio_export_team6_single .
```

### Code Ran to push feature table (from grok)
1. **Initialize Git LFS in the repo** (if not already done; this is safe to run multiple times): <br>
	`git lfs install`
2. **Track the large file with LFS**: This adds an entry to `.gitattributes` for files matching this pattern. If you have other large .txt files, you can broaden it to `"*.txt"`. <br>
	`git lfs track "feature-table.txt"`
3. **Remove the file from Git's index (without deleting it locally)**: This unstages it so you can re-add it under LFS. <br>
	`git rm --cached R_project/single_table_export/feature-table.txt`
4. **Re-add the file** (now it will be handled by LFS): <br>
	`git add R_project/single_table_export/feature-table.txt`
5. **Commit the changes**: If this is amending an existing commit with the large file, use `--amend`. Otherwise, just commit normally. Also commit the `.gitattributes` file. <br>
	`git add .gitattributes` <br>
	`git commit -m "Track large file with Git LFS"`
6. **Push to GitHub**: Force push if you amended history (use with caution if others have pulled your branch). <br>
	`git push origin main`
7. **Migrate the history** to convert the file to LFS in all commits: <br>
	`git lfs migrate import --everything --include="R_project/single_table_export/feature-table.txt"`
8. **Push to GitHub**: Since history was rewritten, use a safe force push: <br>
	`git push origin main --force-with-lease`


### Exact Code & return for terminal
```
Louiss-MacBook-Air:R_project zhulouis$ git lfs install
Updated Git hooks.
Git LFS initialized.
Louiss-MacBook-Air:R_project zhulouis$ git lfs track "feature-table.txt"
Tracking "feature-table.txt"
Louiss-MacBook-Air:R_project zhulouis$ git rm --cached R_project/single_table_export/feature-table.txt
fatal: pathspec 'R_project/single_table_export/feature-table.txt' did not match any files
Louiss-MacBook-Air:R_project zhulouis$ git add R_project/single_table_export/feature-table.txt
warning: could not open directory 'R_project/R_project/single_table_export/': No such file or directory
fatal: pathspec 'R_project/single_table_export/feature-table.txt' did not match any files
Louiss-MacBook-Air:R_project zhulouis$ git add single_table_export/feature-table.txt
Louiss-MacBook-Air:R_project zhulouis$ git add .gitattributes
Louiss-MacBook-Air:R_project zhulouis$ git commit -m "Track large file with Git LFS"
[main f7ac4bf] Track large file with Git LFS
 2 files changed, 4 insertions(+), 61958 deletions(-)
 create mode 100644 R_project/.gitattributes
Louiss-MacBook-Air:R_project zhulouis$ git push origin main
Uploading LFS objects: 100% (1/1), 122 MB | 2.8 MB/s, done.                                                                    
Enumerating objects: 18, done.
Counting objects: 100% (18/18), done.
Delta compression using up to 8 threads
Compressing objects: 100% (14/14), done.
Writing objects: 100% (16/16), 8.32 MiB | 1.08 MiB/s, done.
Total 16 (delta 4), reused 0 (delta 0), pack-reused 0 (from 0)
remote: Resolving deltas: 100% (4/4), completed with 1 local object.
remote: error: Trace: 8c843388fff442627830794acf927ea4a6dc148720e135758ce24ff1956276b9
remote: error: See https://gh.io/lfs for more information.
remote: error: File R_project/single_table_export/feature-table.txt is 116.37 MB; this exceeds GitHub's file size limit of 100.00 MB
remote: error: GH001: Large files detected. You may want to try Git Large File Storage - https://git-lfs.github.com.
To https://github.com/efayola/MICB475_25W1_Team_6.git
 ! [remote rejected] main -> main (pre-receive hook declined)
error: failed to push some refs to 'https://github.com/efayola/MICB475_25W1_Team_6.git'

Louiss-MacBook-Air:R_project zhulouis$ git lfs install
Updated Git hooks.
Git LFS initialized.
Louiss-MacBook-Air:R_project zhulouis$ git lfs track "R_project/single_table_export/feature-table.txt"
Tracking "R_project/single_table_export/feature-table.txt"
Louiss-MacBook-Air:R_project zhulouis$ git add .gitattributes
Louiss-MacBook-Air:R_project zhulouis$ git commit -m "Add LFS tracking for large file"  # Or amend if needed
[main f65aac5] Add LFS tracking for large file
 2 files changed, 5 insertions(+), 61958 deletions(-)
 create mode 100644 R_project/.gitattributes
Louiss-MacBook-Air:R_project zhulouis$ git lfs migrate import --everything --include="R_project/single_table_export/feature-table.txt"
override changes in your working copy?  All uncommitted changes will be lost! [y/N] y
changes in your working copy will be overridden ...
Sorting commits: ..., done.                                                                                                    
Rewriting commits: 100% (147/147), done.                                                                                       
  main  f65aac5a2c538ce501816c4a80e1065210701315 -> 860ac25422e2d176f537e40ccf6eb189a29d0c28
Updating refs: ..., done.                                                                                                      
Checkout: ..., done.                                                                                                           
Louiss-MacBook-Air:R_project zhulouis$ git push origin main --force-with-lease
Uploading LFS objects: 100% (1/1), 122 MB | 0 B/s, done.                                                                       
Enumerating objects: 559, done.
Counting objects: 100% (559/559), done.
Delta compression using up to 8 threads
Compressing objects: 100% (367/367), done.
Writing objects: 100% (559/559), 40.51 MiB | 1.47 MiB/s, done.
Total 559 (delta 294), reused 256 (delta 156), pack-reused 0 (from 0)
remote: Resolving deltas: 100% (294/294), done.
To https://github.com/efayola/MICB475_25W1_Team_6.git
 + ec4ae31...860ac25 main -> main (forced update)
Louiss-MacBook-Air:R_project zhulouis$ 
```
