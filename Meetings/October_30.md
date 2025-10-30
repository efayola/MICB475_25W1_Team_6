## Agenda
* did taxonomy-based filtering, but should we do metadata-based filtering as well?
* check frequency per sample max freuquency when generating rarefaction curve
* check sampling depth: 49216 - 49303
  
## Notes
qiime phylogeny align-to-tree-mafft-fasttree \ --i-sequences single_rep-seqs.qza \ --o-alignment aligned-rep-seqs.qza \ --o-masked-alignment masked-aligned-rep-seqs.qza \ --o-tree unrooted-tree.qza \ --o-rooted-tree rooted-tree.qza 


Plugin error from phylogeny: Command '['mafft', '--preservecase', '--inputorder', '--thread', '1', '/tmp/qiime2/root/data/46e08cd8-4631-4518-9907-208315a3f264/data/dna-sequences.fasta']' returned non-zero exit status 1. Debug info has been saved to /tmp/qiime2-q2cli-err-_mamjrnl.log (qiime2-amplicon-2025.4) root@stu-31774:/data/team6_calf#
