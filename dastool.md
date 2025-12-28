```
DAS_Tool \
  -i concoct_Y256_contigs2bin.tsv,maxbin_Y256_contigs2bin.tsv,metabat_Y256_contigs2bin.tsv \
  -l concoct,maxbin,metabat \
  -c ../../02_spades_ass/spades_Y256_coassembly/Y_S256_contigs.fasta \
  -o Y256_binning_mash \
  --write_bins \
  --write_unbinned \
  -t 40

```
