

```
workinkg with genomes from dastool threshold 0.0

for samples in `cat list` ; do anvi-gen-contigs-database -f  $samples -o $samples.db -n "fovi_H_samples" -T 50 --force-overwrite ; done

anvi-gen-contigs-database -T 40 -f ../00_contigs_bow/doc_contigs.fasta -o doc_contigs.fasta.db -n "mash_north_doc" --force-overwrite
```

