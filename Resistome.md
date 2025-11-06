## The Resistance Gene Identifier (RGI) and  Comprehensive Antibiotic Resistance Database (CARD) v 4.0.0 paper from 2023 db (https://academic.oup.com/nar/article/51/D1/D690/6764414)

## https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst
```

mamba activate rgi
rgi auto_load -i localDB/card.json
for sample in `cat list` ; do  rgi main --input_sequence $sample --output_file  $sample.rgi --local --clean -a DIAMOND -t protein -n 60 --low_quality --include_nudge ; done
/media/mmorenos/disk3/DSG123-selected/001_anvio_metag/Y_genomes
```
```

conda activate abricate_env
abricate --db card ilq_contigs.fasta.db.contigs.fasta >  ilq_card.tsv
```
