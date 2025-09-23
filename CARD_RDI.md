## The Resistance Gene Identifier (RGI) and  Comprehensive Antibiotic Resistance Database (CARD) v 4.0.0 paper from 2023 db (https://academic.oup.com/nar/article/51/D1/D690/6764414)

```
https://card.mcmaster.ca/download

mamba activate rgi
for sample in `cat list` ; do  rgi main --input_sequence $sample --output_file  $sample.rgi --local --clean -a BLAST -t protein -n 60 ; done
```

