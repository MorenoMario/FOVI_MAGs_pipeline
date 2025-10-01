# Anvi’o Pipeline Script Log

This document records the commands used during the processing of metagenomes with Anvi’o and external tools.  
Each section includes the exact command and a brief explanation.

---

## 1. Generate the contigs database

```bash
anvi-gen-contigs-database   -f ilq_contigs.fasta.db.contigs.fasta   -n "fovi_il_chon"   -T 50
```

**Description:**  
Creates the contigs database (`contigs.db`) from the assembly file `ilq_contigs.fasta.db.contigs.fasta`.  

- `-f`: FASTA file containing the contigs.  
- `-n`: name assigned to the database.  
- `-T`: number of threads used (here, 50).  

---

## 2. Generate the contigs database (variant with different name)

```bash
anvi-gen-contigs-database   -f ilq_contigs.fasta.db.contigs.fasta   -n "fovi_IlChon"   -T 50
```

**Description:**  
Same step as above, but the database name is written as `"fovi_IlChon"` (capital **I** and **C**).  
Useful if you want to keep multiple versions of the contigs database with slightly different labels.  

---

## 3. Functional annotation with NCBI COGs (COG24, using DIAMOND)

```bash
nohup anvi-run-ncbi-cogs   -c chon_contigs_fovi.fasta.db   --cog-data-dir /home/mmorenos/github/anvio/anvio/data/misc/COG/   --search-with diamond   -T 20   --cog-version COG24 &
```

```bash
nohup anvi-run-ncbi-cogs   -c ilq_contigs_fovi.fasta.db   --cog-data-dir /home/mmorenos/github/anvio/anvio/data/misc/COG/   --search-with diamond   -T 20   --cog-version COG24 &
```

**Description:**  
Runs functional annotation of gene calls against the **NCBI COG (version 24)** database using **DIAMOND** for faster similarity searches.  

- `nohup ... &`: runs the job in the background, appending logs to `nohup.out`.  
- `-c`: contigs database to annotate.  
- `--cog-data-dir`: path to the COG database directory.  
- `--search-with diamond`: uses DIAMOND instead of BLAST for speed.  
- `-T`: number of threads (20).  
- `--cog-version`: specifies the COG version (`COG24`).  

---

## 4. Functional annotation with CAZymes

```bash
nohup anvi-run-cazymes   -c ilq_contigs_fovi.fasta.db   --cazyme-data-dir /home/mmorenos/github/anvio/anvio/data/misc/CAZyme/   -T 30 &
```

```bash
nohup anvi-run-cazymes   -c chon_contigs_fovi.fasta.db   --cazyme-data-dir /home/mmorenos/github/anvio/anvio/data/misc/CAZyme/   -T 30 &
```

**Description:**  
Runs functional annotation of carbohydrate-active enzymes (CAZymes) using the CAZy database.  

- `nohup ... &`: runs the job in the background, appending logs to `nohup.out`.  
- `-c`: contigs database to annotate.  
- `--cazyme-data-dir`: path to the CAZyme database directory.  
- `-T`: number of threads (30).  

---

## 5. Taxonomic annotation with Kaiju

First, export the nucleotide sequences of predicted genes from each contigs database:

```bash
anvi-get-sequences-for-gene-calls   -c ilq_contigs_fovi.fasta.db   -o ilq_contigs_fovi.gene_calls.fna

anvi-get-sequences-for-gene-calls   -c chon_contigs_fovi.fasta.db   -o chon_contigs_fovi.gene_calls.fna
```

Then, run **Kaiju** against the RefSeq NR protein database:

```bash
nohup kaiju   -t /media/mmorenos/disk1/database/kaijudb_all/refseq_nr_prot_06_17_23/nodes.dmp   -f /media/mmorenos/disk1/database/kaijudb_all/refseq_nr_prot_06_17_23/kaiju_db_refseq_nr.fmi   -i ilq_contigs_fovi.gene_calls.fna   -o ilq_contigs_fovi.gene_nr.out   -z 30   -v &> kaiju.ilq.txt &
```

```bash
nohup kaiju   -t /media/mmorenos/disk1/database/kaijudb_all/refseq_nr_prot_06_17_23/nodes.dmp   -f /media/mmorenos/disk1/database/kaijudb_all/refseq_nr_prot_06_17_23/kaiju_db_refseq_nr.fmi   -i chon_contigs_fovi.gene_calls.fna   -o chon_contigs_fovi.gene_nr.out   -z 30   -v &> kaiju.chon.txt &
```

**Description:**  
Performs taxonomic annotation of predicted gene sequences using **Kaiju** against the RefSeq NR protein database.  

- `anvi-get-sequences-for-gene-calls`: exports nucleotide sequences of predicted genes from the contigs databases.  
- `-t`: taxonomy nodes file (`nodes.dmp`).  
- `-f`: Kaiju database index (`kaiju_db_refseq_nr.fmi`).  
- `-i`: input nucleotide sequences of genes.  
- `-o`: output file with Kaiju assignments.  
- `-z`: number of threads (30).  
- `-v`: verbose mode.  
- `nohup ... &`: runs in the background, logs appended to output files (`kaiju.ilq.txt`, `kaiju.chon.txt`).  

---
