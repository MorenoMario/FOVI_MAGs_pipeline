# Genome-Resolved Metagenomics Refining Pipeline with Collections (Anvi'o)

This repository contains a set of scripts to **refine and curate bins** from metagenomic assemblies using **Anvi'o**.  
The workflow includes contigs database generation, BAM profiling, merging, taxonomy annotation, collection imports, and interactive refinement.

---

## 📂 Workflow Overview

### 1. Create a Genomes Storage Database
```bash
anvi-gen-genomes-storage --external-genomes genome_path --gene-caller prodigal -o mash_H-GENOMES.db
```
Generates a storage database from external genome FASTA paths.

---

### 2. Generate Contigs Databases
```bash
for samples in `cat list` ; do
    anvi-gen-contigs-database -f $samples -o $samples.db -n "mash_H05"
done

anvi-gen-contigs-database -T 40     -f ./datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta     -o /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db     -n "mash_H05" --force-overwrite
```
Creates contigs databases for bins and for the co-assembly.

---

### 3. Extract Contig Names from Genome Paths
```bash
#!/bin/bash
TABLE_PATH="/media/mmorenos/disk3/006_Macro/04_bin/al_binning_mash_DASTool_bins/bins_fasta/genome_path"

while IFS=$'\t' read -r name contigs_db_path; do
  output_file="${name}.txt"
  awk '/^>/{print substr($0, 2) "\t" name}' name="$name" "$contigs_db_path" | tr -cd '\11\12\15\40-\176' > "$output_file"
  echo "File $output_file generated."
done < "$TABLE_PATH"
```
Parses genome names and extracts contig IDs into text files for bin assignment.

---

### 4. Profile BAM Files
```bash
#!/bin/sh
anvi-profile -T 50 -i sample1.sorted.bam -c contigs.fasta.db
anvi-profile -T 50 -i sample2.sorted.bam -c contigs.fasta.db
...
```
Generates coverage profiles from BAM files aligned to the co-assembly.

---

### 5. Merge Profiles and Summarize
```bash
anvi-merge */PROFILE.db   -o mash_H05_samples-merged/   -c contigs.fasta.db

anvi-summarize -c contigs.fasta.db   -p mash_H05_samples-merged/PROFILE.db   -o mash_H05_samples-merged/SUMMARY
```
Combines multiple profiles and produces summaries.

---

### 6. Run HMMs and Taxonomy
```bash
#!/bin/sh
set -euo pipefail
DB="contigs.fasta.db"

anvi-run-hmms -c "$DB" --just-do-it -T 20
anvi-run-scg-taxonomy -c "$DB" --just-do-it -T 20
```
Annotates single-copy genes and assigns taxonomy.

---

### 7. Import Collections
```bash
anvi-import-collection /path/to/H1.txt   -c contigs.fasta.db   -p mash_H05_samples-merged/PROFILE.db   -C "mash_H05" --contigs-mode
```
Imports pre-defined collections of contigs (e.g., bin assignments).

---

### 8. Refine Bins Interactively
```bash
anvi-refine -c contigs.fasta.db   -p mash_H05_samples-merged/PROFILE.db   -C "mash_H05" -b H1
```
Opens the interactive interface to refine bins manually.

---

### 9. Summarize Refined Bins
```bash
anvi-summarize -p MERGED_PROFILE/PROFILE.db   -c contigs.db   -C CONCOCT   -o MERGED_SUMMARY
```
Summarizes refined bins after collections and refinement.

---

## ⚙️ Requirements
- [Anvi'o](https://anvio.org/) (tested with v7+)
- GNU `awk`, `bash`, `coreutils`

---

## 🚀 Usage
1. Generate contigs and genomes storage databases.  
2. Map reads and create BAM profiles.  
3. Merge and summarize profiles.  
4. Import collections of bins.  
5. Refine bins interactively with Anvi'o GUI.  
6. Summarize refined collections.  

---

## 📑 Notes
- Use `--force-overwrite` with caution.  
- Refinement is done interactively after importing collections.  
- Scripts assume a Linux HPC/Server environment.  

---

## 👤 Author
Scripts maintained by **Mario Moreno (mmoreno)**.  
Research focus: genome-resolved metagenomics of kelp-associated microbial communities.
