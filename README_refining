# Genome-Resolved Metagenomics Refining Pipeline (Anvi'o)

This repository contains a set of scripts used to **refine and curate bins** obtained from metagenomic assemblies using **Anvi'o**.  
The workflow covers contigs database generation, profile creation from mapped BAM files, merging, summarization, and taxonomy annotation.

---

## 📂 Workflow Overview

### 1. Generate Contigs Databases from Bins or Samples
```bash
for samples in `cat list` ; do
    anvi-gen-contigs-database -f $samples -o $samples.db -n "mash_H05"
done
```
Creates Anvi'o contigs databases for each FASTA file listed in `list` (typically bins from automated binning).

---

### 2. Generate Contigs Database from Co-Assembly
```bash
anvi-gen-contigs-database -T 40     -f ./datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta     -o /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db     -n "mash_H05" --force-overwrite
```
Builds a contigs database for the co-assembled dataset.

---

### 3. Extract Contig Names from Bin Paths
```bash
#!/bin/bash
TABLE_PATH="/media/mmorenos/disk3/006_Macro/04_bin/al_binning_mash_DASTool_bins/bins_fasta/genome_path"

while IFS=$'\t' read -r name contigs_db_path; do
  output_file="${name}.txt"
  awk '/^>/{print substr($0, 2) "\t" name}' name="$name" "$contigs_db_path" | tr -cd '\11\12\15\40-\176' > "$output_file"
  echo "File $output_file generated."
done < "$TABLE_PATH"
```
Parses a table of genome names and FASTA paths to extract contig IDs, used for refining bin content.

---

### 4. Profile BAM Files Against the Co-Assembly
```bash
#!/bin/sh
anvi-profile -T 50 -i sample1.sorted.bam -c contigs.fasta.db
anvi-profile -T 50 -i sample2.sorted.bam -c contigs.fasta.db
...
```
Generates Anvi'o profiles from BAM files (mapping reads back to the co-assembly).  
These profiles are essential for **bin refinement** (coverage patterns, tetranucleotide frequency).

---

### 5. Merge Profiles and Summarize
```bash
anvi-merge */PROFILE.db   -o mash_H05_samples-merged/   -c contigs.fasta.db

anvi-summarize -c contigs.fasta.db   -p mash_H05_samples-merged/PROFILE.db   -o mash_H05_samples-merged/SUMMARY
```
Combines multiple profiles and produces summary reports, facilitating **interactive refining** of bins.

---

### 6. Run HMMs and Taxonomy Annotation
```bash
#!/bin/sh
set -euo pipefail
DB="contigs.fasta.db"

anvi-run-hmms -c "$DB" --just-do-it -T 20
anvi-run-scg-taxonomy -c "$DB" --just-do-it -T 20
```
Executes HMM search and single-copy gene taxonomy assignment for evaluating completeness and contamination of bins.

---

## ⚙️ Requirements
- [Anvi'o](https://anvio.org/) (tested with v7+)
- GNU `awk`, `bash`, `coreutils`

---

## 🚀 Usage
1. Generate contigs databases for bins or co-assembly.
2. Map reads → create BAM files.
3. Profile BAMs with `anvi-profile`.
4. Merge profiles → summarize results.
5. Run HMMs and taxonomy → use for **refinement in Anvi'o interactive interface**.

---

## 📑 Notes
- Databases will be overwritten if `--force-overwrite` is used.
- Refinement is done interactively after generating the profiles and summaries.
- Use `set -euo pipefail` for safer bash execution.

---

## 👤 Author
Scripts maintained by **Mario Moreno (mmoreno)**.  
Research focus: genome-resolved metagenomics of kelp-associated microbial communities.
