# Genome-Resolved Metagenomics Refining Pipeline with Collections (Anvi'o)

This repository documents a full set of scripts used to **refine bins and collections** from metagenomic assemblies using **Anvi'o**.  
All scripts are kept exactly as executed, without shortening.

---

## 📂 Workflow Overview

### 1. Create a Genomes Storage Database
```bash
anvi-gen-genomes-storage --external-genomes genome_path --gene-caller prodigal -o mash_H-GENOMES.db
```

---

### 2. Generate Contigs Databases
```bash
for samples in `cat list` ; do anvi-gen-contigs-database -f  $samples -o $samples.db -n "mash_H05" ; done
anvi-gen-contigs-database -T 40 -f ./datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta -o /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db -n "mash_H05" --force-overwrite
```

---

### 3. Extract Contig Names from Genome Paths
```bash
#!/bin/bash

# Ruta del archivo de tabla con las rutas de los archivos fasta
TABLE_PATH="/media/mmorenos/disk3/006_Macro/04_bin/al_binning_mash_DASTool_bins/bins_fasta/genome_path"

# Leer la tabla línea por línea
while IFS=$'\t' read -r name contigs_db_path; do
  # Crear el archivo de salida con el nombre del genoma
  output_file="${name}.txt"

  # Usar awk para extraer los nombres de los contigs y guardarlos en el archivo de salida
  awk '/^>/{print substr($0, 2) "\t" name}' name="$name" "$contigs_db_path" | tr -cd '\11\12\15\40-\176' > "$output_file"

  echo "Archivo $output_file generado."
done < "$TABLE_PATH"
```

---

### 4. Generate Contigs Database Again (Force Overwrite)
```bash
anvi-gen-contigs-database -T 40 -f /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta -o /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db -n "mash_H05" --force-overwrite
```

---

### 5. Profile BAM Files
```bash
#!/bin/sh
anvi-profile -T 50 -i /datos1/mmorenos/fovi_data/03_bowtie/H/20250604_PamelaFernandez_DSG_H1_S7.sorted.bam  -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db
anvi-profile -T 50 -i /datos1/mmorenos/fovi_data/03_bowtie/H/20250604_PamelaFernandez_DSG_H1_S8.sorted.bam  -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db
anvi-profile -T 50 -i /datos1/mmorenos/fovi_data/03_bowtie/H/20250604_PamelaFernandez_DSG_H2_S9.sorted.bam  -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db
anvi-profile -T 50 -i /datos1/mmorenos/fovi_data/03_bowtie/H/20250604_PamelaFernandez_DSG_H2_S10.sorted.bam -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db
anvi-profile -T 50 -i /datos1/mmorenos/fovi_data/03_bowtie/H/20250604_PamelaFernandez_DSG_H3_S11.sorted.bam -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db
anvi-profile -T 50 -i /datos1/mmorenos/fovi_data/03_bowtie/H/20250604_PamelaFernandez_DSG_H3_S12.sorted.bam -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db
```

---

### 6. Merge Profiles and Summarize
```bash
anvi-merge -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db -o MERGED -P */PROFILE.db
anvi-summarize -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db -p MERGED/PROFILE.db -o MERGED/SUMMARY

anvi-merge s20250604_PamelaFernandez_DSG*/PROFILE.db -o mash_H05_samples-merged/ -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db 
```

---

### 7. Run HMMs and Taxonomy
```bash
#!/bin/sh
set -euo pipefail
DB="/datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db"

anvi-run-hmms -c "$DB" --just-do-it -T 20
anvi-run-scg-taxonomy -c "$DB" --just-do-it -T 20
```

---

### 8. Import Collections
```bash
anvi-import-collection /datos1/mmorenos/fovi_data/000_anvio/H_genomes/dastool_05/anvio_db/H1.txt -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db -p /datos1/mmorenos/fovi_data/03_bowtie/H/mash_H05_samples-merged/PROFILE.db  -C "mash_H05" --contigs-mode

anvi-import-collection /datos1/mmorenos/fovi_data/000_anvio/H_genomes/dastool_05/anvio_db/H1.txt -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db -p /datos1/mmorenos/fovi_data/03_bowtie/H/mash_H05_samples-merged/PROFILE.db  -C "mash_H05" --contigs-mode
```

---

### 9. Refine Bins Interactively
```bash
anvi-refine -c /datos1/mmorenos/fovi_data/02_spades_ass/00_coassembly_v1/spades_H_coassembly/contigs.fasta.db -p /datos1/mmorenos/fovi_data/03_bowtie/H/mash_H05_samples-merged/PROFILE.db  -C "mash_H05" -b H1

EN SERVIDOR mmorenos@hopto 

anvi-import-collection contigs_txt/H2.txt -c contigs.fasta.db -p merged_profiles/PROFILE.db -C "mash_H05" --contigs-mode  
anvi-refine -c contigs.fasta.db -p merged_profiles/PROFILE.db -C "mash_H05" -b H1
```

---

### 10. Summarize Refined Bins
```bash
$ anvi-summarize -p MERGED_PROFILE/PROFILE.db -c contigs.db -C CONCOCT -o MERGED_SUMMARY
```

---

## ⚙️ Requirements
- [Anvi'o](https://anvio.org/) (tested with v7+)
- GNU `awk`, `bash`, `coreutils`

---

## 🚀 Usage
1. Generate contigs and genomes storage databases.  
2. Extract contig names from genome paths.  
3. Map reads and create BAM profiles.  
4. Merge and summarize profiles.  
5. Run HMMs and taxonomy annotation.  
6. Import collections of bins.  
7. Refine bins interactively with Anvi'o GUI.  
8. Summarize refined collections.  

---

## 👤 Author
Scripts maintained by **Mario Moreno (mmoreno)**.  
Research focus: genome-resolved metagenomics of kelp-associated microbial communities.
