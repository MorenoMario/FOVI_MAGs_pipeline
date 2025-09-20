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
anvi-summarize -p merged_profiles/PROFILE.db -c contigs.fasta.db -C "mash_H05" -o H1
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

# Recommended Workflow for Bin Refinement (Anvi'o)

This guide provides a **quick and effective workflow** for refining bins in Anvi'o, using sequence composition, differential coverage, and SCG/taxonomy checks.

---

## 🔹 Step-by-Step Workflow

### 1. Start with **Seq. Composition + Diff. Coverage**
- Use the combined view of tetranucleotide frequencies (TNFs/GC) and coverage.  
- Visually mark clusters that separate clearly in this view.  
- Check **SCGs (single-copy genes)** and **taxonomy** for those clusters.

---

### 2. Switch to **Sequence Composition**
- If the same clusters remain separated → the difference is real in TNFs/GC → **likely contamination**.  
- If they merge again → the separation was caused by coverage differences → **possible co-occurring populations**.  
  - Decision depends on your biological question:  
    - **Strictly monoclonic MAG?** → separate them.  
    - **Broad functional unit?** → keep them together.

---

### 3. Test **Differential Coverage**
- Confirms whether robust coverage patterns support the separation.  
- Useful to decide between **true strain-level populations** vs. **noise**.

---

### 4. Apply Soft Filters
- Set a **minimum contig length** (2–2.5 kb) temporarily to clarify patterns.  
- Exclude contigs with **ultra-low coverage** (e.g., <1–2×) if they clutter the clustering.

---

### 5. Practical Rule of Decision
- **Remove a cluster** if:
  - It reduces **SCG duplication**, AND
  - It improves the **smooth coverage profile** of the bin.  
- **Keep the cluster** if:
  - Removing it **breaks coverage continuity**, OR
  - It drastically lowers completeness without solving contamination → may represent a **real subpopulation**.

---

## ✅ Summary
- Start with **Seq. Composition + Diff. Coverage** as default.  
- Cross-check with **Sequence Composition** and **Differential Coverage** views.  
- Use **SCG duplication**, **taxonomy inconsistency**, and **coverage smoothness** as the main criteria.  
- Apply filters to clarify the signal, not as hard exclusions.  
- Always align the decision with your **biological goal** (strict MAG vs. population bin).

---
