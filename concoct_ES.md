```mamba activate concoct-1p1```
# CONCOCT Binning Workflow

This document describes how to run **CONCOCT** on co-assembled metagenomes for groups Y and H.

---

## 📌 Requirements
- [CONCOCT](https://github.com/BinPro/CONCOCT)
- Scripts: `cut_up_fasta.py`, `concoct_coverage_table.py`, `merge_cutup_clustering.py`, `extract_fasta_bins.py`
- BAM files mapped to contigs (indexed)

---

## 🛠️ Steps

### 1. Cut contigs into chunks
```bash
cut_up_fasta.py spades_Y_coassembly/contigs.fasta   -c 10000 -o 0 --merge_last -b Y_contigs_10K.bed > Y_contigs_10K.fa

cut_up_fasta.py spades_H_coassembly/contigs.fasta   -c 10000 -o 0 --merge_last -b H_contigs_10K.bed > H_contigs_10K.fa
```

### 2. Build coverage tables
```bash
concoct_coverage_table.py Y_contigs_10K.bed bowtie/Y/*.sorted.bam > Y_coverage_table.tsv
concoct_coverage_table.py H_contigs_10K.bed bowtie/H/*.sorted.bam > H_coverage_table.tsv
```

### 3. Run CONCOCT
```bash
concoct --threads 40 --composition_file Y_contigs_10K.fa   --coverage_file Y_coverage_table.tsv -b concoct_output_Y/

concoct --threads 40 --composition_file H_contigs_10K.fa   --coverage_file H_coverage_table.tsv -b concoct_output_H/
```

### 4. Merge clusters
```bash
merge_cutup_clustering.py concoct_output_Y/clustering_gt1000.csv   > concoct_output_Y/clustering_merged_Y.csv

merge_cutup_clustering.py concoct_output_H/clustering_gt1000.csv   > concoct_output_H/clustering_merged_H.csv
```

### 5. Extract bins
```bash
extract_fasta_bins.py spades_Y_coassembly/contigs.fasta   concoct_output_Y/clustering_merged_Y.csv --output_path bins_Y/

extract_fasta_bins.py spades_H_coassembly/contigs.fasta   concoct_output_H/clustering_merged_H.csv --output_path bins_H/
```

---

## 📂 Output
- `*_contigs_10K.fa` and `*_contigs_10K.bed`
- Coverage table (`*_coverage_table.tsv`)
- Clustering results (`clustering_gt1000.csv`, merged clustering)
- FASTA bins per group in `bins_Y/` and `bins_H/`

---

# 📝 Flujo de trabajo con CONCOCT

Este documento describe cómo ejecutar **CONCOCT** sobre los metagenomas co-ensamblados para los grupos Y y H.

---

## 📌 Requerimientos
- [CONCOCT](https://github.com/BinPro/CONCOCT)
- Scripts: `cut_up_fasta.py`, `concoct_coverage_table.py`, `merge_cutup_clustering.py`, `extract_fasta_bins.py`
- Archivos BAM mapeados a contigs (con índices .bai)

---

## 🛠️ Pasos

### 1. Cortar los contigs en fragmentos
```bash
cut_up_fasta.py spades_Y_coassembly/contigs.fasta   -c 10000 -o 0 --merge_last -b Y_contigs_10K.bed > Y_contigs_10K.fa

cut_up_fasta.py spades_H_coassembly/contigs.fasta   -c 10000 -o 0 --merge_last -b H_contigs_10K.bed > H_contigs_10K.fa
```

### 2. Construir tablas de cobertura
```bash
concoct_coverage_table.py Y_contigs_10K.bed bowtie/Y/*.sorted.bam > Y_coverage_table.tsv
concoct_coverage_table.py H_contigs_10K.bed bowtie/H/*.sorted.bam > H_coverage_table.tsv
```

### 3. Ejecutar CONCOCT
```bash
concoct --threads 40 --composition_file Y_contigs_10K.fa   --coverage_file Y_coverage_table.tsv -b concoct_output_Y/

concoct --threads 40 --composition_file H_contigs_10K.fa   --coverage_file H_coverage_table.tsv -b concoct_output_H/
```

### 4. Fusionar clusters
```bash
merge_cutup_clustering.py concoct_output_Y/clustering_gt1000.csv   > concoct_output_Y/clustering_merged_Y.csv

merge_cutup_clustering.py concoct_output_H/clustering_gt1000.csv   > concoct_output_H/clustering_merged_H.csv
```

### 5. Extraer bins
```bash
extract_fasta_bins.py spades_Y_coassembly/contigs.fasta   concoct_output_Y/clustering_merged_Y.csv --output_path bins_Y/

extract_fasta_bins.py spades_H_coassembly/contigs.fasta   concoct_output_H/clustering_merged_H.csv --output_path bins_H/
```

---

## 📂 Salida
- `*_contigs_10K.fa` y `*_contigs_10K.bed`
- Tablas de cobertura (`*_coverage_table.tsv`)
- Resultados de clustering (`clustering_gt1000.csv`, clustering fusionado)
- FASTA bins por grupo en `bins_Y/` y `bins_H/`
