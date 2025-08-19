# Metagenomic Co-assembly with SPAdes

This document describes the steps performed for the co-assembly of metagenomic data using **SPAdes** in metagenomic mode (*metaSPAdes*), starting from reads previously filtered against a reference (non-Macrohit).

## 1. Input data structure
We worked with a total of 12 samples, organized into two groups:
- **H group**: 6 pairs of R1/R2 files
- **Y group**: 6 pairs of R1/R2 files

Files follow the naming format:
```
<date>_<name>_non-macroHit_R1.fastq.gz
<date>_<name>_non-macroHit_R2.fastq.gz
```

## 2. Co-assembly script
The script automatically detects paired reads for each group and runs **metaSPAdes** to generate two co-assemblies (one per group).


```
/home/mmoreno/.local/share/mamba/envs/spades315/bin/spades.py     --only-assembler        -k      21,33,55,77,99,127      --threads       30      --memory
        200     -o      /datos1/mmorenos/fovi_data/01_bbduk_ref/spades_Y_coassembly     --pe1-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y1_S1_non-macroHit_R1.fastq.gz     --pe1-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y1_S1_non-macroHit_R2.fastq.gz     --pe2-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y2_S2_non-macroHit_R1.fastq.gz     --pe2-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y2_S2_non-macroHit_R2.fastq.gz     --pe3-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y3_S3_non-macroHit_R1.fastq.gz     --pe3-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y3_S3_non-macroHit_R2.fastq.gz     --pe4-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y4_S4_non-macroHit_R1.fastq.gz     --pe4-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y4_S4_non-macroHit_R2.fastq.gz     --pe5-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y5_S5_non-macroHit_R1.fastq.gz     --pe5-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y5_S5_non-macroHit_R2.fastq.gz     --pe6-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y6_S6_non-macroHit_R1.fastq.gz     --pe6-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y6_S6_non-macroHit_R2.fastq.gz 

/home/mmoreno/.local/share/mamba/envs/spades315/bin/spades.py     --only-assembler        -k      21,33,55,77,99,127      --threads       30      --memory
        200     -o      /datos1/mmorenos/fovi_data/01_bbduk_ref/spades_H_coassembly2    --pe1-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S7_non-macroHit_R1.fastq.gz     --pe1-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S7_non-macroHit_R2.fastq.gz     --pe2-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S8_non-macroHit_R1.fastq.gz     --pe2-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S8_non-macroHit_R2.fastq.gz     --pe3-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S10_non-macroHit_R1.fastq.gz    --pe3-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S10_non-macroHit_R2.fastq.gz    --pe4-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S9_non-macroHit_R1.fastq.gz     --pe4-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S9_non-macroHit_R2.fastq.gz     --pe5-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S11_non-macroHit_R1.fastq.gz    --pe5-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S11_non-macroHit_R2.fastq.gz    --pe6-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S12_non-macroHit_R1.fastq.gz    --pe6-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S12_non-macroHit_R2.fastq.gz```

## 3. Generated results
The script generates two output folders:
```
spades_coassemblies/spades_H_coassembly/
spades_coassemblies/spades_Y_coassembly/
```
Each folder contains:
- `contigs.fasta` → Assembled contigs
- `scaffolds.fasta` → Assembled scaffolds
- `spades.log` → Full execution log

## 4. Recommended next steps
1. Map original reads back to contigs to calculate coverage.
2. Filter contigs by length (≥ 1–1.5 kb).
3. Perform binning with **MetaBAT2**, **VAMB**, or another tool.
4. Annotate rRNA genes (optional) using **Barrnap** to check for ribosomal gene presence.
