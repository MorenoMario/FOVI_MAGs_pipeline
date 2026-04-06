# FOVI MAGs Pipeline

> **A complete, reproducible metagenome-assembled genome (MAG) reconstruction pipeline for environmental metagenomics.**

This repository documents the full bioinformatics workflow used to assemble, bin, refine, and annotate metagenome-assembled genomes (MAGs) from paired-end Illumina metagenomic reads. The pipeline was developed for the **FOVI** project and covers everything from raw read quality control to functional annotation and resistome analysis.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Pipeline Architecture](#2-pipeline-architecture)
3. [Requirements](#3-requirements)
4. [Directory Structure](#4-directory-structure)
5. [Stage 1 — Read Quality Control (Skewer)](#5-stage-1--read-quality-control-skewer)
6. [Stage 2 — Host Read Removal (BBDuk)](#6-stage-2--host-read-removal-bbduk)
7. [Stage 3 — Co-Assembly (metaSPAdes)](#7-stage-3--co-assembly-metaspades)
8. [Stage 4 — Read Mapping (Bowtie2)](#8-stage-4--read-mapping-bowtie2)
9. [Stage 5 — Metagenomic Binning](#9-stage-5--metagenomic-binning)
   - [5a. MetaBAT2](#5a-metabat2)
   - [5b. MaxBin2](#5b-maxbin2)
   - [5c. CONCOCT](#5c-concoct)
   - [5d. DAS_Tool (bin refinement)](#5d-das_tool-bin-refinement)
10. [Stage 6 — MAG Refinement & Quality (Anvi'o + CheckM2)](#10-stage-6--mag-refinement--quality-anvio--checkm2)
11. [Stage 7 — Functional Annotation (COG14 + Pfam)](#11-stage-7--functional-annotation-cog14--pfam)
12. [Stage 8 — Resistome Analysis](#12-stage-8--resistome-analysis)
13. [Output Summary](#13-output-summary)
14. [Citation](#14-citation)

---

## 1. Overview

This pipeline reconstructs MAGs from co-assembled Illumina short-read metagenomes. It is structured as a series of modular stages, each of which can be run independently once the required inputs from the previous stage are available.

**Sample groups:** The project contains two sets of metagenomes referred to as **H samples** (6 paired-end libraries) and **Y samples** (6 paired-end libraries). Both groups are co-assembled independently and then processed through the full pipeline.

**Key outputs at a glance:**

| Output | Tool | Format |
|--------|------|--------|
| Quality-filtered reads | Skewer / BBDuk | `.fastq.gz` |
| Co-assembled contigs | metaSPAdes | `contigs.fasta` |
| Coverage profiles | Bowtie2 + jgi | `.bam`, depth matrix |
| Initial bins | MetaBAT2, MaxBin2, CONCOCT | `.fa` per bin |
| Refined bins | DAS_Tool | `.fa` per bin |
| Curated MAGs | Anvi'o | `.fa` per MAG |
| MAG quality | CheckM2 | `.tsv` report |
| Functional annotation | DIAMOND + HMMER | Anvi'o function tables |
| Resistance genes | RGI, ABRicate, ARGS-OAP | `.txt` / `.tsv` |

---

## 2. Pipeline Architecture

```
Raw reads (Illumina paired-end)
        │
        ▼
┌───────────────────────┐
│ 1. Quality Control    │  Skewer (adapter trimming + quality filtering)
└───────────┬───────────┘
            │
            ▼
┌───────────────────────┐
│ 2. Host Read Removal  │  BBDuk (k-mer based read filtering)
└───────────┬───────────┘
            │
            ▼
┌───────────────────────┐
│ 3. Co-Assembly        │  metaSPAdes (k=21,33,55,77,99,127)
└───────────┬───────────┘
            │
       ┌────┴────┐
       │         │
       ▼         ▼
┌──────────┐  ┌──────────────────────────────────┐
│ contigs  │  │ 4. Read Mapping (Bowtie2)         │
│ .fasta   │  │    → depth/coverage profiles      │
└──────────┘  └────────────────┬─────────────────┘
                               │
              ┌────────────────┼────────────────┐
              │                │                │
              ▼                ▼                ▼
        ┌──────────┐   ┌──────────┐   ┌──────────────┐
        │ MetaBAT2 │   │ MaxBin2  │   │   CONCOCT    │
        └────┬─────┘   └────┬─────┘   └──────┬───────┘
             │              │                 │
             └──────────────┼─────────────────┘
                            │
                            ▼
                   ┌─────────────────┐
                   │  5. DAS_Tool    │  Bin dereplication & scoring
                   └────────┬────────┘
                            │
                            ▼
                   ┌─────────────────┐
                   │  6. Anvi'o      │  Interactive refinement & curation
                   │  + CheckM2      │  Completeness / contamination
                   └────────┬────────┘
                            │
               ┌────────────┴───────────┐
               │                        │
               ▼                        ▼
   ┌───────────────────┐    ┌────────────────────────┐
   │ 7. Annotation     │    │ 8. Resistome            │
   │  COG14 + Pfam v32 │    │  RGI / ABRicate /       │
   │  (DIAMOND, HMMER) │    │  ARGS-OAP               │
   └───────────────────┘    └────────────────────────┘
```

---

## 3. Requirements

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| [Skewer](https://github.com/relipmoc/skewer) | ≥0.2.2 | Adapter trimming |
| [BBTools (BBDuk)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | ≥38 | Host read filtering |
| [SPAdes / metaSPAdes](https://github.com/ablab/spades) | ≥3.15 | Metagenomic assembly |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) | ≥2.4 | Read mapping |
| [SAMtools](http://www.htslib.org/) | ≥1.12 | BAM processing |
| [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | ≥2.15 | Binning |
| [MaxBin2](https://sourceforge.net/projects/maxbin2/) | ≥2.2.7 | Binning |
| [CONCOCT](https://github.com/BinPro/CONCOCT) | 1.1.0 | Binning |
| [DAS_Tool](https://github.com/cmks/DAS_Tool) | ≥1.1.4 | Bin dereplication |
| [Anvi'o](https://anvio.org/) | ≥7.0 | Refinement & visualization |
| [CheckM2](https://github.com/chklovski/CheckM2) | ≥1.0 | MAG quality |
| [DIAMOND](https://github.com/bbuchfink/diamond) | ≥2.0 | Protein search |
| [HMMER](http://hmmer.org/) | ≥3.3 | Profile HMM search |
| [RGI](https://github.com/arpcard/rgi) | ≥6.0 | Resistance genes |
| [ABRicate](https://github.com/tseemann/abricate) | ≥1.0 | Resistance screening |
| [ARGS-OAP](https://github.com/biofuture/Ublastx_stageone) | ≥3.0 | ARG profiling |

### Databases

| Database | Version | Used by |
|----------|---------|---------|
| COG (eggNOG) | 2014 (COG14) | DIAMOND annotation |
| Pfam-A | v32 | HMMER annotation |
| CARD | v4.0.0 | RGI |
| argannot, ncbi, resfinder, plasmidfinder | Latest | ABRicate |

### Compute Resources

| Stage | Threads | Memory |
|-------|---------|--------|
| Skewer / BBDuk | 16 (default) | ~16 GB |
| metaSPAdes | 30 | **200 GB** |
| Bowtie2 | 30 | ~50 GB |
| MetaBAT2 / MaxBin2 | 30–60 | ~50 GB |
| CONCOCT | 40 | ~50 GB |
| Anvi'o profiling | 50 | ~100 GB |
| DIAMOND / HMMER | 30–50 | ~30 GB |
| RGI | 60 | ~50 GB |

> **Note:** metaSPAdes co-assembly requires at least **200 GB RAM**. All long-running jobs should be submitted via `nohup` or a job scheduler (SLURM/SGE).

---

## 4. Directory Structure

A recommended project layout:

```
FOVI/
├── 00_raw/                    # Raw paired-end FASTQ files
│   ├── H_samples/
│   │   ├── H1_R1.fastq.gz
│   │   └── H1_R2.fastq.gz
│   └── Y_samples/
├── 01_trimmed/                # Skewer adapter-trimmed reads
├── 02_filtered/               # BBDuk host-filtered reads
├── 03_assembly/               # metaSPAdes output
│   ├── H_coassembly/
│   │   ├── contigs.fasta
│   │   └── scaffolds.fasta
│   └── Y_coassembly/
├── 04_mapping/                # Bowtie2 BAM files + depth matrices
├── 05_binning/
│   ├── metabat2/
│   ├── maxbin2/
│   ├── concoct/
│   └── dastool/
├── 06_refinement/             # Anvi'o databases + CheckM2 results
├── 07_annotation/             # COG14 + Pfam results
└── 08_resistome/              # RGI + ABRicate + ARGS-OAP outputs
```

---

## 5. Stage 1 — Read Quality Control (Skewer)

**Tool:** [Skewer](https://github.com/relipmoc/skewer)
**Script:** `adapter_qual_filt_eng.md`

Adapter sequences are trimmed and quality filtering is applied to raw paired-end reads.

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-x` | adapter FASTA | Adapter sequences to trim |
| `-q` | 20–30 | Phred quality threshold (3' end) |
| `-l` | 30 | Minimum read length after trimming |
| `-t` | 16 | Threads |
| `-z` | — | Output gzip-compressed FASTQ |
| `-o` | output dir | Output directory |

### Example command

```bash
for sample in /path/to/raw/*_R1.fastq.gz; do
  base=$(basename ${sample} _R1.fastq.gz)
  skewer \
    -x adapters.fa \
    -q 25 \
    -l 30 \
    -t 16 \
    -z \
    -o /path/to/01_trimmed/${base} \
    ${sample} \
    /path/to/raw/${base}_R2.fastq.gz
done
```

### Output

- `*-trimmed-pair1.fastq.gz` — Trimmed forward reads
- `*-trimmed-pair2.fastq.gz` — Trimmed reverse reads
- `*-trimmed.log` — Trimming statistics

---

## 6. Stage 2 — Host Read Removal (BBDuk)

**Tool:** [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)
**Script:** `filter_bbduk_README_en.md`

Reads matching the host genome (*Macrocystis pyrifera* in this project) are removed using k-mer-based matching.

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `k` | 27 | K-mer size |
| `hdist` | 1 | Hamming distance (allows 1 mismatch) |
| `ref` | host FASTA | Reference genome for filtering |
| `-Xmx` | 150g | Maximum Java heap memory |
| `t` | 16 | Threads |

### Example command

```bash
for sample in /path/to/01_trimmed/*-trimmed-pair1.fastq.gz; do
  base=$(basename ${sample} -trimmed-pair1.fastq.gz)
  bbduk.sh \
    -Xmx150g \
    in1=${sample} \
    in2=/path/to/01_trimmed/${base}-trimmed-pair2.fastq.gz \
    out1=/path/to/02_filtered/${base}_R1_filtered.fastq.gz \
    out2=/path/to/02_filtered/${base}_R2_filtered.fastq.gz \
    ref=macrocystis_pyrifera.fasta \
    k=27 \
    hdist=1 \
    t=16
done
```

### Output

- `*_R1_filtered.fastq.gz` / `*_R2_filtered.fastq.gz` — Host-depleted reads

---

## 7. Stage 3 — Co-Assembly (metaSPAdes)

**Tool:** [SPAdes](https://github.com/ablab/spades) (`--meta` mode)
**Script:** `coensamble_spades_EN.md`

All filtered reads from each sample group (H and Y separately) are co-assembled into a single metagenome.

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--meta` | — | Metagenomic assembly mode |
| `-k` | 21,33,55,77,99,127 | K-mer sizes |
| `-t` | 30 | Threads |
| `-m` | 200 | Memory limit (GB) |
| `--pe-1 1..N` | R1 files | Forward reads from N samples |
| `--pe-2 1..N` | R2 files | Reverse reads from N samples |

### Example command (H co-assembly, 6 samples)

```bash
spades.py --meta \
  -k 21,33,55,77,99,127 \
  -t 30 \
  -m 200 \
  --pe1-1 H1_R1.fastq.gz --pe1-2 H1_R2.fastq.gz \
  --pe2-1 H2_R1.fastq.gz --pe2-2 H2_R2.fastq.gz \
  --pe3-1 H3_R1.fastq.gz --pe3-2 H3_R2.fastq.gz \
  --pe4-1 H4_R1.fastq.gz --pe4-2 H4_R2.fastq.gz \
  --pe5-1 H5_R1.fastq.gz --pe5-2 H5_R2.fastq.gz \
  --pe6-1 H6_R1.fastq.gz --pe6-2 H6_R2.fastq.gz \
  -o 03_assembly/H_coassembly/
```

### Output

- `contigs.fasta` — Assembled contigs
- `scaffolds.fasta` — Scaffolded sequences
- `spades.log` — Assembly log

> **Tip:** Use `nohup` to run in the background: `nohup spades.py ... > spades.out 2>&1 &`

---

## 8. Stage 4 — Read Mapping (Bowtie2)

**Tool:** [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) + SAMtools
**Script:** `bowtieY256.md`

Filtered reads from each sample are mapped back to the co-assembled contigs to generate per-sample coverage profiles required for binning.

### Steps

1. Build Bowtie2 index from the co-assembly contigs
2. Map each sample's reads to the index
3. Convert, sort, and index the BAM files
4. Generate a depth/coverage matrix using `jgi_summarize_bam_contig_depths`

### Example commands

```bash
# 1. Build index
bowtie2-build contigs.fasta contigs_index

# 2. Map each sample
for sample in H1 H2 H3 H4 H5 H6; do
  bowtie2 \
    -x contigs_index \
    -1 02_filtered/${sample}_R1_filtered.fastq.gz \
    -2 02_filtered/${sample}_R2_filtered.fastq.gz \
    --sensitive-local \
    -p 30 \
    -S 04_mapping/${sample}.sam

  # Convert to sorted BAM
  samtools view -bS 04_mapping/${sample}.sam | \
    samtools sort -o 04_mapping/${sample}_sorted.bam
  samtools index 04_mapping/${sample}_sorted.bam
  rm 04_mapping/${sample}.sam
done

# 3. Generate depth matrix
jgi_summarize_bam_contig_depths \
  --outputDepth 04_mapping/depth_matrix.txt \
  04_mapping/*_sorted.bam
```

### Output

- `*_sorted.bam` + `.bai` index — Per-sample sorted BAM files
- `depth_matrix.txt` — Coverage matrix (contigs × samples), used as input for all binners

---

## 9. Stage 5 — Metagenomic Binning

Three complementary binning algorithms are run independently. Results are later combined by DAS_Tool.

### 5a. MetaBAT2

**Script:** `00_bin_Metabat2_EN.md`

MetaBAT2 bins contigs using differential coverage and tetranucleotide frequency. Three minimum contig size thresholds are tested.

```bash
# Run MetaBAT2 at three minimum contig thresholds
for min_contig in 1500 2000 2500; do
  metabat2 \
    -i 03_assembly/H_coassembly/contigs.fasta \
    -a 04_mapping/depth_matrix.txt \
    -o 05_binning/metabat2/min${min_contig}/bin \
    -m ${min_contig} \
    -t 30 \
    --saveCls
done
```

### 5b. MaxBin2

**Script:** `00_MaxBin2_EN_ES.md`

MaxBin2 uses marker gene expectation-maximization to assign contigs to bins.

```bash
# Generate depth file for MaxBin2 (single-column format)
awk '{print $1"\t"$3}' 04_mapping/depth_matrix.txt > 04_mapping/maxbin_depth.txt

# Run MaxBin2 at three minimum contig sizes
for min_contig in 1500 2000 2500; do
  run_MaxBin.pl \
    -contig 03_assembly/H_coassembly/contigs.fasta \
    -abund 04_mapping/maxbin_depth.txt \
    -out 05_binning/maxbin2/min${min_contig}/bin \
    -min_contig_length ${min_contig} \
    -prob_threshold 0.9 \
    -thread 60
done
```

### 5c. CONCOCT

**Script:** `concoct_ES.md`
**Environment:** `concoct-1p1` conda environment

CONCOCT clusters contigs by splitting them into sub-fragments to improve coverage uniformity.

```bash
# Activate the CONCOCT environment
conda activate concoct-1p1

# Step 1: Fragment contigs (10,000 bp chunks)
cut_up_fasta.py contigs.fasta \
  -c 10000 -o 0 \
  --merge_last \
  -b 05_binning/concoct/contigs_10K.bed \
  > 05_binning/concoct/contigs_10K.fasta

# Step 2: Generate coverage table
concoct_coverage_table.py \
  05_binning/concoct/contigs_10K.bed \
  04_mapping/*_sorted.bam \
  > 05_binning/concoct/coverage_table.tsv

# Step 3: Run CONCOCT
concoct \
  --composition_file 05_binning/concoct/contigs_10K.fasta \
  --coverage_file 05_binning/concoct/coverage_table.tsv \
  -b 05_binning/concoct/output/ \
  -t 40

# Step 4: Merge sub-contig clusters back
merge_cutup_clustering.py \
  05_binning/concoct/output/clustering_gt1000.csv \
  > 05_binning/concoct/output/clustering_merged.csv

# Step 5: Extract bin FASTA files
mkdir 05_binning/concoct/bins/
extract_fasta_bins.py \
  contigs.fasta \
  05_binning/concoct/output/clustering_merged.csv \
  --output_path 05_binning/concoct/bins/
```

### 5d. DAS_Tool (bin refinement)

**Script:** `dastool.md`

DAS_Tool selects the best non-redundant set of bins across all three binners using a scoring function based on single-copy gene frequencies.

```bash
# Prepare scaffold-to-bin tables for each binner
# (convert bin FASTA files into 2-column TSV: scaffold_id TAB bin_id)
Fasta_to_Scaffolds2Bin.sh -i 05_binning/metabat2/min1500/ -e fa > metabat_scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i 05_binning/maxbin2/min1500/  -e fasta > maxbin_scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i 05_binning/concoct/bins/     -e fa > concoct_scaffolds2bin.tsv

# Run DAS_Tool (test multiple score thresholds)
for score_thresh in 0.5 0.4 0.3; do
  DAS_Tool \
    -i metabat_scaffolds2bin.tsv,maxbin_scaffolds2bin.tsv,concoct_scaffolds2bin.tsv \
    -l MetaBAT2,MaxBin2,CONCOCT \
    -c contigs.fasta \
    -o 05_binning/dastool/score${score_thresh}/DASTool \
    --score_threshold ${score_thresh} \
    --search_engine blast \
    --write_bins 1 \
    -t 10 \
    --debug
done
```

> **Score threshold guidance:** Start with `0.5` (strict); lower to `0.4` or `0.3` if too few bins pass. Each run produces its own output directory for comparison.

---

## 10. Stage 6 — MAG Refinement & Quality (Anvi'o + CheckM2)

**Scripts:** `anvi-refine_H_samples.md`, `anvi-refine_Y_samples.md`, `README_refining.md`
**Tools:** [Anvi'o ≥7.0](https://anvio.org/), [CheckM2](https://github.com/chklovski/CheckM2)

Anvi'o is used for interactive visual curation of bins. Each bin is inspected, and contigs are manually reassigned if needed based on GC content, differential coverage, and single-copy gene (SCG) duplication patterns.

### Step-by-step workflow

```bash
# 1. Create a contigs database
anvi-gen-contigs-database \
  -f contigs.fasta \
  -o contigs.db \
  -T 40

# 2. Annotate HMMs (bacterial + archaeal single-copy genes)
anvi-run-hmms -c contigs.db -T 40

# 3. Assign SCG taxonomy
anvi-run-scg-taxonomy -c contigs.db -T 40

# 4. Profile each sample's BAM file
for sample in H1 H2 H3 H4 H5 H6; do
  anvi-profile \
    -i 04_mapping/${sample}_sorted.bam \
    -c contigs.db \
    -o 06_refinement/profiles/${sample}/ \
    -T 50 \
    --min-contig-length 1000
done

# 5. Merge all profiles
anvi-merge \
  06_refinement/profiles/*/PROFILE.db \
  -c contigs.db \
  -o 06_refinement/merged_profile/ \
  --enforce-hierarchical-clustering

# 6. Import DAS_Tool bins as a collection
anvi-import-collection \
  metabat_scaffolds2bin.tsv \
  -c contigs.db \
  -p 06_refinement/merged_profile/PROFILE.db \
  --collection-name DAS_Tool_bins \
  --contigs-mode

# 7. Launch interactive refinement for each bin
anvi-refine \
  -c contigs.db \
  -p 06_refinement/merged_profile/PROFILE.db \
  -C DAS_Tool_bins \
  -b BIN_NAME

# 8. Summarize refined collection
anvi-summarize \
  -c contigs.db \
  -p 06_refinement/merged_profile/PROFILE.db \
  -C REFINED_COLLECTION \
  -o 06_refinement/summary/

# 9. Assess quality with CheckM2
checkm2 predict \
  --input 06_refinement/summary/bin_by_bin/ \
  --output-directory 06_refinement/checkm2/ \
  -x fa \
  --threads 30
```

### MAG quality thresholds (MIMAG standard)

| Category | Completeness | Contamination |
|----------|-------------|---------------|
| High Quality (HQ) | ≥90% | <5% |
| Medium Quality (MQ) | ≥50% | <10% |
| Low Quality (LQ) | <50% | — |

> **Refinement tips:**
> - Inspect GC content vs. coverage plots in the Anvi'o interactive interface
> - Split bins with >10% contamination by removing outlier contigs
> - Use the "redundancy" panel to identify duplicated single-copy genes
> - The `anvi-refine_H_samples.md` and `anvi-refine_Y_samples.md` scripts contain sample-group-specific details

---

## 11. Stage 7 — Functional Annotation (COG14 + Pfam)

**Script:** `COG14_Pfam_anvio_pipeline_plasmids.md`
**Tools:** DIAMOND, HMMER (hmmscan), Anvi'o

Predicted protein sequences from each MAG are searched against COG14 (protein clusters) and Pfam-A (domain profiles).

### Step-by-step workflow

```bash
# 1. Export gene calls (amino acid sequences) from Anvi'o
anvi-get-sequences-for-gene-calls \
  -c contigs.db \
  --get-aa-sequences \
  -o 07_annotation/gene-calls.faa

# 2. Run DIAMOND against COG14
diamond blastp \
  -q 07_annotation/gene-calls.faa \
  -d /path/to/COG14/db.dmnd \
  -o 07_annotation/diamond_cog14.txt \
  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  -e 1e-5 \
  -p 30

# 3. Format and fix encoding issues
# (ensure UTF-8, remove non-ASCII artifacts before import)
iconv -f latin1 -t utf-8 07_annotation/diamond_cog14.txt \
  > 07_annotation/diamond_cog14_utf8.txt

# 4. Import COG14 annotations into Anvi'o
anvi-import-functions \
  -c contigs.db \
  -i 07_annotation/diamond_cog14_utf8.txt \
  --annotation-source COG14

# 5. Run HMMER against Pfam-A v32
hmmscan \
  --tblout 07_annotation/pfam_results.tbl \
  --cut_ga \
  --cpu 30 \
  /path/to/Pfam-A.hmm \
  07_annotation/gene-calls.faa

# 6. Format Pfam output for Anvi'o and import
# (use the parsing script provided in the repository)
anvi-import-functions \
  -c contigs.db \
  -i 07_annotation/pfam_formatted.txt \
  --annotation-source Pfam_v32
```

> **Note for plasmid sequences:** An identical annotation workflow applies to plasmid contigs. See `COG14_Pfam_anvio_pipeline_plasmids.md` for the plasmid-specific version.

---

## 12. Stage 8 — Resistome Analysis

**Script:** `Resistome.md`
**Tools:** RGI, ABRicate, ARGS-OAP

Three complementary tools are used to identify antibiotic resistance genes (ARGs) to maximize sensitivity.

### Tool 1: RGI (Resistance Gene Identifier)

Uses the CARD database (v4.0.0) with DIAMOND as the search engine.

```bash
# Load CARD database
rgi load \
  --card_json /path/to/card.json \
  --local

# Run RGI on MAG protein sequences
rgi main \
  --input_sequence 07_annotation/gene-calls.faa \
  --output_file 08_resistome/rgi_output \
  --input_type protein \
  --alignment_tool DIAMOND \
  --num_threads 60 \
  --local \
  --clean
```

### Tool 2: ABRicate

Screens nucleotide sequences against multiple databases simultaneously.

```bash
for db in argannot card ncbi ncbibetalactamase plasmidfinder resfinder; do
  abricate \
    --db ${db} \
    --minid 80 \
    --mincov 80 \
    contigs.fasta \
    > 08_resistome/abricate_${db}.tab
done

# Summarize all ABRicate results
abricate --summary 08_resistome/abricate_*.tab > 08_resistome/abricate_summary.tab
```

### Tool 3: ARGS-OAP

Provides normalized ARG abundance profiles.

```bash
# Stage one: extract candidate ARG reads
python /path/to/args_oap/stage_one.py \
  -i 02_filtered/ \
  -o 08_resistome/args_oap/ \
  -t 30

# Stage two: classify and quantify ARGs
python /path/to/args_oap/stage_two.py \
  -i 08_resistome/args_oap/ \
  -o 08_resistome/args_oap_results/ \
  -t 30
```

---

## 13. Output Summary

After completing the full pipeline, the following key files are available for downstream analysis:

```
06_refinement/
├── checkm2/quality_report.tsv          # Completeness & contamination per MAG
└── summary/bin_by_bin/                 # Individual MAG FASTA files

07_annotation/
├── diamond_cog14_utf8.txt              # COG14 annotations
└── pfam_results.tbl                    # Pfam domain annotations
    (Both also imported into contigs.db for Anvi'o queries)

08_resistome/
├── rgi_output.txt                      # CARD-based ARGs (RGI)
├── abricate_*.tab                      # Multi-database ARG screening
├── abricate_summary.tab                # Combined ABRicate summary
└── args_oap_results/                   # Normalized ARG abundances
```

---

## 14. Citation

If you use this pipeline, please cite the tools used in each stage:

- **Skewer:** Jiang et al. (2014) BMC Bioinformatics 15:182
- **BBDuk:** Bushnell B. (2014) sourceforge.net/projects/bbmap
- **metaSPAdes:** Nurk et al. (2017) Genome Research 27:824–834
- **Bowtie2:** Langmead & Salzberg (2012) Nature Methods 9:357–359
- **MetaBAT2:** Kang et al. (2019) PeerJ 7:e7359
- **MaxBin2:** Wu et al. (2016) Bioinformatics 32:605–607
- **CONCOCT:** Alneberg et al. (2014) Nature Methods 11:1144–1146
- **DAS_Tool:** Sieber et al. (2018) Nature Microbiology 3:836–843
- **Anvi'o:** Eren et al. (2021) Nature Microbiology 6:727–740
- **CheckM2:** Chklovski et al. (2023) Nature Methods 20:1203–1212
- **DIAMOND:** Buchfink et al. (2021) Nature Methods 18:366–368
- **HMMER:** Eddy (2011) PLOS Computational Biology 7:e1002195
- **RGI / CARD:** Alcock et al. (2023) Nucleic Acids Research 51:D690–D699
- **ABRicate:** Seemann T. (2020) github.com/tseemann/abricate
- **ARGS-OAP:** Yin et al. (2018) Water Research 138:70–77

---

*Pipeline developed by Mario Moreno — FOVI Project*
*Repository: [github.com/MorenoMario/FOVI_MAGs_pipeline](https://github.com/MorenoMario/FOVI_MAGs_pipeline)*
