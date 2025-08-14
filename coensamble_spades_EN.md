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

Key parameters used:
- **metaSPAdes** enabled with `--meta`
- K-mers: `21,33,55,77,99,127`
- Threads: `30`
- Memory: `200 GB`
- `--only-assembler` to skip error correction in this stage
- Temporary directory on fast storage with `--tmp-dir` (optional)

## 3. Running the script
Example execution:
```bash
./coassemble_spades_v2.sh   -i /path/to/01_bbduk_ref   -o /path/to/spades_coassemblies   -t 30 -m 200   --kmers 21,33,55,77,99,127   --meta   --tmp-dir /path/to/tmp
```

## 4. Generated results
The script generates two output folders:
```
spades_coassemblies/spades_H_coassembly/
spades_coassemblies/spades_Y_coassembly/
```
Each folder contains:
- `contigs.fasta` → Assembled contigs
- `scaffolds.fasta` → Assembled scaffolds
- `spades.log` → Full execution log

## 5. Recommended next steps
1. Map original reads back to contigs to calculate coverage.
2. Filter contigs by length (≥ 1–1.5 kb).
3. Perform binning with **MetaBAT2**, **VAMB**, or another tool.
4. Annotate rRNA genes (optional) using **Barrnap** to check for ribosomal gene presence.
