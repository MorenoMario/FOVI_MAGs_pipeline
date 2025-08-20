# MetaBAT2: Binning of co-assemblies (Groups Y and H)

This workflow runs **MetaBAT2** on the Y and H co‑assemblies, using coverage files
generated from BAMs (`jgi_summarize_bam_contig_depths`) or by merging
`*_bam_depth.txt` files. The goal is to obtain MAG bins under different minimum
contig length parameters.

## Requirements

- `metabat2`
- `samtools`
- `jgi_summarize_bam_contig_depths` (part of MetaBAT2 package)
- Previous assemblies (e.g., `contigs.fasta` from SPAdes)
- BAM files mapped against the co‑assemblies (generated with Bowtie2)

## Example script (for Y)

```bash
#!/usr/bin/env bash
set -euo pipefail

# =========================
# MetaBAT2 on Y co-assembly
#  - Generates unified coverage (jgi) from BAMs
#  - Alternative: merges existing *bam_depth.txt files
#  - Runs MetaBAT2 with multiple -m values
# =========================

# --- Config ---
BAM_DIR="/datos1/mmorenos/fovi_data/03_bowtie/Y"
DEPTH_DIR="${BAM_DIR}"
CONTIGS="/datos1/mmorenos/fovi_data/02_spades_ass/spades_Y_coassembly/contigs.fasta"

OUT_ROOT="${BAM_DIR}/metabat2_Y"
THREADS=30
MINLEN_LIST=("1500" "2000" "2500")

# Coverage output names
JGI_DEPTH="${OUT_ROOT}/Y_coverage_jgi.tsv"
MERGED_DEPTH="${OUT_ROOT}/Y_coverage_merged.tsv"

# --- Preparation ---
mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/bins"

# --- Collect BAMs ---
mapfile -t BAMS < <(ls -1 "${BAM_DIR}"/*.sorted.bam 2>/dev/null | sort)

# --- Generate coverage with JGI ---
if (( ${#BAMS[@]} > 0 )); then
  jgi_summarize_bam_contig_depths --outputDepth "${JGI_DEPTH}" --minContigLength 1000 "${BAMS[@]}"
  DEPTH_FOR_METABAT="${JGI_DEPTH}"
else
  DEPTH_FOR_METABAT="${MERGED_DEPTH}"
fi

# --- Run MetaBAT2 with different -m values ---
for M in "${MINLEN_LIST[@]}"; do
  OUT_PREFIX="${OUT_ROOT}/bins/metabat2_Y_m${M}"
  metabat2 -i "${CONTIGS}" -a "${DEPTH_FOR_METABAT}" -o "${OUT_PREFIX}" -m "${M}" -t "${THREADS}" --seed 1
done
```

## Version for H

Change `BAM_DIR`, `CONTIGS`, and `OUT_ROOT` replacing `Y` with `H`.
