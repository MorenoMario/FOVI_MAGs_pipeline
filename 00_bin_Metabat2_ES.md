# MetaBAT2: Binning de co-ensambles metagenómicos (Grupo Y y H)

Este flujo de trabajo ejecuta **MetaBAT2** sobre los co‑ensambles Y y H, utilizando
archivos de cobertura generados a partir de BAMs (`jgi_summarize_bam_contig_depths`)
o fusionando archivos `*_bam_depth.txt`. El objetivo es obtener bins de MAGs bajo
diferentes parámetros de longitud mínima de contigs.

## Requisitos

- `metabat2`
- `samtools`
- `jgi_summarize_bam_contig_depths` (incluido en MetaBAT2 package)
- Ensambles previos (ej. `contigs.fasta` de SPAdes)
- Archivos BAM mapeados contra los co‑ensambles (generados con Bowtie2)

## Script ejemplo (para Y)

```bash
#!/usr/bin/env bash
set -euo pipefail

# =========================
# MetaBAT2 sobre co-ensamble Y
#  - Genera cobertura unificada (jgi) desde BAMs
#  - Alternativa: fusiona *bam_depth.txt existentes
#  - Ejecuta MetaBAT2 con múltiples -m
# =========================

# --- Config ---
BAM_DIR="/datos1/mmorenos/fovi_data/03_bowtie/Y"
DEPTH_DIR="${BAM_DIR}"   # donde están los *bam_depth.txt (si existieran)
CONTIGS="/datos1/mmorenos/fovi_data/02_spades_ass/spades_Y_coassembly/contigs.fasta"

OUT_ROOT="${BAM_DIR}/metabat2_Y"
THREADS=30
MINLEN_LIST=("1500" "2000" "2500")

# Nombre de las salidas de cobertura
JGI_DEPTH="${OUT_ROOT}/Y_coverage_jgi.tsv"
MERGED_DEPTH="${OUT_ROOT}/Y_coverage_merged.tsv"

# --- Preparación ---
mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/bins"

# --- Recolectar BAMs ---
mapfile -t BAMS < <(ls -1 "${BAM_DIR}"/*.sorted.bam 2>/dev/null | sort)

# --- Generar cobertura con JGI ---
if (( ${#BAMS[@]} > 0 )); then
  jgi_summarize_bam_contig_depths --outputDepth "${JGI_DEPTH}" --minContigLength 1000 "${BAMS[@]}"
  DEPTH_FOR_METABAT="${JGI_DEPTH}"
else
  DEPTH_FOR_METABAT="${MERGED_DEPTH}"
fi

# --- Ejecutar MetaBAT2 con distintos -m ---
for M in "${MINLEN_LIST[@]}"; do
  OUT_PREFIX="${OUT_ROOT}/bins/metabat2_Y_m${M}"
  metabat2 -i "${CONTIGS}" -a "${DEPTH_FOR_METABAT}" -o "${OUT_PREFIX}" -m "${M}" -t "${THREADS}" --seed 1
done
```

## Versión para H

Cambiar `BAM_DIR`, `CONTIGS` y `OUT_ROOT` reemplazando `Y` por `H`.
