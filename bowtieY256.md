bowtieY256.md

```
#!/usr/bin/env bash
set -euo pipefail

# =========================
# CONFIG
# =========================
THREADS=30
INPUT_DIR="/datos1/mmorenos/fovi_data/01_bbduk_ref/00_non-macroHit"

# Co-ensamble (referencia)
REF_Y256="/datos1/mmorenos/fovi_data/02_spades_ass/spades_Y256_coassembly/Y_S256_contigs.fasta"

# Output base (puedes cambiarlo)
OUT_BASE="/datos1/mmorenos/fovi_data/03_mapping"
SUB="Y_256"
OUT="${OUT_BASE}/${SUB}"

# Bowtie2 options
BT2_OPTS="--threads ${THREADS} --reorder --very-sensitive-local -q --no-unal"

# Index prefix
IDX_PREFIX="${OUT}/index/Y256_contigs"

# =========================
# PREP
# =========================
mkdir -p "${OUT}/logs" "${OUT}/index" "${OUT}/bam"

echo "============================================================"
echo "[INFO] Bowtie2 mapping to Y256 co-assembly"
echo "[INFO] INPUT_DIR : ${INPUT_DIR}"
echo "[INFO] REF_Y256   : ${REF_Y256}"
echo "[INFO] OUT        : ${OUT}"
echo "[INFO] THREADS    : ${THREADS}"
echo "============================================================"

# =========================
# CHECKS
# =========================
command -v bowtie2 >/dev/null || { echo "ERROR: bowtie2 no está en PATH"; exit 1; }
command -v bowtie2-build >/dev/null || { echo "ERROR: bowtie2-build no está en PATH"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools no está en PATH (requerido para BAM)"; exit 1; }
[[ -s "${REF_Y256}" ]] || { echo "ERROR: falta REF_Y256: ${REF_Y256}"; exit 1; }

# =========================
# BUILD INDEX (si falta)
# =========================
need_build() {
  local p="$1"
  for s in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
    [[ -s "${p}.${s}" ]] || return 0
  done
  return 1
}

if need_build "${IDX_PREFIX}"; then
  echo "[INFO] Construyendo índice Bowtie2..."
  bowtie2-build --threads "${THREADS}" "${REF_Y256}" "${IDX_PREFIX}" \
    > "${OUT}/logs/build_Y256.log" 2>&1
  echo "[INFO] Índice listo: ${IDX_PREFIX}.*"
else
  echo "[INFO] Índice ya existe: ${IDX_PREFIX}.*"
fi

# =========================
# DISCOVER INPUTS (Y1..Y6)
# =========================
shopt -s nullglob

mapfile -t R1S < <(ls -1 "${INPUT_DIR}"/20250604_PamelaFernandez_DSG_Y{1,2,3,4,5,6}_S*_non-macroHit_R1.fastq.gz | sort)
(( ${#R1S[@]} )) || { echo "ERROR: no encontré R1 en ${INPUT_DIR}"; exit 1; }

echo "[INFO] Encontradas ${#R1S[@]} librerías R1"

# =========================
# MAPPING LOOP
# =========================
for R1 in "${R1S[@]}"; do
  base="$(basename "${R1}" _non-macroHit_R1.fastq.gz)"
  R2="${INPUT_DIR}/${base}_non-macroHit_R2.fastq.gz"

  if [[ ! -s "${R2}" ]]; then
    echo "[WARN] falta R2 para ${base}. Se omite."
    continue
  fi

  LOG="${OUT}/logs/${base}.bowtie2.log"
  BAM="${OUT}/bam/${base}.sorted.bam"

  # Read group (útil para downstream)
  RG="--rg-id ${base} --rg SM:${base} --rg PL:ILLUMINA"

  echo "------------------------------------------------------------"
  echo "[INFO] Mapping: ${base}"
  echo "[INFO] R1: ${R1}"
  echo "[INFO] R2: ${R2}"
  echo "[INFO] BAM: ${BAM}"
  echo "------------------------------------------------------------"

  bowtie2 ${BT2_OPTS} ${RG} -x "${IDX_PREFIX}" -1 "${R1}" -2 "${R2}" 2> "${LOG}" \
    | samtools view -@ "${THREADS}" -bS - \
    | samtools sort -@ "${THREADS}" -o "${BAM}" -

  samtools index -@ "${THREADS}" "${BAM}"
done

echo "============================================================"
echo "[DONE] Salida final:"
echo "  - BAMs + BAI : ${OUT}/bam/"
echo "  - logs       : ${OUT}/logs/"
echo "  - index      : ${OUT}/index/"
echo "============================================================"

´´´```

