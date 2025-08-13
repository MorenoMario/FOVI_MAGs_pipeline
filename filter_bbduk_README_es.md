# Filtrado de Secuencias de Alga con BBDuk

Este script filtra lecturas de secuenciación que provienen de un genoma de referencia (por ejemplo, un alga como *Macrocystis pyrifera*) usando **BBDuk**.  
Permite guardar las lecturas coincidentes y no coincidentes en archivos separados.

## Uso

```bash
./filter_host_bbduk.sh -i ./00_skewer -o ./01_ski_ref -r GCA_031763025.fna -t 16 -m 150g --force
```

## Estructura esperada de entrada

Ejemplo:
```
./00_skewer/20250604_PamelaFernandez_DSG_H1_S7/20250604_PamelaFernandez_DSG_H1_S7-trimmed-pair1.fastq.gz
./00_skewer/20250604_PamelaFernandez_DSG_H1_S7/20250604_PamelaFernandez_DSG_H1_S7-trimmed-pair2.fastq.gz
...
```

## Script completo

```bash
#!/usr/bin/env bash
set -euo pipefail

# ===============================
# BBDuk host/alga filter (paired)
# ===============================
# Uso:
#   ./filter_host_bbduk.sh -i INPUT_DIR -o OUTPUT_DIR -r REFERENCE.fa [-t 16] [-m 150g] [--force]
#
# Salidas por muestra:
#   <sample>_macroHit_R1.fastq.gz      # lecturas que MATCHEAN al alga (R1)
#   <sample>_macroHit_R2.fastq.gz      # ... (R2)
#   <sample>_non-macroHit_R1.fastq.gz  # lecturas que NO matchean (R1)
#   <sample>_non-macroHit_R2.fastq.gz  # ... (R2)
#   logs/<sample>.bbduk.stats.txt
#   manifest_bbduk.tsv (acumulativo)

INPUT_DIR=""
OUTPUT_DIR=""
REFERENCE=""
THREADS=8
MEMORY="64g"
FORCE=false

print_usage() {
  cat <<EOF
Uso: $(basename "$0") -i INPUT_DIR -o OUTPUT_DIR -r REFERENCE.fa [-t THREADS] [-m MEM] [--force]

Opciones:
  -i   Directorio con FASTQ(.gz) de entrada (puede tener subcarpetas)
  -o   Directorio de salida
  -r   FASTA de referencia del alga (ej. GCA_031763025.fna)
  -t   Hilos (por defecto: 8)
  -m   Memoria para Java (por defecto: 64g)  # se pasa como -Xmx
  --force  Reprocesar aunque existan salidas
EOF
}

# Parseo de --force antes de getopts
long_opts=()
for arg in "$@"; do
  case "$arg" in
    --force) FORCE=true ;;
    *) long_opts+=("$arg") ;;
  esac
done
set -- "${long_opts[@]}"

# Flags cortos
while getopts ":i:o:r:t:m:h" opt; do
  case $opt in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) REFERENCE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MEMORY="$OPTARG" ;;
    h) print_usage; exit 0 ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; print_usage; exit 1 ;;
    :) echo "La opción -$OPTARG requiere un argumento." >&2; print_usage; exit 1 ;;
  esac
done || true

# Validaciones
[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$REFERENCE" ]] && { print_usage; exit 1; }
[[ ! -d "$INPUT_DIR" ]] && { echo "No existe INPUT_DIR: $INPUT_DIR" >&2; exit 1; }
[[ ! -f "$REFERENCE" ]] && { echo "No existe REFERENCE: $REFERENCE" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR/logs"
MANIFEST="$OUTPUT_DIR/manifest_bbduk.tsv"
if [[ ! -f "$MANIFEST" ]]; then
  echo -e "sample	r1_in	r2_in	matched_r1	matched_r2	nonmatched_r1	nonmatched_r2	stats" > "$MANIFEST"
fi

# === Funciones ===
derive_r2_from_pair1() {
  local r1="$1"
  local r2=""
  if [[ "$r1" == *"-trimmed-pair1.fastq.gz" ]]; then
    r2="${r1/-trimmed-pair1.fastq.gz/-trimmed-pair2.fastq.gz}"
  elif [[ "$r1" == *"-trimmed-pair1.fastq" ]]; then
    r2="${r1/-trimmed-pair1.fastq/-trimmed-pair2.fastq}"
  fi
  [[ -f "$r2" ]] && { echo "$r2"; return 0; }
  echo ""
  return 1
}

sample_prefix() {
  local p="$1"
  local b
  b="$(basename "$p")"
  b="${b%-trimmed-pair1.fastq.gz}"
  b="${b%-trimmed-pair1.fastq}"
  echo "$b"
}

# === Buscar R1 ===
shopt -s nullglob
mapfile -t R1_FILES < <(find "$INPUT_DIR" -type f \
  \( -name "*-trimmed-pair1.fastq.gz" -o -name "*-trimmed-pair1.fastq" \) | sort)

if (( ${#R1_FILES[@]} == 0 )); then
  echo "No se encontraron R1 tipo *-trimmed-pair1.fastq(.gz) en: $INPUT_DIR" >&2
  exit 1
fi

# === Loop principal ===
for R1 in "${R1_FILES[@]}"; do
  R2="$(derive_r2_from_pair1 "$R1" || true)"
  if [[ -z "$R2" ]]; then
    echo "✗ No se encontró R2 para: $R1. Saltando."
    continue
  fi

  SAMPLE="$(sample_prefix "$R1")"
  OUT_PREFIX="$OUTPUT_DIR/${SAMPLE}"
  mkdir -p "$OUTPUT_DIR"

  MATCHED1="${OUT_PREFIX}_macroHit_R1.fastq.gz"
  MATCHED2="${OUT_PREFIX}_macroHit_R2.fastq.gz"
  UNMATCHED1="${OUT_PREFIX}_non-macroHit_R1.fastq.gz"
  UNMATCHED2="${OUT_PREFIX}_non-macroHit_R2.fastq.gz"
  STATS="$OUTPUT_DIR/logs/${SAMPLE}.bbduk.stats.txt"

  if [[ "$FORCE" = false && -s "$MATCHED1" && -s "$MATCHED2" && -s "$UNMATCHED1" && -s "$UNMATCHED2" ]]; then
    echo "• Ya existen salidas para $SAMPLE. Usa --force para reprocesar. Saltando."
    continue
  fi

  echo "→ Filtrando (BBDuk) contra referencia: $SAMPLE"
  set +e
  bbduk.sh -Xmx"$MEMORY" \
    in="$R1" in2="$R2" \
    outm="$MATCHED1" outm2="$MATCHED2" \
    outu="$UNMATCHED1" outu2="$UNMATCHED2" \
    ref="$REFERENCE" \
    k=27 hdist=1 \
    threads="$THREADS" \
    overwrite=t stats="$STATS" \
    >> "$OUTPUT_DIR/logs/${SAMPLE}.log" 2>&1
  status=$?
  set -e

  if [[ $status -ne 0 ]]; then
    echo "✗ BBDuk falló para $SAMPLE (ver $OUTPUT_DIR/logs/${SAMPLE}.log)."
    continue
  fi

  if [[ ! -s "$MATCHED1" || ! -s "$MATCHED2" || ! -s "$UNMATCHED1" || ! -s "$UNMATCHED2" ]]; then
    echo "✗ Faltan salidas para $SAMPLE (revisa logs)."
    continue
  fi

  echo -e "${SAMPLE}	${R1}	${R2}	${MATCHED1}	${MATCHED2}	${UNMATCHED1}	${UNMATCHED2}	${STATS}" >> "$MANIFEST"
  echo "✓ Listo: $SAMPLE"
done

echo "== Terminado =="
echo "Manifiesto BBDuk: $MANIFEST"

```
