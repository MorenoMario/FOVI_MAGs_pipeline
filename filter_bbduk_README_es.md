# Filtrado de lecturas contra referencia de *Macrocystis pyrifera* (BBDuk)

Este paso filtra lecturas **paired-end** contra el genoma de *Macrocystis pyrifera* usando `bbduk.sh`. Soporta:
1) **Conservar ambos conjuntos** (hits y no-hits) y 2) **Conservar solo no-hits** (descartar hits).  
El script busca **recursivamente** dentro de `INPUT_DIR` archivos `*-trimmed-pair1.fastq.gz` (o `.fastq`), por lo que funciona cuando cada muestra está en su propia carpeta, p.ej.:

```
/media/mmorenos/disk3/DSG123-selected/rawReads/00_skewer/20250604_PamelaFernandez_DSG_H1_S7/
  ├─ 20250604_PamelaFernandez_DSG_H1_S7-trimmed-pair1.fastq.gz
  └─ 20250604_PamelaFernandez_DSG_H1_S7-trimmed-pair2.fastq.gz
```

---

## Requisitos

- Linux/macOS con `bash`
- [BBTools / BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) en el `PATH`
- `gzip`

Comprobar:
```bash
bbduk.sh --version
```

---

## Script (`filter_host_bbduk.sh`)

```bash
#!/usr/bin/env bash
set -euo pipefail

# ===============================
# BBDuk host/alga filter (paired)
# ===============================
# Usage / Uso:
#   ./filter_host_bbduk.sh -i INPUT_DIR -o OUTPUT_DIR -r REFERENCE.fa [-t 16] [-m 150g] [-k 27] [-d 1] [--force] [--no-keep-hits]
#
# Outputs per sample / Salidas por muestra:
#   <sample>_macroHit_R1.fastq.gz      # reads that MATCH the algal reference (R1) / lecturas que MATCHEAN al alga (R1)
#   <sample>_macroHit_R2.fastq.gz      # ... (R2)
#   <sample>_non-macroHit_R1.fastq.gz  # reads that DO NOT match (R1) / lecturas que NO matchean (R1)
#   <sample>_non-macroHit_R2.fastq.gz  # ... (R2)
#   logs/<sample>.bbduk.stats.txt
#   manifest_bbduk.tsv (accumulative / acumulativo)

INPUT_DIR=""
OUTPUT_DIR=""
REFERENCE=""
THREADS=8
MEMORY="64g"
KMER=27
HDIST=1
FORCE=false
KEEP_HITS=true  # set to false with --no-keep-hits

print_usage() {
  cat <<'EOF'
Usage: filter_host_bbduk.sh -i INPUT_DIR -o OUTPUT_DIR -r REFERENCE.fa [-t THREADS] [-m MEM] [-k KMER] [-d HDIST] [--force] [--no-keep-hits]

Options:
  -i    Input directory with Skewer outputs (FASTQ(.gz))
  -o    Output directory
  -r    Algal reference FASTA (e.g., GCA_031763025.fna)
  -t    Threads (default: 8)
  -m    Java memory for BBTools (default: 64g)  # passed as -Xmx
  -k    k-mer size for matching (default: 27)
  -d    hdist (Hamming distance) for k-mers (default: 1)
  --force         Reprocess even if outputs exist
  --no-keep-hits  Do not write matched reads (send to /dev/null)

Notas/Notes:
- The script searches *recursively* within INPUT_DIR for files named "*-trimmed-pair1.fastq.gz" (or .fastq).
- BBduk auto-detects gzip via file extension.
- Place -Xmx before other BBduk parameters.
EOF
}

# --- Parse long flags first ---
long_opts=()
for arg in "$@"; do
  case "$arg" in
    --force) FORCE=true ;;
    --no-keep-hits) KEEP_HITS=false ;;
    --) shift; break ;;  # end of long options
    --*) echo "Invalid long option: $arg" >&2; print_usage; exit 1 ;;
    *) long_opts+=("$arg") ;;
  esac
done
set -- "${long_opts[@]}"

# --- Short options ---
while getopts ":i:o:r:t:m:k:d:h" opt; do
  case "$opt" in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) REFERENCE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MEMORY="$OPTARG" ;;
    k) KMER="$OPTARG" ;;
    d) HDIST="$OPTARG" ;;
    h) print_usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; print_usage; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; print_usage; exit 1 ;;
  esac
done || true

# --- Validations ---
[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$REFERENCE" ]] && { print_usage; exit 1; }
[[ ! -d "$INPUT_DIR" ]] && { echo "INPUT_DIR does not exist: $INPUT_DIR" >&2; exit 1; }
[[ ! -f "$REFERENCE" ]] && { echo "REFERENCE does not exist: $REFERENCE" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR/logs"
MANIFEST="$OUTPUT_DIR/manifest_bbduk.tsv"
if [[ ! -f "$MANIFEST" ]]; then
  echo -e "sample\tr1_in\tr2_in\tmatched_r1\tmatched_r2\tnonmatched_r1\tnonmatched_r2\tstats\tkmer\thdist\tkeep_hits" > "$MANIFEST"
fi

# --- Helpers ---
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

# --- Discover R1 recursively ---
shopt -s nullglob
mapfile -t R1_FILES < <( \
  find "$INPUT_DIR" -type f \( -name "*-trimmed-pair1.fastq.gz" -o -name "*-trimmed-pair1.fastq" \) \
  | sort
)

if (( ${#R1_FILES[@]} == 0 )); then
  echo "No R1 files like *-trimmed-pair1.fastq(.gz) found under: $INPUT_DIR" >&2
  exit 1
fi

# --- Main loop ---
for R1 in "${R1_FILES[@]}"; do
  R2="$(derive_r2_from_pair1 "$R1" || true)"
  if [[ -z "$R2" ]]; then
    echo "✗ No matching R2 for: $R1. Skipping."
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
  LOG="$OUTPUT_DIR/logs/${SAMPLE}.log"

  # Resume check
  if [[ "$FORCE" = false ]]; then
    if $KEEP_HITS; then
      if [[ -s "$MATCHED1" && -s "$MATCHED2" && -s "$UNMATCHED1" && -s "$UNMATCHED2" ]]; then
        echo "• Outputs exist for $SAMPLE. Use --force to reprocess. Skipping."
        continue
      fi
    else
      if [[ -s "$UNMATCHED1" && -s "$UNMATCHED2" ]]; then
        echo "• Non-hit outputs exist for $SAMPLE. Use --force to reprocess. Skipping."
        continue
      fi
    fi
  fi

  echo "→ BBDuk filtering: $SAMPLE (k=$KMER, hdist=$HDIST, keep_hits=$KEEP_HITS)"
  set +e
  if $KEEP_HITS; then
    bbduk.sh -Xmx"$MEMORY" \
      in="$R1" in2="$R2" \
      outm="$MATCHED1" outm2="$MATCHED2" \
      outu="$UNMATCHED1" outu2="$UNMATCHED2" \
      ref="$REFERENCE" \
      k="$KMER" hdist="$HDIST" \
      t="$THREADS" overwrite=t stats="$STATS" \
      >> "$LOG" 2>&1
  else
    bbduk.sh -Xmx"$MEMORY" \
      in="$R1" in2="$R2" \
      outm=/dev/null outm2=/dev/null \
      outu="$UNMATCHED1" outu2="$UNMATCHED2" \
      ref="$REFERENCE" \
      k="$KMER" hdist="$HDIST" \
      t="$THREADS" overwrite=t stats="$STATS" \
      >> "$LOG" 2>&1
  fi
  status=$?
  set -e

  if [[ $status -ne 0 ]]; then
    echo "✗ BBDuk failed for $SAMPLE (see $LOG)."
    continue
  fi

  # Minimal checks
  if $KEEP_HITS; then
    if [[ ! -s "$MATCHED1" || ! -s "$MATCHED2" || ! -s "$UNMATCHED1" || ! -s "$UNMATCHED2" ]]; then
      echo "✗ Missing outputs for $SAMPLE (see logs)."
      continue
    fi
    echo -e "${SAMPLE}\t${R1}\t${R2}\t${MATCHED1}\t${MATCHED2}\t${UNMATCHED1}\t${UNMATCHED2}\t${STATS}\t${KMER}\t${HDIST}\ttrue" >> "$MANIFEST"
  else
    if [[ ! -s "$UNMATCHED1" || ! -s "$UNMATCHED2" ]]; then
      echo "✗ Missing non-hit outputs for $SAMPLE (see logs)."
      continue
    fi
    echo -e "${SAMPLE}\t${R1}\t${R2}\tNA\tNA\t${UNMATCHED1}\t${UNMATCHED2}\t${STATS}\t${KMER}\t${HDIST}\tfalse" >> "$MANIFEST"
  fi

  echo "✓ Done: $SAMPLE"
done

echo "== Finished =="
echo "Manifest: $MANIFEST"

```

> **Permisos**: `chmod +x filter_host_bbduk.sh`

---

## Uso

```bash
./filter_host_bbduk.sh -i INPUT_DIR -o OUTPUT_DIR -r REFERENCE.fa [-t THREADS] [-m MEM] [-k KMER] [-d HDIST] [--force] [--no-keep-hits]
```

**Opciones clave:**
- `-i`  Directorio raíz que contiene las carpetas por muestra con los `*-trimmed-pair1.fastq.gz` / `*-trimmed-pair2.fastq.gz`
- `-o`  Directorio de salida
- `-r`  FASTA de referencia del alga (ej.: `GCA_031763025.fna`)
- `-t`  Hilos (por defecto: 8)
- `-m`  Memoria Java para BBduk (por defecto: `64g`) → se pasa como `-Xmx`
- `-k`  k-mer (por defecto: 27)
- `-d`  hdist para k-mers (por defecto: 1)
- `--no-keep-hits`  No escribir las lecturas que matchean (se envían a `/dev/null`)
- `--force`  Reprocesar aunque existan salidas

---

## Ejemplos

**A) Conservar ambos conjuntos (hits & no-hits):**
```bash
./filter_host_bbduk.sh   -i /media/mmorenos/disk3/DSG123-selected/rawReads/00_skewer/   -o /media/mmorenos/disk3/DSG123-selected/rawReads/01_bbduk_ref/   -r /media/mmorenos/disk3/refs/GCA_031763025.fna   -t 32 -m 150g -k 27 -d 1
```

**B) Conservar solo no-hits (descartar hits):**
```bash
./filter_host_bbduk.sh   -i /media/mmorenos/disk3/DSG123-selected/rawReads/00_skewer/   -o /media/mmorenos/disk3/DSG123-selected/rawReads/01_bbduk_ref_nonhits/   -r /media/mmorenos/disk3/refs/GCA_031763025.fna   -t 32 -m 150g -k 27 -d 1   --no-keep-hits
```

**C) Ejecutar en segundo plano con `nohup`:**
```bash
nohup ./filter_host_bbduk.sh   -i /media/mmorenos/disk3/DSG123-selected/rawReads/00_skewer/   -o /media/mmorenos/disk3/DSG123-selected/rawReads/01_bbduk_ref/   -r /media/mmorenos/disk3/refs/GCA_031763025.fna   -t 60 -m 150g   > filter_bbduk.out 2>&1 &
tail -f filter_bbduk.out
```

---

## Salidas

Por muestra se generan:
- `<sample>_macroHit_R1.fastq.gz`, `<sample>_macroHit_R2.fastq.gz` (si no usas `--no-keep-hits`)
- `<sample>_non-macroHit_R1.fastq.gz`, `<sample>_non-macroHit_R2.fastq.gz`
- `logs/<sample>.bbduk.stats.txt`, `logs/<sample>.log`
- `manifest_bbduk.tsv` con columnas:  
  `sample, r1_in, r2_in, matched_r1, matched_r2, nonmatched_r1, nonmatched_r2, stats, kmer, hdist, keep_hits`

---

## Notas

- `-Xmx` debe ir **antes** del resto de parámetros de BBduk.
- Si aparecen falsos positivos, prueba `-k 31` o `-d 0`. Si necesitas más sensibilidad, `-k 23` (con más riesgo de falso positivo).
- El script es reanudable: usa `--force` para reprocesar una muestra.
