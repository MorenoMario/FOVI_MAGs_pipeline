# Skewer PE Trimming – Illumina (TruSeq/Nextera)

Script robusto para **recorte de adaptadores** y **(opcional) trimming por calidad** en lecturas Illumina **paired-end**.  
Detecta automáticamente nombres tipo `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` (formato BCL2FASTQ) y otros formatos comunes.

- Genera **salidas comprimidas** (`.fastq.gz`)
- Crea **logs por muestra** y un **manifest.tsv** para encadenar pasos posteriores
- Reanudable: salta muestras ya procesadas (a menos que uses `--force`)

---

## Requisitos

- Linux/macOS con `bash`
- [`skewer`](https://github.com/relipmoc/skewer) ≥ 0.2.2 en el `PATH`
- `gzip`, `coreutils`

Comprueba:
```bash
skewer --version
```

---

## Instalación

Coloca estos archivos en tu repo/proyecto:

```
.
├── trim_skewer.sh
└── adapt.fa
```

Da permisos de ejecución:

```bash
chmod +x trim_skewer.sh
```

---

## Uso

```bash
./trim_skewer.sh -i INPUT_DIR -o OUTPUT_DIR -a ADAPTERS.fa [-t THREADS] [-q QUAL] [--force]
```

**Opciones:**
- `-i`  Directorio con FASTQ(.gz) de entrada  
- `-o`  Directorio de salida  
- `-a`  Archivo FASTA con adaptadores (para `skewer -x`)  
- `-t`  Hilos (por defecto: 8)  
- `-q`  **Umbral de calidad Phred** para trimming (por defecto **0** = sin recorte por calidad)  
- `--force`  Reprocesa aunque existan salidas

**Notas:**
- Si **no** especificas `-q`, el script usa **0** (solo adaptadores).
- Si deseas recorte por calidad, usa por ejemplo `-q 20` o `-q 30`.

---

## Ejemplos

**Solo cortar adaptadores (sin trimming de calidad):**
```bash
./trim_skewer.sh -i ./fastq -o ./00_skewer -a adapt.fa -t 32 -q 0
```

**Adaptadores + calidad ≥ 20 (recomendado en general):**
```bash
./trim_skewer.sh -i ./fastq -o ./00_skewer -a adapt.fa -t 32 -q 20
```

**Ejecutar en segundo plano con `nohup`:**
```bash
nohup ./trim_skewer.sh -i ./fastq -o ./00_skewer -a adapt.fa -t 60 -q 20 --force > trim_skewer.out 2>&1 &
tail -f trim_skewer.out
```

---

## Salidas

Para cada muestra detectada (por ejemplo `SAMPLE`):

```
OUTPUT_DIR/
├── logs/
│   └── SAMPLE.log
└── SAMPLE/
    ├── SAMPLE-trimmed-pair1.fastq.gz
    ├── SAMPLE-trimmed-pair2.fastq.gz
    ├── SAMPLE-trimmed-unpaired1.fastq.gz   # sólo si Skewer las genera
    └── SAMPLE-trimmed-unpaired2.fastq.gz   # sólo si Skewer las genera
```

Además, se genera/actualiza:

```
OUTPUT_DIR/manifest.tsv
```

**`manifest.tsv`** columnas:
```
sample  r1_in   r2_in   r1_out  r2_out  unpaired1  unpaired2  log
```

Este manifiesto facilita encadenar el siguiente paso (e.g., BBDuk, Bowtie2, FastQC…).

---

## Archivo de adaptadores (`adapt.fa`)

Incluye adaptadores universales de **Illumina TruSeq** y **Nextera**, además de colas homopoliméricas opcionales (útiles en plataformas 2‑color que generan colas de G).

```fasta
>Illumina_TruSeq_Adapter_Read1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_TruSeq_Adapter_Read2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# (Opcional) Nextera / tagmentación
>Nextera_Adapter_Read1
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>Nextera_Adapter_Read2
CTGTCTCTTATACACATCTGACGCTGCCGACGA

# (Opcional) Nextera mosaic end (ME) mínima
>Nextera_ME_min
CTGTCTCTTATACACATCT

# (Opcional) Poly-G tail (NextSeq/Novaseq 2-color)
>Poly_G_30
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

# (Opcional) Poly-A/T
>Poly_A_30
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>Poly_T_30
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

> Si tu preparación fue **small RNA** u otra con adaptadores específicos, agrega las secuencias correspondientes aquí.

---

## Script (`trim_skewer.sh`)

```bash
#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
ADAPTERS=""
THREADS=8
QUAL=0
FORCE=false

print_usage() {
  cat <<EOF
Uso: $(basename "$0") -i INPUT_DIR -o OUTPUT_DIR -a ADAPTERS.fa [-t THREADS] [-q QUAL] [--force]

Opciones:
  -i   Directorio con FASTQ(.gz) de entrada
  -o   Directorio de salida
  -a   Archivo FASTA con adaptadores (Skewer -x)
  -t   Hilos (por defecto: 8)
  -q   Calidad mínima Phred para trimming (0 = sin recorte, por defecto: 0)
  --force  Reprocesar aunque existan salidas
EOF
}

# Parseo de --force antes de getopts
long_opts=()
for arg in "$@"; do
  case "${arg}" in
    --force) FORCE=true ;;
    *) long_opts+=("$arg") ;;
  esac
done
set -- "${long_opts[@]}"

# Flags cortos
while getopts ":i:o:a:t:q:h" opt; do
  case $opt in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    a) ADAPTERS="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    q) QUAL="$OPTARG" ;;
    h) print_usage; exit 0 ;;
    \?) echo "Opción inválida: -$OPTARG" >&2; print_usage; exit 1 ;;
    :) echo "La opción -$OPTARG requiere un argumento." >&2; print_usage; exit 1 ;;
  esac
done || true

# Validaciones
[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$ADAPTERS" ]] && { print_usage; exit 1; }
[[ ! -d "$INPUT_DIR" ]] && { echo "No existe INPUT_DIR: $INPUT_DIR" >&2; exit 1; }
[[ ! -f "$ADAPTERS" ]] && { echo "No existe archivo de adaptadores: $ADAPTERS" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR/logs"
MANIFEST="$OUTPUT_DIR/manifest.tsv"
if [[ ! -f "$MANIFEST" ]]; then
  echo -e "sample	r1_in	r2_in	r1_out	r2_out	unpaired1	unpaired2	log" > "$MANIFEST"
fi

# Intenta derivar R2 desde R1 para varios patrones (incluye *_R1_001.fastq.gz)
derive_r2() {
  local r1="$1"
  local r2=""
  for pat in "_1.fq.gz:_2.fq.gz"              "_1.fastq.gz:_2.fastq.gz"              "_R1.fq.gz:_R2.fq.gz"              "_R1.fastq.gz:_R2.fastq.gz"              "_R1_001.fastq.gz:_R2_001.fastq.gz"; do
    IFS=":" read -r a b <<< "$pat"
    if [[ "$r1" == *"$a" ]]; then
      r2="${r1/$a/$b}"
      [[ -f "$r2" ]] && { echo "$r2"; return 0; }
    fi
  done
  echo ""
  return 1
}

# Nombre de muestra sin sufijos R1/R2/_001
sample_name() {
  local path="$1"
  local base
  base="$(basename "$path")"
  base="${base%.fastq.gz}"
  base="${base%.fq.gz}"
  base="${base/_R1_001/}"
  base="${base/_R2_001/}"
  base="${base/_R1/}"
  base="${base/_R2/}"
  base="${base/_1/}"
  base="${base/_2/}"
  base="${base/R1/}"
  base="${base/R2/}"
  echo "$base"
}

# Descubrir R1
shopt -s nullglob
mapfile -t R1_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f   \( -name "*_1.fq.gz" -o -name "*_1.fastq.gz"      -o -name "*_R1.fq.gz" -o -name "*_R1.fastq.gz"      -o -name "*_R1_001.fastq.gz" \)   | sort)

if (( ${#R1_FILES[@]} == 0 )); then
  echo "No se encontraron archivos R1 en: $INPUT_DIR" >&2
  exit 1
fi

# Procesamiento
for R1 in "${R1_FILES[@]}"; do
  R2="$(derive_r2 "$R1" || true)"
  if [[ -z "$R2" || ! -f "$R2" ]]; then
    echo "✗ No se encontró el R2 correspondiente a: $R1. Saltando."
    continue
  fi

  SAMPLE="$(sample_name "$R1")"
  OUT_DIR="$OUTPUT_DIR/$SAMPLE"
  mkdir -p "$OUT_DIR"

  OUT_PREFIX="$OUT_DIR/$SAMPLE"
  LOG_FILE="$OUTPUT_DIR/logs/${SAMPLE}.log"

  R1_OUT="${OUT_PREFIX}-trimmed-pair1.fastq.gz"
  R2_OUT="${OUT_PREFIX}-trimmed-pair2.fastq.gz"
  U1_OUT="${OUT_PREFIX}-trimmed-unpaired1.fastq.gz"
  U2_OUT="${OUT_PREFIX}-trimmed-unpaired2.fastq.gz"

  if [[ "$FORCE" = false && -s "$R1_OUT" && -s "$R2_OUT" ]]; then
    echo "• Ya existen salidas para $SAMPLE. Usa --force para reprocesar. Saltando."
    continue
  fi

  echo "→ Procesando $SAMPLE"
  {
    echo "== $(date) =="
    echo "R1: $R1"
    echo "R2: $R2"
    echo "Adaptadores: $ADAPTERS"
    echo "Hilos: $THREADS"
    echo "Calidad mínima (q): $QUAL"
    echo "Salida: $OUT_PREFIX"
  } >> "$LOG_FILE"

  set +e
  skewer -x "$ADAPTERS" -t "$THREADS" -q "$QUAL" -m pe -z -o "$OUT_PREFIX" "$R1" "$R2" >> "$LOG_FILE" 2>&1
  status=$?
  set -e

  if [[ $status -ne 0 ]]; then
    echo "✗ Skewer falló para $SAMPLE (ver $LOG_FILE)."
    continue
  fi

  if [[ ! -s "$R1_OUT" || ! -s "$R2_OUT" ]]; then
    echo "✗ No se generaron archivos pareados para $SAMPLE (ver $LOG_FILE)."
    continue
  fi

  echo -e "${SAMPLE}\t${R1}\t${R2}\t${R1_OUT}\t${R2_OUT}\t${U1_OUT:-}\t${U2_OUT:-}\t${LOG_FILE}" >> "$MANIFEST"

  echo "✓ Listo: $SAMPLE"
done

echo "== Terminado =="
echo "Manifiesto: $MANIFEST"
```

---

## FAQ rápida

- **¿Skewer corta por calidad por defecto?**  
  No. Solo si pasas `-q <score>`. En este script, el **default es 0** (sin recorte de calidad) para evitar sobre‑filtrado.

- **¿Qué valor de `-q` usar?**  
  - `-q 0`: solo adaptadores  
  - `-q 20`: recorte moderado  
  - `-q 30`: recorte estricto (Phred 30)

- **¿Cómo verifico el progreso en segundo plano?**  
  Lanza con `nohup ... > run.out 2>&1 &` y mira con `tail -f run.out`.  
  También puedes inspeccionar `OUTPUT_DIR/logs/<sample>.log`.
