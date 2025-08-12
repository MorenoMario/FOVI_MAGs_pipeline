# Skewer PE Trimming – Illumina (TruSeq/Nextera)

Robust script for **adapter trimming** and **(optional) quality trimming** on Illumina **paired-end** reads.  
Automatically detects file names like `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` (BCL2FASTQ output) and other common formats.

- Produces **compressed outputs** (`.fastq.gz`)
- Creates **per-sample logs** and a **manifest.tsv** for chaining downstream steps
- Resumable: skips samples already processed (unless you use `--force`)

---

## Requirements

- Linux/macOS with `bash`
- [`skewer`](https://github.com/relipmoc/skewer) ≥ 0.2.2 in the `PATH`
- `gzip`, `coreutils`

Check:
```bash
skewer --version
```

---

## Installation

Place these files in your repo/project:

```
.
├── trim_skewer.sh
└── adapt.fa
```

Make the script executable:

```bash
chmod +x trim_skewer.sh
```

---

## Usage

```bash
./trim_skewer.sh -i INPUT_DIR -o OUTPUT_DIR -a ADAPTERS.fa [-t THREADS] [-q QUAL] [--force]
```

**Options:**
- `-i`  Directory with FASTQ(.gz) input files  
- `-o`  Output directory  
- `-a`  Adapter FASTA file (for `skewer -x`)  
- `-t`  Threads (default: 8)  
- `-q`  **Minimum Phred quality** for trimming (default **0** = no quality trimming)  
- `--force`  Reprocess even if output files already exist

**Notes:**
- If you **don’t** specify `-q`, the script uses **0** (adapter trimming only).
- For quality trimming, use e.g., `-q 20` or `-q 30`.

---

## Examples

**Adapter trimming only (no quality trimming):**
```bash
./trim_skewer.sh -i ./fastq -o ./00_skewer -a adapt.fa -t 32 -q 0
```

**Adapters + quality ≥ 20 (recommended in most cases):**
```bash
./trim_skewer.sh -i ./fastq -o ./00_skewer -a adapt.fa -t 32 -q 20
```

**Run in background with `nohup`:**
```bash
nohup ./trim_skewer.sh -i ./fastq -o ./00_skewer -a adapt.fa -t 60 -q 20 --force > trim_skewer.out 2>&1 &
tail -f trim_skewer.out
```

---

## Output

For each detected sample (e.g., `SAMPLE`):

```
OUTPUT_DIR/
├── logs/
│   └── SAMPLE.log
└── SAMPLE/
    ├── SAMPLE-trimmed-pair1.fastq.gz
    ├── SAMPLE-trimmed-pair2.fastq.gz
    ├── SAMPLE-trimmed-unpaired1.fastq.gz   # only if Skewer produces them
    └── SAMPLE-trimmed-unpaired2.fastq.gz   # only if Skewer produces them
```

Additionally, a manifest is generated/updated:

```
OUTPUT_DIR/manifest.tsv
```

**`manifest.tsv`** columns:
```
sample  r1_in   r2_in   r1_out  r2_out  unpaired1  unpaired2  log
```

This manifest makes it easy to chain the next steps (e.g., BBDuk, Bowtie2, FastQC…).

---

## Adapter file (`adapt.fa`)

Includes universal adapters for **Illumina TruSeq** and **Nextera**, plus optional homopolymer tails (useful in 2‑color platforms generating G-tails).

```fasta
>Illumina_TruSeq_Adapter_Read1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_TruSeq_Adapter_Read2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# (Optional) Nextera / tagmentation
>Nextera_Adapter_Read1
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>Nextera_Adapter_Read2
CTGTCTCTTATACACATCTGACGCTGCCGACGA

# (Optional) Nextera mosaic end (ME) minimal
>Nextera_ME_min
CTGTCTCTTATACACATCT

# (Optional) Poly-G tail (NextSeq/Novaseq 2-color)
>Poly_G_30
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

# (Optional) Poly-A/T
>Poly_A_30
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>Poly_T_30
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

> If your library prep was **small RNA** or other with specific adapters, add them here.

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
Usage: $(basename "$0") -i INPUT_DIR -o OUTPUT_DIR -a ADAPTERS.fa [-t THREADS] [-q QUAL] [--force]

Options:
  -i   Directory with FASTQ(.gz) input
  -o   Output directory
  -a   Adapter FASTA file (Skewer -x)
  -t   Threads (default: 8)
  -q   Minimum Phred quality for trimming (0 = no quality trimming, default: 0)
  --force  Reprocess even if outputs exist
EOF
}

# Parse --force before getopts
long_opts=()
for arg in "$@"; do
  case "${arg}" in
    --force) FORCE=true ;;
    *) long_opts+=("$arg") ;;
  esac
done
set -- "${long_opts[@]}"

# Short flags
while getopts ":i:o:a:t:q:h" opt; do
  case $opt in
    i) INPUT_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    a) ADAPTERS="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    q) QUAL="$OPTARG" ;;
    h) print_usage; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; print_usage; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; print_usage; exit 1 ;;
  esac
done || true

# Validations
[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$ADAPTERS" ]] && { print_usage; exit 1; }
[[ ! -d "$INPUT_DIR" ]] && { echo "INPUT_DIR does not exist: $INPUT_DIR" >&2; exit 1; }
[[ ! -f "$ADAPTERS" ]] && { echo "Adapter file does not exist: $ADAPTERS" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR/logs"
MANIFEST="$OUTPUT_DIR/manifest.tsv"
if [[ ! -f "$MANIFEST" ]]; then
  echo -e "sample	r1_in	r2_in	r1_out	r2_out	unpaired1	unpaired2	log" > "$MANIFEST"
fi

# Try to derive R2 from R1 (supports *_R1_001.fastq.gz)
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

# Clean sample name
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

# Discover R1
shopt -s nullglob
mapfile -t R1_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f   \( -name "*_1.fq.gz" -o -name "*_1.fastq.gz"      -o -name "*_R1.fq.gz" -o -name "*_R1.fastq.gz"      -o -name "*_R1_001.fastq.gz" \)   | sort)

if (( ${#R1_FILES[@]} == 0 )); then
  echo "No R1 files found in: $INPUT_DIR" >&2
  exit 1
fi

# Processing
for R1 in "${R1_FILES[@]}"; do
  R2="$(derive_r2 "$R1" || true)"
  if [[ -z "$R2" || ! -f "$R2" ]]; then
    echo "✗ Could not find matching R2 for: $R1. Skipping."
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
    echo "• Outputs already exist for $SAMPLE. Use --force to reprocess. Skipping."
    continue
  fi

  echo "→ Processing $SAMPLE"
  {
    echo "== $(date) =="
    echo "R1: $R1"
    echo "R2: $R2"
    echo "Adapters: $ADAPTERS"
    echo "Threads: $THREADS"
    echo "Quality threshold (q): $QUAL"
    echo "Output prefix: $OUT_PREFIX"
  } >> "$LOG_FILE"

  set +e
  skewer -x "$ADAPTERS" -t "$THREADS" -q "$QUAL" -m pe -z -o "$OUT_PREFIX" "$R1" "$R2" >> "$LOG_FILE" 2>&1
  status=$?
  set -e

  if [[ $status -ne 0 ]]; then
    echo "✗ Skewer failed for $SAMPLE (see $LOG_FILE)."
    continue
  fi

  if [[ ! -s "$R1_OUT" || ! -s "$R2_OUT" ]]; then
    echo "✗ No paired outputs generated for $SAMPLE (see $LOG_FILE)."
    continue
  fi

  echo -e "${SAMPLE}	${R1}	${R2}	${R1_OUT}	${R2_OUT}	${U1_OUT:-}	${U2_OUT:-}	${LOG_FILE}" >> "$MANIFEST"

  echo "✓ Done: $SAMPLE"
done

echo "== Finished =="
echo "Manifest: $MANIFEST"
```

---

## Quick FAQ

- **Does Skewer trim by quality by default?**  
  No. Only if you pass `-q <score>`. In this script, the **default is 0** (no quality trimming).

- **Which `-q` value should I use?**  
  - `-q 0`: adapters only  
  - `-q 20`: moderate trimming  
  - `-q 30`: strict trimming (Phred 30)

- **How to monitor progress in background?**  
  Run with `nohup ... > run.out 2>&1 &` and check with `tail -f run.out`.  
  You can also inspect `OUTPUT_DIR/logs/<sample>.log`.
