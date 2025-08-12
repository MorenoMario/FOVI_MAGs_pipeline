# Read filtering against *Macrocystis pyrifera* reference (BBDuk)

This script filters **paired-end** reads against the *Macrocystis pyrifera* reference genome using `bbduk.sh`.  
It supports two modes:

1. **Keep both sets**: reads matching the reference (hits) and reads not matching (no-hits).
2. **Keep only no-hits**: discard reads matching the reference.

---

## Requirements

- Linux/macOS with `bash`
- [BBTools / BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) in the `PATH`
- `gzip`

Check:
```bash
bbduk.sh --version
```

---

## Usage examples

### 1️⃣ Keep both sets (hits and no-hits)
```bash
./filter_host_bbduk.sh   -i ./00_skewer/   -o ./01_bbduk_ref/   -r ./GCA_031763025.fna   -t 32   -m 150g
```

This will generate:
- `<sample>_macroHit_R1.fastq.gz` / `_R2.fastq.gz` (hits)
- `<sample>_non-macroHit_R1.fastq.gz` / `_R2.fastq.gz` (no-hits)
- `logs/<sample>.bbduk.stats.txt`
- `manifest_bbduk.tsv`

### 2️⃣ Keep only no-hits
Edit the script to set:
```
outm=/dev/null
outm2=/dev/null
```
or create a copy of the script with this configuration.

Example:
```bash
./filter_host_bbduk.sh   -i ./00_skewer/   -o ./01_bbduk_ref/   -r ./GCA_031763025.fna   -t 32   -m 150g   --force
```

### 3️⃣ Run in background with `nohup`
```bash
nohup ./filter_host_bbduk.sh   -i ./00_skewer/   -o ./01_bbduk_ref/   -r ./GCA_031763025.fna   -t 60   -m 150g   --force   > filter_bbduk.out 2>&1 &
tail -f filter_bbduk.out
```

---

## Inputs in different folders

If your `*-trimmed-pair1.fastq.gz` and `*-trimmed-pair2.fastq.gz` files are inside subfolders, e.g.:

```
/media/mmorenos/disk3/DSG123-selected/rawReads/00_skewer/20250604_PamelaFernandez_DSG_H1_S7/
```

Make sure the script searches recursively. In the R1 search section, change:

```bash
mapfile -t R1_FILES < <(find "$INPUT_DIR" -type f \
  \( -name "*-trimmed-pair1.fastq.gz" -o -name "*-trimmed-pair1.fastq" \) | sort)
```

This way, it will process files located in any subfolder of `INPUT_DIR`.

