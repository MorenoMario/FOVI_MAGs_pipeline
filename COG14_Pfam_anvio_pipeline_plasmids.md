# Anotación COG14 + Pfam (AA de anvi’o) e importación a contigs.db

Este documento describe un flujo reproducible para:

1. Exportar **proteínas (AA)** desde múltiples `*.db` de anvi’o.
2. Anotar AA contra **COG14** usando **DIAMOND** (salida `m6`).
3. Convertir `m6` → **tabla de funciones** (formato `anvi-import-functions`) usando los archivos RAW de COG14.
4. Convertir a **UTF-8**.
5. Importar funciones COG14 en cada `*.db`.
6. Anotar AA contra **Pfam v32** usando **hmmscan** (HMMER).
7. Convertir `domtblout` → tabla de funciones e importar a cada `*.db`.

> Importante:
> - Este pipeline **no activa ni desactiva ambientes**. Asume que tú ya activaste el entorno correcto (`anvio-8`, `hmmer`, `diamond`, etc.).
> - Ajusta `THREADS` según tus recursos.
> - La clave es que el **header** del FASTA exportado desde anvi’o sea el `gene_callers_id` (numérico); eso permite que `anvi-import-functions` asocie correctamente.

---

## 0) Rutas de trabajo (tus rutas actuales)

### DBs de anvi’o
- `DB_DIR=/media/mmorenos/disk3/010_fovi_metagenome/db_anvio8/cathy_db`
  - contiene: `chon.db  H_S1_6.db  ilq2.db  Y_S256.db`

### FASTA de proteínas (AA) exportados
- `FAA_DIR=/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa`
  - ejemplo: `chon.faa  H_S1_6.faa  ilq2.faa  Y_S256.faa`

### COG14 (DIAMOND + RAW)
- `COG_DMND=/media/mmorenos/disk3/012_eduCastroLab/00_old_an/COG_2014/COG14/DB_DIAMOND/COG.dmnd`
- `COG_RAW_DIR=/media/mmorenos/disk3/012_eduCastroLab/00_old_an/COG_2014/COG14/RAW_DATA_FROM_NCBI`

### Pfam v32
- `PFAM_DIR=/media/mmorenos/disk3/012_eduCastroLab/00_old_an/Pfam_v32`
- `PFAM_HMM=${PFAM_DIR}/Pfam-A.hmm`

---

## 1) Exportar AA desde TODOS los DBs (anvi-get-sequences-for-gene-calls)

Guarda como: `01_export_faa_from_dbs.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

DB_DIR="/media/mmorenos/disk3/010_fovi_metagenome/db_anvio8/cathy_db"
FAA_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa"

mkdir -p "${FAA_DIR}"

for db in "${DB_DIR}"/*.db; do
  [[ -e "$db" ]] || { echo "[ERROR] No hay *.db en ${DB_DIR}" >&2; exit 1; }
  sample="$(basename "$db" .db)"
  out="${FAA_DIR}/${sample}.faa"

  echo "============================================="
  echo "[INFO] Sample: ${sample}"
  echo "[INFO] DB    : ${db}"
  echo "[INFO] FAA   : ${out}"
  echo "============================================="

  anvi-get-sequences-for-gene-calls -c "$db" --get-aa-sequences -o "$out"

  # sanity check
  if [[ ! -s "$out" ]]; then
    echo "[ERROR] FAA vacío: $out" >&2
    exit 1
  fi
done

echo "[OK] Exportación de AA finalizada."
```

---

## 2) COG14 con DIAMOND (genera m6 por muestra)

Guarda como: `02_diamond_cog14.sh`

> Nota: tu error anterior fue `Invalid output field: mismatches`.
> En DIAMOND el campo válido es **`mismatch`** (singular). Este script ya lo corrige.

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

FAA_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa"
COG_DMND="/media/mmorenos/disk3/012_eduCastroLab/00_old_an/COG_2014/COG14/DB_DIAMOND/COG.dmnd"

OUT_DIR="${FAA_DIR}/diamond_cog14_out"
TMP_BASE="${OUT_DIR}/tmp"
LOG_DIR="${OUT_DIR}/logs"

THREADS=10

mkdir -p "${OUT_DIR}" "${TMP_BASE}" "${LOG_DIR}"

for faa in "${FAA_DIR}"/*.faa; do
  sample="$(basename "${faa}" .faa)"
  m6="${OUT_DIR}/${sample}.cog14.m6.tsv"
  log="${LOG_DIR}/${sample}.diamond_cog14.log"

  echo "============================================="
  echo "[INFO] Sample   : ${sample}"
  echo "[INFO] FAA      : ${faa}"
  echo "[INFO] DB (dmnd): ${COG_DMND}"
  echo "[INFO] Out m6   : ${m6}"
  echo "[INFO] Threads  : ${THREADS}"
  echo "============================================="

  mkdir -p "${TMP_BASE}/${sample}"

  diamond blastp     -q "${faa}"     -d "${COG_DMND}"     -o "${m6}"     -p "${THREADS}"     --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore     --max-target-seqs 1     --evalue 1e-5     -t "${TMP_BASE}/${sample}"     >"${log}" 2>&1

  [[ -s "${m6}" ]] || { echo "[ERROR] m6 vacío: ${m6}" >&2; exit 1; }
done

echo "[OK] DIAMOND COG14 finalizado."
```

---

## 3) Convertir m6 → funciones (formato anvi-import-functions)

Guarda como: `03_make_COG14_functions_from_m6.sh`

- Entrada: `diamond_cog14_out/*.cog14.m6.tsv`
- Salida: `cog14_functions_out/<sample>.cog14_functions.tsv`

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

M6_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/diamond_cog14_out"
OUT_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/cog14_functions_out"

RAW_DIR="/media/mmorenos/disk3/012_eduCastroLab/00_old_an/COG_2014/COG14/RAW_DATA_FROM_NCBI"
COG_NAMES="${RAW_DIR}/cognames2003-2014.tab"
COG_CSV="${RAW_DIR}/cog2003-2014.csv"

GI_COG_DESC="${OUT_DIR}/GI_TO_COG14_DESC.tsv"

mkdir -p "${OUT_DIR}"

echo "============================================="
echo "[INFO] M6_DIR      : ${M6_DIR}"
echo "[INFO] OUT_DIR     : ${OUT_DIR}"
echo "[INFO] RAW_DIR     : ${RAW_DIR}"
echo "[INFO] COG_NAMES   : ${COG_NAMES}"
echo "[INFO] COG_CSV     : ${COG_CSV}"
echo "[INFO] GI_COG_DESC : ${GI_COG_DESC}"
echo "============================================="

for f in "${COG_NAMES}" "${COG_CSV}"; do
  [[ -f "$f" ]] || { echo "[ERROR] No se encuentra: $f" >&2; exit 1; }
done

# 1) Cache GI -> (COG_ID, desc)
if [[ ! -s "${GI_COG_DESC}" ]]; then
  echo "[INFO] Generando GI->COG_ID+desc (cache)..."
  awk '
  BEGIN{ OFS="\t"; }
  FNR==1 && FILENAME ~ /cognames/ { FS="\t"; }
  FNR==1 && FILENAME ~ /cog2003-2014.csv/ { FS=","; }

  FILENAME ~ /cognames2003-2014.tab/ {
      cog_id=$1; desc=$3; cog_desc[cog_id]=desc; next;
  }
  FILENAME ~ /cog2003-2014.csv/ {
      gi=$1; cog_id=$7;
      if (cog_id != "" && (cog_id in cog_desc)) {
          print gi, cog_id, cog_desc[cog_id];
      }
      next;
  }' "${COG_NAMES}" "${COG_CSV}" > "${GI_COG_DESC}"

  echo "[OK] GI_COG_DESC generado: ${GI_COG_DESC}"
  head -n 3 "${GI_COG_DESC}" || true
else
  echo "[INFO] Usando cache existente: ${GI_COG_DESC}"
fi

# 2) Para cada m6: generar funciones
m6s=( "${M6_DIR}"/*.cog14.m6.tsv )
(( ${#m6s[@]} > 0 )) || { echo "[ERROR] No encontré m6 en ${M6_DIR}" >&2; exit 1; }

for M6_FILE in "${m6s[@]}"; do
  sample="$(basename "${M6_FILE}" .cog14.m6.tsv)"
  FUN_FILE="${OUT_DIR}/${sample}.cog14_functions.tsv"

  echo "---------------------------------------------"
  echo "[INFO] Sample  : ${sample}"
  echo "[INFO] M6_FILE : ${M6_FILE}"
  echo "[INFO] FUN_FILE: ${FUN_FILE}"

  {
    printf "gene_callers_id\tsource\taccession\tfunction\te_value\n";

    awk -v OFS="\t" '
    BEGIN{ FS="\t"; }
    NR==FNR {
      gi=$1; cog_id=$2; desc=$3;
      gi_map[gi]=cog_id "\t" desc;
      next;
    }
    {
      n=split($0, f, /[ \t]+/);
      if (n < 12) next;

      qseqid=f[1];
      sseqid=f[2];
      evalue=f[11];

      gi="";
      if (match(sseqid, /^gi\|([0-9]+)\|/, m)) {
        gi=m[1];
      }

      if (gi != "" && (gi in gi_map)) {
        split(gi_map[gi], tmp, "\t");
        cog_id=tmp[1];
        desc=tmp[2];
        print qseqid, "COG14_FUNCTION", cog_id, desc, evalue;
      }
    }' "${GI_COG_DESC}" "${M6_FILE}"
  } > "${FUN_FILE}"

  [[ -s "${FUN_FILE}" ]] || { echo "[ERROR] FUN_FILE vacío: ${FUN_FILE}" >&2; exit 1; }
done

echo "[OK] Funciones COG14 generadas en: ${OUT_DIR}"
```

---

## 4) Convertir funciones a UTF-8 (por muestra)

Guarda como: `04_utf8_cog14_functions.sh`

- Entrada: `cog14_functions_out/*.cog14_functions.tsv`
- Salida: `cog14_functions_out_utf8/<sample>.cog14_functions.utf8.tsv`

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

IN_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/cog14_functions_out"
OUT_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/cog14_functions_out_utf8"

mkdir -p "${OUT_DIR}"

for FUN_FILE in "${IN_DIR}"/*.cog14_functions.tsv; do
  sample="$(basename "${FUN_FILE}" .cog14_functions.tsv)"
  FUN_FILE_UTF8="${OUT_DIR}/${sample}.cog14_functions.utf8.tsv"

  echo "============================================="
  echo "[INFO] Sample : ${sample}"
  echo "[INFO] Input  : ${FUN_FILE}"
  echo "[INFO] Output : ${FUN_FILE_UTF8}"
  echo "============================================="

  if command -v iconv >/dev/null 2>&1; then
    iconv -f WINDOWS-1252 -t UTF-8 -c "${FUN_FILE}" > "${FUN_FILE_UTF8}"
  else
    python - << EOF
inp  = "${FUN_FILE}"
outp = "${FUN_FILE_UTF8}"
with open(inp, "rb") as f_in, open(outp, "w", encoding="utf-8") as f_out:
    for bline in f_in:
        f_out.write(bline.decode("latin-1", errors="ignore"))
EOF
  fi

  [[ -s "${FUN_FILE_UTF8}" ]] || { echo "[ERROR] UTF8 vacío: ${FUN_FILE_UTF8}" >&2; exit 1; }
done

echo "[OK] UTF-8 generado en: ${OUT_DIR}"
```

---

## 5) Importar COG14 en TODOS los DBs (anvi-import-functions)

Guarda como: `05_import_cog14_to_anvio.sh`

- Entrada: `cog14_functions_out_utf8/<sample>.cog14_functions.utf8.tsv`
- DBs: `DB_DIR/*.db`

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

DB_DIR="/media/mmorenos/disk3/010_fovi_metagenome/db_anvio8/cathy_db"
FUN_UTF8_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/cog14_functions_out_utf8"

LOG_DIR="${FUN_UTF8_DIR}/logs_import_anvio"
mkdir -p "${LOG_DIR}"

echo "============================================="
echo "[INFO] Importando COG14 en anvi'o (por muestra)"
echo "[INFO] DB_DIR      : ${DB_DIR}"
echo "[INFO] FUN_UTF8_DIR: ${FUN_UTF8_DIR}"
echo "[INFO] LOG_DIR     : ${LOG_DIR}"
echo "============================================="

for CONTIGS_DB in "${DB_DIR}"/*.db; do
  [[ -e "$CONTIGS_DB" ]] || { echo "[ERROR] No hay *.db en ${DB_DIR}" >&2; exit 1; }
  sample="$(basename "${CONTIGS_DB}" .db)"

  FUN_FILE_UTF8="${FUN_UTF8_DIR}/${sample}.cog14_functions.utf8.tsv"
  LOG_OUT="${LOG_DIR}/${sample}.anvi_import_cog14.out"
  LOG_ERR="${LOG_DIR}/${sample}.anvi_import_cog14.err"

  echo "-------------------------------------------------"
  echo "[INFO] Sample     : ${sample}"
  echo "[INFO] CONTIGS_DB : ${CONTIGS_DB}"
  echo "[INFO] FUN_FILE   : ${FUN_FILE_UTF8}"

  if [[ ! -s "${FUN_FILE_UTF8}" ]]; then
    echo "[SKIP] No existe o está vacío: ${FUN_FILE_UTF8}" | tee -a "${LOG_OUT}"
    continue
  fi

  anvi-import-functions -c "${CONTIGS_DB}" -i "${FUN_FILE_UTF8}"     >"${LOG_OUT}" 2>"${LOG_ERR}"

  echo "[OK] Importado: ${sample}"
done

echo "[OK] Importación COG14 completada."
```

---

## 6) Pfam v32 con hmmscan (genera domtblout por muestra)

Guarda como: `06_hmmscan_pfam.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

FAA_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa"
PFAM_HMM="/media/mmorenos/disk3/012_eduCastroLab/00_old_an/Pfam_v32/Pfam-A.hmm"

OUT_DIR="${FAA_DIR}/hmmscan_pfam_out"
LOG_DIR="${OUT_DIR}/logs"
THREADS=10

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

[[ -f "${PFAM_HMM}" ]] || { echo "[ERROR] No existe PFAM_HMM: ${PFAM_HMM}" >&2; exit 1; }

for faa in "${FAA_DIR}"/*.faa; do
  sample="$(basename "${faa}" .faa)"
  domtbl="${OUT_DIR}/${sample}.pfam.domtblout"
  log="${LOG_DIR}/${sample}.pfam.hmmscan.log"

  echo "============================================="
  echo "[INFO] Sample : ${sample}"
  echo "[INFO] FAA    : ${faa}"
  echo "[INFO] HMM    : ${PFAM_HMM}"
  echo "[INFO] domtbl : ${domtbl}"
  echo "[INFO] Threads: ${THREADS}"
  echo "============================================="

  hmmscan --cpu "${THREADS}"     --domtblout "${domtbl}"     "${PFAM_HMM}" "${faa}"     >"${log}" 2>&1

  [[ -s "${domtbl}" ]] || { echo "[ERROR] domtbl vacío: ${domtbl}" >&2; exit 1; }
done

echo "[OK] hmmscan Pfam finalizado."
```

---

## 7) Convertir domtblout → funciones (formato anvi-import-functions)

Guarda como: `07_parse_pfam_domtblout_to_functions.sh`

- Entrada: `hmmscan_pfam_out/*.pfam.domtblout`
- Salida: `pfam_functions_out_utf8/<sample>.pfam_functions.utf8.tsv`

> Nota: Este parser toma el **mejor hit por query** según menor `i-Evalue` (columna 13 en domtblout).
> Usa:
> - `accession` = accession del Pfam (ej. `PF00005.27`)
> - `function`  = nombre del target (ej. `ABC_tran`)
> - `source`    = `PFAM32`

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

DOMTBL_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/hmmscan_pfam_out"
OUT_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/pfam_functions_out_utf8"

mkdir -p "${OUT_DIR}"

doms=( "${DOMTBL_DIR}"/*.pfam.domtblout )
(( ${#doms[@]} > 0 )) || { echo "[ERROR] No encontré domtblout en ${DOMTBL_DIR}" >&2; exit 1; }

for domtbl in "${doms[@]}"; do
  sample="$(basename "${domtbl}" .pfam.domtblout)"
  out="${OUT_DIR}/${sample}.pfam_functions.utf8.tsv"

  echo "============================================="
  echo "[INFO] Sample : ${sample}"
  echo "[INFO] Input  : ${domtbl}"
  echo "[INFO] Output : ${out}"
  echo "============================================="

  awk '
  BEGIN{
    OFS="\t";
    print "gene_callers_id\tsource\taccession\tfunction\te_value";
  }
  $0 ~ /^#/ { next; }
  {
    target_name=$1;
    target_acc=$2;
    query=$4;
    ievalue=$13;

    if (!(query in best) || ievalue+0 < best_e[query]+0) {
      best[query]=target_name "\t" target_acc "\t" ievalue;
      best_e[query]=ievalue;
    }
  }
  END{
    for (q in best) {
      split(best[q], a, "\t");
      print q, "PFAM32", a[2], a[1], a[3];
    }
  }' "${domtbl}"   | sort -k1,1n -k5,5g > "${out}"

  [[ -s "${out}" ]] || { echo "[ERROR] Output vacío: ${out}" >&2; exit 1; }
done

echo "[OK] Funciones Pfam (UTF-8) generadas en: ${OUT_DIR}"
```

---

## 8) Importar Pfam en TODOS los DBs (anvi-import-functions)

Guarda como: `08_import_pfam_to_anvio.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

DB_DIR="/media/mmorenos/disk3/010_fovi_metagenome/db_anvio8/cathy_db"
PFAM_FUN_DIR="/media/mmorenos/disk3/010_fovi_metagenome/plasmidomics/cathy_faa/pfam_functions_out_utf8"

LOG_DIR="${PFAM_FUN_DIR}/logs_import_anvio"
mkdir -p "${LOG_DIR}"

echo "============================================="
echo "[INFO] Importando Pfam en anvi'o (por muestra)"
echo "[INFO] DB_DIR      : ${DB_DIR}"
echo "[INFO] PFAM_FUN_DIR: ${PFAM_FUN_DIR}"
echo "[INFO] LOG_DIR     : ${LOG_DIR}"
echo "============================================="

for CONTIGS_DB in "${DB_DIR}"/*.db; do
  [[ -e "$CONTIGS_DB" ]] || { echo "[ERROR] No hay *.db en ${DB_DIR}" >&2; exit 1; }
  sample="$(basename "${CONTIGS_DB}" .db)"

  FUN_FILE_UTF8="${PFAM_FUN_DIR}/${sample}.pfam_functions.utf8.tsv"
  LOG_OUT="${LOG_DIR}/${sample}.anvi_import_pfam.out"
  LOG_ERR="${LOG_DIR}/${sample}.anvi_import_pfam.err"

  echo "-------------------------------------------------"
  echo "[INFO] Sample     : ${sample}"
  echo "[INFO] CONTIGS_DB : ${CONTIGS_DB}"
  echo "[INFO] FUN_FILE   : ${FUN_FILE_UTF8}"

  if [[ ! -s "${FUN_FILE_UTF8}" ]]; then
    echo "[SKIP] No existe o está vacío: ${FUN_FILE_UTF8}" | tee -a "${LOG_OUT}"
    continue
  fi

  anvi-import-functions -c "${CONTIGS_DB}" -i "${FUN_FILE_UTF8}"     >"${LOG_OUT}" 2>"${LOG_ERR}"

  echo "[OK] Importado: ${sample}"
done

echo "[OK] Importación Pfam completada."
```

---

## 9) Orden sugerido de ejecución

```bash
bash 01_export_faa_from_dbs.sh
bash 02_diamond_cog14.sh
bash 03_make_COG14_functions_from_m6.sh
bash 04_utf8_cog14_functions.sh
bash 05_import_cog14_to_anvio.sh

bash 06_hmmscan_pfam.sh
bash 07_parse_pfam_domtblout_to_functions.sh
bash 08_import_pfam_to_anvio.sh
```

---

## 10) Verificación rápida en anvi’o

Para cada DB:

```bash
anvi-display-contigs-stats -c /media/mmorenos/disk3/010_fovi_metagenome/db_anvio8/cathy_db/chon.db
```

y revisa que existan annotation sources tipo `COG14_FUNCTION` y `PFAM32` (o según el nombre que importaste en `source`).
