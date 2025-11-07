## The Resistance Gene Identifier (RGI) and  Comprehensive Antibiotic Resistance Database (CARD) v 4.0.0 paper from 2023 db (https://academic.oup.com/nar/article/51/D1/D690/6764414)

## https://github.com/arpcard/rgi/blob/master/docs/rgi_main.rst
```

mamba activate rgi
rgi auto_load -i localDB/card.json
for sample in `cat list` ; do  rgi main --input_sequence $sample --output_file  $sample.rgi --local --clean -a DIAMOND -t protein -n 60 --low_quality --include_nudge ; done
/media/mmorenos/disk3/DSG123-selected/001_anvio_metag/Y_genomes
```
```

conda activate abricate_env
abricate --db card ilq_contigs.fasta.db.contigs.fasta >  ilq_card.tsv

#!/usr/bin/env bash
set -euo pipefail

# === Configura aquí las bases de datos que quieres usar ===
DBS=(argannot card ncbi ncbibetalactamase plasmidfinder resfinder)

# Sufijo exacto que se removerá para obtener el prefijo (p.ej., "ilq")
SUFFIX="_contigs.fasta.db.contigs.fasta"

# Si no se pasan argumentos, usar todos los FASTA que calcen el patrón
shopt -s nullglob
if [[ $# -gt 0 ]]; then
  INPUTS=("$@")
else
  INPUTS=(*"$SUFFIX")
fi

if [[ ${#INPUTS[@]} -eq 0 ]]; then
  echo "No encontré archivos con el patrón *$SUFFIX en el directorio actual." >&2
  exit 1
fi

# Obtener lista de DBs instaladas en abricate (primer campo de --list)
# No aborta si falta alguna: solo avisa y sigue con las demás.
mapfile -t ABR_DBS < <(abricate --list | awk 'NR>1 {print $1}')

for fasta in "${INPUTS[@]}"; do
  base=$(basename "$fasta")
  # Quita el sufijo para obtener el prefijo (ej.: ilq_contigs... -> ilq)
  if [[ "$base" != *"$SUFFIX" ]]; then
    echo "ADVERTENCIA: '$base' no termina en $SUFFIX; lo salto." >&2
    continue
  fi
  prefix=${base%"$SUFFIX"}

  for db in "${DBS[@]}"; do
    if printf '%s\n' "${ABR_DBS[@]}" | grep -qx "$db"; then
      echo "[$(date +%T)] Ejecutando: sample=${prefix} | db=${db}"
      abricate --db "$db" "$fasta" > "${prefix}_${db}.tsv"
    else
      echo "ADVERTENCIA: la base '$db' no aparece en 'abricate --list'. Me la salto." >&2
    fi
  done
done

echo "Listo: salidas *.tsv generadas por cada muestra x base."

```

```

https://github.com/xinehc/args_oap
nohup args_oap stage_one -i ./ -o output -f fa -t 60 &
*recordar que se analiza a partir de reads en formato fasta
```
