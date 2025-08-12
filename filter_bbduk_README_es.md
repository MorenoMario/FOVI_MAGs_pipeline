# Filtrado de lecturas contra referencia de *Macrocystis pyrifera* (BBDuk)

Este script filtra lecturas **paired-end** contra el genoma de referencia del alga *Macrocystis pyrifera* usando `bbduk.sh`.  
Permite dos modos de trabajo:

1. **Guardar ambos conjuntos**: lecturas que matchean (hits) y las que no (no-hits).
2. **Guardar solo no-hits**: descartar lecturas que matchean a la referencia.

---

## Requisitos

- Linux/macOS con `bash`
- [BBTools / BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) en el `PATH`
- `gzip`

Verifica:
```bash
bbduk.sh --version
```

---

## Ejemplos de uso

### 1️⃣ Guardar ambos conjuntos (hits y no-hits)
```bash
./filter_host_bbduk.sh   -i ./00_skewer/   -o ./01_bbduk_ref/   -r ./GCA_031763025.fna   -t 32   -m 150g
```

Esto generará:
- `<sample>_macroHit_R1.fastq.gz` / `_R2.fastq.gz` (hits)
- `<sample>_non-macroHit_R1.fastq.gz` / `_R2.fastq.gz` (no-hits)
- `logs/<sample>.bbduk.stats.txt`
- `manifest_bbduk.tsv`

### 2️⃣ Guardar solo no-hits
Modifica el script para establecer:
```
outm=/dev/null
outm2=/dev/null
```
o crea una copia del script con esa configuración.

Ejemplo:
```bash
./filter_host_bbduk.sh   -i ./00_skewer/   -o ./01_bbduk_ref/   -r ./GCA_031763025.fna   -t 32   -m 150g   --force
```

### 3️⃣ Ejecutar en segundo plano con `nohup`
```bash
nohup ./filter_host_bbduk.sh   -i ./00_skewer/   -o ./01_bbduk_ref/   -r ./GCA_031763025.fna   -t 60   -m 150g   --force   > filter_bbduk.out 2>&1 &
tail -f filter_bbduk.out
```

---

## Entradas en carpetas distintas

Si tus `*-trimmed-pair1.fastq.gz` y `*-trimmed-pair2.fastq.gz` están dentro de subcarpetas (por ejemplo):

```
/media/mmorenos/disk3/DSG123-selected/rawReads/00_skewer/20250604_PamelaFernandez_DSG_H1_S7/
```

Debes permitir que `find` busque de forma recursiva. En el script, ajusta la búsqueda de R1:

```bash
mapfile -t R1_FILES < <(find "$INPUT_DIR" -type f \
  \( -name "*-trimmed-pair1.fastq.gz" -o -name "*-trimmed-pair1.fastq" \) | sort)
```

Así procesará archivos dentro de cualquier subcarpeta de `INPUT_DIR`.

