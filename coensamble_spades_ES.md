
# Co-ensamblaje metagenómico con SPAdes

Este documento describe los pasos realizados para el co-ensamblaje de datos metagenómicos utilizando **SPAdes** en modo metagenómico (*metaSPAdes*), partiendo de lecturas previamente filtradas contra una referencia (no-Macrohit).

## 1. Estructura de datos de entrada
Se trabajó con 12 muestras en total, organizadas en dos grupos:
- Grupo **H**: 6 pares de archivos R1/R2
- Grupo **Y**: 6 pares de archivos R1/R2

Los archivos tienen el siguiente formato de nombre:
```
<fecha>_<nombre>_non-macroHit_R1.fastq.gz
<fecha>_<nombre>_non-macroHit_R2.fastq.gz
```

## 2. Script de co-ensamblaje
El script detecta automáticamente los pares de lecturas para cada grupo y ejecuta **metaSPAdes** para generar dos co-ensambles (uno por grupo).

Parámetros clave utilizados:
- **metaSPAdes** activado con `--meta`
- K-mers: `21,33,55,77,99,127`
- Hilos: `30`
- Memoria: `200 GB`
- `--only-assembler` para omitir corrección de errores en esta etapa
- Directorio temporal en disco rápido con `--tmp-dir` (opcional)

## 3. Ejecución del script
Ejemplo de ejecución:
```bash
./coassemble_spades_v2.sh   -i /ruta/a/01_bbduk_ref   -o /ruta/a/spades_coassemblies   -t 30 -m 200   --kmers 21,33,55,77,99,127   --meta   --tmp-dir /ruta/a/tmp
```

## 4. Resultados generados
El script genera dos carpetas de salida:
```
spades_coassemblies/spades_H_coassembly/
spades_coassemblies/spades_Y_coassembly/
```
Cada carpeta contiene:
- `contigs.fasta` → Contigs ensamblados
- `scaffolds.fasta` → Scaffolds ensamblados
- `spades.log` → Registro completo de ejecución

## 5. Siguientes pasos recomendados
1. Mapear las lecturas originales contra los contigs para obtener coberturas.
2. Filtrar contigs por longitud (≥ 1–1.5 kb).
3. Realizar binning con **MetaBAT2**, **VAMB** u otra herramienta.
4. Anotar genes rRNA (opcional) usando **Barrnap** para verificar la presencia de genes ribosomales.
