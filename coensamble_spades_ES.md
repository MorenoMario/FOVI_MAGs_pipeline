
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
/home/mmoreno/.local/share/mamba/envs/spades315/bin/spades.py     --only-assembler        -k      21,33,55,77,99,127      --threads       30      --memory
        200     -o      /datos1/mmorenos/fovi_data/01_bbduk_ref/spades_Y_coassembly     --pe1-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y1_S1_non-macroHit_R1.fastq.gz     --pe1-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y1_S1_non-macroHit_R2.fastq.gz     --pe2-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y2_S2_non-macroHit_R1.fastq.gz     --pe2-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y2_S2_non-macroHit_R2.fastq.gz     --pe3-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y3_S3_non-macroHit_R1.fastq.gz     --pe3-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y3_S3_non-macroHit_R2.fastq.gz     --pe4-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y4_S4_non-macroHit_R1.fastq.gz     --pe4-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y4_S4_non-macroHit_R2.fastq.gz     --pe5-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y5_S5_non-macroHit_R1.fastq.gz     --pe5-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y5_S5_non-macroHit_R2.fastq.gz     --pe6-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y6_S6_non-macroHit_R1.fastq.gz     --pe6-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_Y6_S6_non-macroHit_R2.fastq.gz 

/home/mmoreno/.local/share/mamba/envs/spades315/bin/spades.py     --only-assembler        -k      21,33,55,77,99,127      --threads       30      --memory
        200     -o      /datos1/mmorenos/fovi_data/01_bbduk_ref/spades_H_coassembly2    --pe1-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S7_non-macroHit_R1.fastq.gz     --pe1-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S7_non-macroHit_R2.fastq.gz     --pe2-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S8_non-macroHit_R1.fastq.gz     --pe2-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H1_S8_non-macroHit_R2.fastq.gz     --pe3-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S10_non-macroHit_R1.fastq.gz    --pe3-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S10_non-macroHit_R2.fastq.gz    --pe4-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S9_non-macroHit_R1.fastq.gz     --pe4-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H2_S9_non-macroHit_R2.fastq.gz     --pe5-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S11_non-macroHit_R1.fastq.gz    --pe5-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S11_non-macroHit_R2.fastq.gz    --pe6-1 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S12_non-macroHit_R1.fastq.gz    --pe6-2 /datos1/mmorenos/fovi_data/01_bbduk_ref/20250604_PamelaFernandez_DSG_H3_S12_non-macroHit_R2.fastq.gz

```

## 3. Resultados generados
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
