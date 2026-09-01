
nohup anvi-pan-genome -g fovi-GENOMES.db  --project-name "fovi_genome_pan" -o fovi_pan_blast-
```

PAN.db  --num-threads 55 --minbit 0.5 --mcl-inflation 2 --I-know-this-is-not-a-good-idea --force-overwrite & anvi-compute-functional-enrichment-in-pan -p ../pangenome_QC311/fovi_pan_blast-PAN.db -g ../pangenome_QC311/fovi-GENOMES.db --category-variable group  --annotation-source KOfam -o kofam-fovi.txt --force-overwrite --display-db-calls
 anvi-compute-functional-enrichment-in-pan -p ../pangenome_QC311/fovi_pan_blast-PAN.db -g ../pangenome_QC311/fovi-GENOMES.db --category-variable group  --annotation-source KOfam -o kofam-fovi.txt --force-overwrite --display-db-calls

```

| Functional annotation source | source                     | used                                                        |
| ---------------------------- | -------------------------- | ----------------------------------------------------------- |
| COG20_FUNCTION               | NCBI COGs                  | Función de proteínas                                        |
| COG20_CATEGORY               | NCBI COGs                  | Categorías funcionales COG                                  |
| COG20_PATHWAY                | NCBI COGs                  | Vías asociadas a COG                                        |
| KOfam                        | KEGG Ortholog / KOfam HMMs | KO, metabolismo                                             |
| KEGG_Module                  | KEGG                       | Módulos metabólicos                                         |
| KEGG_Class                   | KEGG                       | Clasificación funcional KEGG                                |
| PFAM                         | Pfam                       | Familias de proteínas                                       |
| CAZyme                       | dbCAN                      | Enzimas activas sobre carbohidratos                         |
| kegg-functions               | KEGG/KOfam                 | Anotaciones funcionales usadas por anvi-estimate-metabolism |
