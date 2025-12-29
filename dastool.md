```
nohup DAS_Tool --write_bins --write_unbinned --search_engine blastp -i concoct.new_Y256.contigs2bin.tsv,maxbin_Y256_contigs2bin.tsv,metabat_Y256_contigs2bin.tsv -l concoct,maxbin,metabat -c ../../../02_spades_ass/spades_Y256_coassembly/Y_S256_contigs.fasta -o Y256_binning_mash_5  -t 10 --score_threshold 0.5 &> dastool.5.out.txt &

nohup DAS_Tool --write_bins --write_unbinned --search_engine blastp -i concoct.new_Y256.contigs2bin.tsv,maxbin_Y256_contigs2bin.tsv,metabat_Y256_contigs2bin.tsv -l concoct,maxbin,metabat -c ../../../02_spades_ass/spades_Y256_coassembly/Y_S256_contigs.fasta -o Y256_binning_mash_4  -t 10 --score_threshold 0.4 &> dastool.4.out.txt &

nohup DAS_Tool --write_bins --write_unbinned --search_engine blastp -i concoct.new_Y256.contigs2bin.tsv,maxbin_Y256_contigs2bin.tsv,metabat_Y256_contigs2bin.tsv -l concoct,maxbin,metabat -c ../../../02_spades_ass/spades_Y256_coassembly/Y_S256_contigs.fasta -o Y256_binning_mash_3  -t 10 --score_threshold 0.3 &> dastool.3.out.txt &

```
