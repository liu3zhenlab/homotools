perl ~/software/simuG/simuG.pl -r Zm00001eb001720.fasta -snp_count 20  -indel_count 5 -prefix allele1
perl ~/software/simuG/simuG.pl -r Zm00001eb001720.fasta -snp_count 200  -indel_count 5 -prefix allele2
perl ~/software/simuG/simuG.pl -r allele2.simseq.genome.fa -snp_count 5 -prefix allele3
