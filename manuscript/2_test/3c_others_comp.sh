#!/bin/bash
genolist="
CML333
Ki3
P39
"
query=B73/B73.1.Zm00001eb001720.fasta

for geno in $genolist; do
	ref_db=../1_data/${geno}/${geno}.fasta
	prefix=${geno}comp
	perl ../../homocomp \
		--prefix $prefix \
		--dbacc $geno \
		--query $query \
		--db $ref_db
done
