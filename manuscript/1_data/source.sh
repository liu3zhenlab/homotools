#!/bin/bash

genolist="
A188
B73
CML333
Ki3
P39
"

for geno in $genolist; do
	if [ ! -d $geno ]; then
		mkdir $geno
	fi
	cp ../../wikiexample/1_data/genome/${geno}.fasta $geno/
	cp ../../wikiexample/1_data/gtf/${geno}.gtf $geno/
	cp ../../wikiexample/1_data/protein/${geno}.protein.fasta $geno/
	cp ../../wikiexample/1_data/cds/${geno}.cds.fasta $geno/
	cp ../../wikiexample/1_data/cdna/${geno}.cdna.fasta $geno/
done

# make blast+ database
for geno in $genolist; do
	pushd $geno
	makeblastdb -in ${geno}.fasta -dbtype nucl
	popd
done

