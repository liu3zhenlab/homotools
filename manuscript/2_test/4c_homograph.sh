#!/bin/bash

datadir=4i_fas
if [ ! -d $datadir ]; then
	mkdir $datadir
fi

pushd $datadir
cp ../*comp/*4.target.fas .
# simplify sequence names:
for fas in *fas; do
	sed -i 's/_.*//g' $fas
done
popd

perl ../../homograph \
	--ref B73/B73.1.Zm00001eb001720.fasta \
	--genebed B73/B73.4.pos.adjusted.gtf.bed/Zm00001eb001720_T001.adjusted.bed \
	--fasdir $datadir \
	--prefix 4o_hg

