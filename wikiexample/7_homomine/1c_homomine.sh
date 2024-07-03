#!/bin/bash
#conda activate homotools
qrygene=Zm00001eb000510
outdir=hmout2
if [ ! -d $outdir ]; then
	mkdir $outdir
fi
pushd $outdir

datadir=../../1_data/homomine
perl ~/scripts2/homotools/homomine \
  --qrygene $qrygene \
  --qrydir $datadir/B73 \
  --qrybase B73 \
  --tgtdir $datadir/A188 \
  --tgtbase A188 \
  --cleanup
popd

