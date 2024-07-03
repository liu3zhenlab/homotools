#!/bin/bash
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00

qrygene=Zm00001eb214410

#. "/homes/liu3zhen/anaconda3/etc/profile.d/conda.sh"
#export PATH="/homes/liu3zhen/anaconda3/bin:$PATH"
conda activate homotools

outdir=hmout

if [ ! -d $outdir ]; then
	mkdir $outdir
fi

pushd $outdir
perl ~/scripts2/homotools/homomine \
	--qrygene $qrygene \
	--qrydir ~/references/maizeCurGenomes/DBs/B73-5.0 \
	--qrybase B73-5.0 \
	--tgtdir ~/references/maizeCurGenomes/DBs/A188v1 \
	--tgtbase A188v1 \
	--cleanup
popd

