#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16

. "/homes/liu3zhen/anaconda3/etc/profile.d/conda.sh"
conda activate homotools
datadir=/homes/liu3zhen/scripts2/homotools/wikiexample/1_data/homograph

perl /homes/liu3zhen/scripts2/homotools/homograph \
	--genebed $datadir/0_ref/B73.bed \
	--msatool clustalo \
	--threads 16 \
	--ref $datadir/0_ref/B73.fasta \
	--fasdir $datadir/1_fas

