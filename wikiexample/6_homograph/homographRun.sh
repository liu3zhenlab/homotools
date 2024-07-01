#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00

datadir=../1_data/homograph

perl /homes/liu3zhen/scripts2/homotools/homograph \
	--ref $datadir/0_ref/B73.fasta \
	--genebed $datadir/0_ref/B73.bed \
	--fasdir $datadir/1_fas \
	--threads 4 \
	--prefix hgrun

