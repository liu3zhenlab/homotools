#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --time=1-00:00:00

module load R

# input data
gene=Zm00035ab147480
ext5=2000
ext3=2000

source_ref=/homes/liu3zhen/references/NAM1.0/Ms71-1.0/Ms71-1.0.fasta
source_gtf=/homes/liu3zhen/references/maizeCurGenomes/gtf/Ms71-1.0.gtf
source_seq_dir=MAD69
te_gff=/homes/liu3zhen/references/NAM1.0/TEs/Zm-Ms71-REFERENCE-NAM-1.0.TE.gff3

# fas extraction
/homes/liu3zhen/scripts2/homotools/geneseq --fas $source_ref \
	--gtf $source_gtf \
	--othergff $te_gff \
	--gene $gene \
	--prefix ${source_seq_dir} \
	--ext5 $ext5 \
	--ext3 $ext3

