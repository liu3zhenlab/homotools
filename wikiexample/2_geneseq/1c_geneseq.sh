# input data
gene=Zm00001eb001720
source_ref=../1_data/genome/B73.fasta
source_gtf=../1_data/gtf/B73.gtf
source_seq_dir=1o_geneseq

# fas extraction
../../geneseq \
	--fas $source_ref \
	--gtf $source_gtf \
	--gene $gene \
	--prefix $source_seq_dir

