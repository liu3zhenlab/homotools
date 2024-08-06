# input data
source_ref=../1_data/B73/B73.fasta
source_gtf=../1_data/B73/B73.gtf
prefix=B73
gene=Zm00001eb001720

# fas extraction
perl ../../geneseq \
    --fas $source_ref \
    --gtf $source_gtf \
    --gene $gene \
    --prefix $prefix
