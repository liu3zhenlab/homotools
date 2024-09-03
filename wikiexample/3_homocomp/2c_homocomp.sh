#!/bin/bash

query=../2_geneseq/1o_geneseq/1o_geneseq.1.Zm00001eb001720.fasta
query_bed=../2_geneseq/1o_geneseq/1o_geneseq.4.pos.adjusted.gtf.bed/Zm00001eb001720_T001.adjusted.bed
ref_db=../1_data/blastDB/A188.fasta
ref_name=A188
target_gtf=../1_data/gtf/A188.gtf
prefix=2o_homocomp

perl ../../homocomp \
	--prefix $prefix \
	--query $query \
	--qryadd $query_bed \
	--db $ref_db \
	--dbacc $ref_name \
	--tgtf $target_gtf

