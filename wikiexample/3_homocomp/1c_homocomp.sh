#!/bin/bash

query=../2_geneseq/1o_geneseq/1o_geneseq.1.Zm00001eb001720.fasta
ref_db=../1_data/blastDB/A188
prefix=1o_homocomp

perl ../../homocomp \
	--prefix $prefix \
	--query $query \
	--db $ref_db \

