#!/bin/bash
query=B73/B73.1.Zm00001eb001720.fasta
ref_db=../1_data/A188/A188.fasta
prefix=A188comp

perl ../../homocomp \
	--prefix $prefix \
	--dbacc A188 \
	--query $query \
	--db $ref_db
