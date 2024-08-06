#!/bin/bash
perl ../../homostack \
	--seq B73/B73.1.Zm00001eb001720.fasta \
	--annot B73/B73.4.pos.adjusted.gtf.bed/Zm00001eb001720_T001.adjusted.bed \
	--seq CML333comp/CML333comp.4.target.fas --annot none \
	--seq Ki3comp/Ki3comp.4.target.fas --annot none \
	--seq A188comp/A188comp.4.target.fas --annot none \
	--seq P39comp/P39comp.4.target.fas --annot none \
	--prefix hsOut


