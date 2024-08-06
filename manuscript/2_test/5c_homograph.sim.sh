#!/bin/bash
perl /homes/liu3zhen/scripts2/homotools/homograph \
	--ref ../1_data/gene/Zm00001eb001720.fasta \
	--fasdir ../1_data/gene/1_simulation \
	--cdhitpara "-g 1 -s 0.98 -c 0.98 -r 0" \
	--prefix 5o_sim

