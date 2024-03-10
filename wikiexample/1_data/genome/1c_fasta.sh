#!/bin/bash
# extract sequence
bedtools getfasta -nameOnly -fi ~/references/maizeCurGenomes/genomes/B73-5.0.fasta -bed 1i_region.bed -fo B73.fasta
bedtools getfasta -nameOnly -fi ~/references/maizeCurGenomes/genomes/P39-1.0.fasta -bed 1i_region.bed -fo P39.fasta
bedtools getfasta -nameOnly -fi ~/references/maizeCurGenomes/genomes/Ki3-1.0.fasta -bed 1i_region.bed -fo Ki3.fasta
bedtools getfasta -nameOnly -fi ~/references/maizeCurGenomes/genomes/CML333-1.0.fasta -bed 1i_region.bed -fo CML333.fasta

sed 's/^chr//g' 1i_region.bed > 1i_A188.region.bed
bedtools getfasta -nameOnly -fi ~/references/maizeCurGenomes/genomes/A188Ref1.fasta -bed 1i_A188.region.bed -fo A188.fasta


