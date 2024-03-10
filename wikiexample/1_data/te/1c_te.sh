awk '$1=="chr1" && $4<10000000' ~/references/NAM1.0/TEs/Zm-B73-REFERENCE-NAM-5.0.TE.gff3 > B73.te.gff3
awk '$1=="chr1" && $4<10000000' ~/references/NAM1.0/TEs/Zm-CML333-REFERENCE-NAM-1.0.TE.gff3 > CML333.te.gff3
awk '$1=="chr1" && $4<10000000' ~/references/NAM1.0/TEs/Zm-Ki3-REFERENCE-NAM-1.0.TE.gff3 > Ki3.te.gff3
awk '$1=="chr1" && $4<10000000' ~/references/NAM1.0/TEs/Zm-P39-REFERENCE-NAM-1.0.TE.gff3 > P39.te.gff3
awk '$1==1 && $4<10000000' /homes/liu3zhen/references/A188Ref1/repeats/A188Ref1.fasta.mod.EDTA.TEanno.gff > A188.te.gff3
sed -i 's/^1/chr1/g' A188.te.gff3
