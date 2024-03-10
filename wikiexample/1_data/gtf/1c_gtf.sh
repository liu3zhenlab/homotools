awk '$1=="chr1" && $4<10000000' ~/references/maizeCurGenomes/gtf/B73-5.0.gtf > B73.gtf
awk '$1=="chr1" && $4<10000000' ~/references/maizeCurGenomes/gtf/CML333-1.0.gtf > CML333.gtf
awk '$1=="chr1" && $4<10000000' ~/references/maizeCurGenomes/gtf/Ki3-1.0.gtf > Ki3.gtf
awk '$1=="chr1" && $4<10000000' ~/references/maizeCurGenomes/gtf/P39-1.0.gtf > P39.gtf
awk '$1==1 && $4<10000000' ~/references/maizeCurGenomes/gtf/A188Ref1.gtf > A188.gtf
sed -i 's/^1/chr1/g' A188.gtf
