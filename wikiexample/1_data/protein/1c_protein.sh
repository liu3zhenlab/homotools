#!/bin/bash

sed 's/.*transcript_id \"//g' ../gtf/A188.gtf | sed 's/\".*//g' | sort | uniq > A188.proteins
seqtk subseq ~/references/A188Ref1/confident/A188Ref1a1.confident.proteins.fasta A188.proteins > A188.protein.fasta
rm A188.proteins

sed 's/.*transcript_id \"//g' ../gtf/B73.gtf | sed 's/\".*//g' | sort | uniq | sed 's/_T/_P/g' > B73.proteins
seqtk subseq ~/references/maizeCurGenomes/protein/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa B73.proteins > B73.protein.fasta
rm B73.proteins

sed 's/.*transcript_id \"//g' ../gtf/P39.gtf | sed 's/\".*//g' | sort | uniq | sed 's/_T/_P/g' > P39.proteins
seqtk subseq ~/references/maizeCurGenomes/protein/Zm-P39-REFERENCE-NAM-1.0_Zm00040ab.1.protein.fa P39.proteins > P39.protein.fasta
rm P39.proteins

sed 's/.*transcript_id \"//g' ../gtf/CML333.gtf | sed 's/\".*//g' | sort | uniq | sed 's/_T/_P/g' > CML333.proteins
seqtk subseq ~/references/maizeCurGenomes/protein/Zm-CML333-REFERENCE-NAM-1.0_Zm00026ab.1.protein.fa CML333.proteins > CML333.protein.fasta
rm CML333.proteins

sed 's/.*transcript_id \"//g' ../gtf/Ki3.gtf | sed 's/\".*//g' | sort | uniq | sed 's/_T/_P/g' > Ki3.proteins
seqtk subseq ~/references/maizeCurGenomes/protein/Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1.protein.fa Ki3.proteins > Ki3.protein.fasta
rm Ki3.proteins
