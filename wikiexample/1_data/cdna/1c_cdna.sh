#!/bin/bash

sed 's/.*transcript_id \"//g' ../gtf/A188.gtf | sed 's/\".*//g' | sort | uniq > A188.transcripts
seqtk subseq ~/references/A188Ref1/confident/A188Ref1a1.confident.transcripts.fasta A188.transcripts > A188.cdna.fasta
rm A188.transcripts

sed 's/.*transcript_id \"//g' ../gtf/B73.gtf | sed 's/\".*//g' | sort | uniq > B73.transcripts
seqtk subseq ~/references/maizeCurGenomes/cdna_b/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cdna.fa B73.transcripts > B73.cdna.fasta
rm B73.transcripts

sed 's/.*transcript_id \"//g' ../gtf/P39.gtf | sed 's/\".*//g' | sort | uniq > P39.transcripts
seqtk subseq ~/references/maizeCurGenomes/cdna_b/Zm-P39-REFERENCE-NAM-1.0_Zm00040ab.1.cdna.fa P39.transcripts > P39.cdna.fasta
rm P39.transcripts

sed 's/.*transcript_id \"//g' ../gtf/CML333.gtf | sed 's/\".*//g' | sort | uniq > CML333.transcripts
seqtk subseq ~/references/maizeCurGenomes/cdna_b/Zm-CML333-REFERENCE-NAM-1.0_Zm00026ab.1.cdna.fa CML333.transcripts > CML333.cdna.fasta
rm CML333.transcripts

sed 's/.*transcript_id \"//g' ../gtf/Ki3.gtf | sed 's/\".*//g' | sort | uniq > Ki3.transcripts
seqtk subseq ~/references/maizeCurGenomes/cdna_b/Zm-Ki3-REFERENCE-NAM-1.0_Zm00029ab.1.cdna.fa Ki3.transcripts > Ki3.cdna.fasta
rm Ki3.transcripts
