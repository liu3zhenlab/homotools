#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long
##############################
# vcfeff.ann.parser.pl
Sanzhen Liu
# 2/14/2023
##############################


#"ANN[*].ALLELE" (alias GENOTYPE)
#"ANN[*].EFFECT" (alias ANNOTATION): Effect in Sequence ontology terms (e.g. 'missense_variant', 'synonymous_variant', 'stop_gained', etc.)
#"ANN[*].IMPACT:" { HIGH, MODERATE, LOW, MODIFIER }
#"ANN[*].GENE:" Gene name (e.g. 'PSD3')
#"ANN[*].GENEID:" Gene ID
#"ANN[*].
#"ANN[*].FEATUREID" (alias TRID: Transcript ID)
#"ANN[*].BIOTYPE:" Biotype, as described by the annotations (e.g. 'protein_coding')
#"ANN[*].RANK:" Exon or Intron rank (i.e. exon number in a transcript)
#"ANN[*].HGVS_C" (alias HGVS_DNA, CODON): Variant in HGVS (DNA) notation
#"ANN[*].HGVS_P" (alias HGVS, HGVS_PROT, AA): Variant in HGVS (protein) notation
#"ANN[*].CDNA_POS" (alias POS_CDNA)
#"ANN[*].CDNA_LEN" (alias LEN_CDNA)
#"ANN[*].CDS_POS" (alias POS_CDS)
#"ANN[*].CDS_LEN" (alias LEN_CDS)
#"ANN[*].AA_POS" (alias POS_AA)
#"ANN[*].AA_LEN" (alias LEN_AA)
#"ANN[*].DISTANCE"
#"ANN[*].ERRORS" (alias WARNING, INFOS


