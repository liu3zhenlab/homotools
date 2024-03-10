#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

##############################
# vcfeff.ann.parser.pl
# Sanzhen Liu
# 2/14/2023
##############################

my $anncol = 8;
my $format = "long";
sub prompt {                                                                      
    print <<EOF;
    Usage: perl $0 --vcfeff <vcf_snpeff> [options]
    [Options]
    --vcfeff <file>  : input file in a VCF format with snpEff annotation; required
    --anncol <num>   : column number for snpEff annotation ($anncol)
    --format <str>   : output format (long or short)
    --comment        : keep VCF comments if specified
    --singleletter   : convert 3 letter amino acid notation to single letter notation; off by default
    --help           : help information
EOF
exit;
}
my %opts = ();
&GetOptions(\%opts,
			"vcfeff=s",
			"anncol=i",
			"format=s",
			"comment",
			"singleletter",
			"help");


### parameters
&prompt if exists $opts{help} or !%opts;

my $vcfeff;
if (!exists $opts{vcfeff}) { 
	print STDERR "--vcfeff is required\n";
} else {
	$vcfeff = $opts{vcfeff};
}

if (exists $opts{format}) {
	if (($opts{format} ne "long") and ($opts{format} ne "short")) {
		print STDERR "--format must be either long or short\n";
		exit;
	}
	$format = $opts{format};
}

### read data
# ANN=C|missense_variant&splice_region_variant|MODERATE|EXON_Zm00001eb320870_1001_1200|Zm00001eb320870|transcrip
# t|Zm00001eb320870_T002|protein_coding|4/8|c.172G>C|p.Val58Leu|380/942|172/360|58/119||

# header
if ($format eq "short") {
	print "seq\tpos\tREF\tALT\tgene\ttranscript\teffect\timpact\tprotVar\tprotPos\n";
} else {
	#print "\t$gene\t$transcript\t$effect\t$impact\t$dna_notation\t$prot_notation\t$cdna_pos\t$cds_pos\t$prot_pos\t$distance\n";
	print "seq\tpos\tREF\tALT\tgene\ttranscript\teffect\timpact\tdnaVar\tprotVar\tcdnaPos\tcdsPos\tproPos\tdistance2gene\n";
}

open(IN, $vcfeff) || die;
while (<IN>) {
	chomp;
	if (/^#/) {
		if (exists $opts{comment}) {
			print "$_\n";
		}
	} else {
		my @line = split(/\t/, $_);
		my $annot = $line[$anncol - 1]; # annotation
		my @annot = split(/\|/, $annot);
		for (my $i=0; $i<=$#annot; $i++) {
			if ($annot[$i] eq "") {
				$annot[$i] = ".";
			}
		}
		my $effect = $annot[1];
		my $impact = $annot[2];
		my $gene = $annot[4];
		$gene = defined $gene ? $gene : '.';
		my $geneid = $annot[5];
		$geneid = defined $geneid ? $geneid : '.';
		my $transcript = $annot[6];
		$transcript = defined $transcript ? $transcript : '.';
		my $rank = $annot[8];
		$rank = defined $rank ? $rank : '.';
		my $dna_notation = $annot[9];
		$dna_notation = defined $dna_notation ? $dna_notation : '.';
		$dna_notation =~ s/^c\.//;
		my $prot_notation = ".";
		my $cdna_pos = ".";
		my $cds_pos = ".";
		my $prot_pos = ".";
		my $distance = ".";
		if ($#annot > 9) {
			$prot_notation = $annot[10];
			$prot_notation =~ s/^p\.//;
			if (exists $opts{singleletter}) {
				if ($prot_notation =~ /^([a-zA-Z]{3})([0-9]+)([a-zA-Z]{3})/) {
					my $ori_aa = $1;
					my $aa_pos = $2;
					my $new_aa = $3;
					my $ori_aa_single = convert_to_single_letter($ori_aa);
					my $new_aa_single = convert_to_single_letter($new_aa);
					$prot_notation = $ori_aa_single.$aa_pos.$new_aa_single;
				}
			}
		}
		if ($#annot > 10) {
			$cdna_pos = $annot[11];
		}
		if ($#annot > 11) {
			$cds_pos = $annot[12];
		}
		if ($#annot > 12) {
			$prot_pos = $annot[13];
		}
		if ($#annot > 13) {
			$distance = $annot[14];
		}
		print join("\t", @line[0..3]);
		if ($format eq "short") {
			print "\t$gene\t$transcript\t$effect\t$impact\t$prot_notation\t$prot_pos\n";
		} else {
			print "\t$gene\t$transcript\t$effect\t$impact\t$dna_notation\t$prot_notation\t$cdna_pos\t$cds_pos\t$prot_pos\t$distance\n";
		}
	}
}

close IN;

sub convert_to_single_letter {
# to convert three-letter aa notation to single-letter notation
# developed by ChatGPT and modified by Sanzhen Liu on 12/11/2023
	my $three_letter_code = shift;
	# Define the conversion table
	my %conversion_table = (
	'ALA' => 'A', 'ARG' => 'R', 'ASN' => 'N', 'ASP' => 'D',
	'CYS' => 'C', 'GLN' => 'Q', 'GLU' => 'E', 'GLY' => 'G',
	'HIS' => 'H', 'ILE' => 'I', 'LEU' => 'L', 'LYS' => 'K',
	'MET' => 'M', 'PHE' => 'F', 'PRO' => 'P', 'SER' => 'S',
	'THR' => 'T', 'TRP' => 'W', 'TYR' => 'Y', 'VAL' => 'V'
	);
	# Convert three letters to uppercase
	$three_letter_code = uc($three_letter_code);
	# Perform the conversion
	my $single_letter_code = $conversion_table{$three_letter_code};
	return defined $single_letter_code ? $single_letter_code : 'X'; # 'X' for unknown
}

#"ANN[*].CDNA_POS" (alias POS_CDNA)
# 87 #"ANN[*].CDNA_LEN" (alias LEN_CDNA)
# 88 #"ANN[*].CDS_POS" (alias POS_CDS)
# 89 #"ANN[*].CDS_LEN" (alias LEN_CDS)
# 90 #"ANN[*].AA_POS" (alias POS_AA)
# 91 #"ANN[*].AA_LEN" (alias LEN_AA)
# 92 #"ANN[*].DISTANCE"
#
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


