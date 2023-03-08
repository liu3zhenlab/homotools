#!/usr/bin/perl -w
# ====================================================================================================
# File: msa2var.pl
# Author: Sanzhen Liu
# Date: 2/8/2023
# ====================================================================================================

use strict;
use warnings;
use Getopt::Long;

my $help;
my $pos_adj = 0;
my $base_tolerant = 0;
my $absence_char = "-";

sub prompt {
    print <<EOF;
    Usage: perl $0 <Input Fasta Files> [options]
    --msa <file>      : multipe sequence alignment output in a fasta format; required
    --ref <str>       : reference sequence name; first sequence by default
    --posadj <num>    : position adjust by adding this number ($pos_adj)
	--basefuzzy <num> : number of tolerant basepairs for INDEL mergeing ($base_tolerant)
    --abschar <char>  : "-" or "." for absence sequences ($absence_char)
	--help     : help information
EOF
exit;
}

# read parameters:
my ($msa, $ref);
&GetOptions("msa=s"       => \$msa,
            "ref=s"       => \$ref,
			"posadj=i"    => \$pos_adj,
			"basefuzzy=i" => \$base_tolerant,
			"abschar=s"   => \$absence_char,
			"help|h"      => \$help) || &prompt;

if ($help) { &prompt; }

if (!defined $msa) {
	print STDERR "--msa is required\n";
	exit;
}

if (($absence_char ne "-") and ($absence_char ne ".")) {
	print STDERR "--abschar must be - or .\n";
	exit;
}

### read alignment fasta data
my ($seq, $seqname, $len, %seqhash, %lenhash);
my (@seqnames, @seqlen);
open(IN, $msa) || die;
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		if (defined $seqname) {
			$seqhash{$seqname} = $seq;
			$lenhash{$seqname} = $len;
			push(@seqlen, $len);
		}
    	$seqname = $1;
		push(@seqnames, $seqname);
		$seq = "";
		$len = 0;
 	 } else {
		$seq .= $_;
		$len += length($_);
	}
}
# last element:
$seqhash{$seqname} = $seq;
$lenhash{$seqname} = $len;
push(@seqlen, $len);
close IN;

# sequence length
for (my $i=0; $i<$#seqlen; $i++) {
	if ($seqlen[$i+1] != $seqlen[$i]) {
		print STDERR "Not all sequences in $msa have the same length. Wrong FASTA formate. Exit!\n";
		exit;
	}
}

# check refname
if (!defined $ref) {
	$ref = $seqnames[0];
	print STDERR "By default, the 1st sequence was considered as the reference: $ref\n";
} elsif (!exists $seqhash{$ref}) {
	print STDERR "$ref is not in $msa. Exit!\n";
	exit;
}

### process aln fasta data
my (%ref_base, %ref_pos);
my %alt_base;
my %all_uniq_base;
foreach my $eachseq (@seqnames) {
	my $count = 0;
	if ($eachseq eq $ref) {
		my $ref_pos = 0;
		my @ref_base = split(//, $seqhash{$ref});
		foreach (@ref_base) {
			if ($_ ne $absence_char) {
				$ref_pos++;
			}
			$count++;
			$ref_base{$count} = $_;
			$ref_pos{$count} = $ref_pos;
			$all_uniq_base{$count}{$_}++;
		}
	} else {
		my @alt_base = split(//, $seqhash{$eachseq});
		foreach (@alt_base) {
			$count++;
			$alt_base{$count}{$eachseq} = $_;
			$all_uniq_base{$count}{$_}++;
		}
	}
}

### output genotyping data on polymorphic sites
# header
print "refname\trefpos\treftype\tref";
foreach my $eseqname (@seqnames) {
	if ($eseqname ne $ref) {
		print "\t$eseqname";
	}
}
print "\n";

# genotyping data
my $seqlen = $seqlen[0];
my %polymorphic_site_ref;
my %polymorphic_site_alt;
for (my $pos=1; $pos<=$seqlen; $pos++) {
	my %cur_bases = %{$all_uniq_base{$pos}};
	my @cur_bases = keys %cur_bases;
	if ($#cur_bases > 0) { # polymorphic sites
		my $ref_base_feature = ".";
		if ($ref_base{$pos} eq $absence_char) {
			$ref_base_feature = "ins";
		}
		#print "$ref\t$ref_pos{$pos}\t$ref_base_feature\t$ref_base{$pos}";
		foreach my $eseqname (@seqnames) {
			if ($eseqname ne $ref) {
				my $cur_alt_base = $alt_base{$pos}{$eseqname};
				#print "\t$cur_alt_base";
				if ($cur_alt_base ne $ref_base{$pos}) {
					$polymorphic_site_ref{$eseqname}{$ref_pos{$pos}} .= $ref_base{$pos};
					$polymorphic_site_alt{$eseqname}{$ref_pos{$pos}} .= $cur_alt_base;
				}
			}
		}
	}
}

# genotyping by sample
foreach my $eseqname (@seqnames) {
	if ($eseqname ne $ref) {
		if (exists $polymorphic_site_ref{$eseqname}) {
			my %epolymorph_ref = %{$polymorphic_site_ref{$eseqname}};
			my %epolymorph_alt = %{$polymorphic_site_alt{$eseqname}};
			my @poly_refpos = sort {$a <=> $b} keys %epolymorph_ref;
			my $previous_pos = $poly_refpos[0];
			my $previous_ref = "";
			my $previous_alt = "";
			my $cur_ref = "";
			my $cur_alt = "";
			for (my $i=0; $i<=$#poly_refpos; $i++) {
				my $cur_ref_base = $epolymorph_ref{$poly_refpos[$i]};
				my $cur_alt_base = $epolymorph_alt{$poly_refpos[$i]};
				$cur_ref .= $cur_ref_base;
				$cur_alt .= $cur_alt_base;
				
				if (exists $epolymorph_ref{$poly_refpos[$i] + 1}) { # neighboring polymorphism exists
					$previous_ref = $cur_ref;
					$previous_alt = $cur_alt;
				} else {
					$cur_ref =~ s/^\-+$/\./g;
					$cur_ref =~ s/\-+//g;
					$cur_alt =~ s/^\-+$/\./g;
					$cur_alt =~ s/\-+//g;
					print "$eseqname\t$ref\t$previous_pos\t$cur_ref\t$cur_alt\n";
					$previous_ref = "";
					$previous_alt = "";
					$cur_ref = "";
					$cur_alt = "";
					if ($i < $#poly_refpos) {
						$previous_pos = $poly_refpos[$i + 1];
					}
				}
			}
		}
	}
}


#########################################################

