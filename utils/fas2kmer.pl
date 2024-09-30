#!/usr/bin/perl -w
# ==============================================================
# File: fas2kmer.pl
# Author: Sanzhen Liu
# Creation: 7/18/2021
# Modification: 9/25/2024
# ==============================================================

use strict;
use warnings;
use Getopt::Long;

my ($seq, $seq_name);
my ($sort, $help);

# default
my $ksize = 13;
my $prefix = "kout";

sub prompt {
    print <<EOF;
Usage: perl $0 [options] <fasta>
  --ksize|k <num> : length of mer ($ksize)
  --min|m <num>   : minimum count of a k-mer ($min_count)
  --prefix|p <str>: prefix of outputs ($prefix)
  --help|h        : help information
EOF
exit;
}
# read the parameters:
&GetOptions("ksize|k=i"   => \$ksize,
			"min|m=i"     => \$min_count,
			"prefix|p=s"  => \$prefix,
			"help|h"      => \$help) || &prompt;

if ($help) { &prompt; }
if (@ARGV<1) { &prompt; }

my $allout = $prefix.".".$ksize.".mer.all";
if ($mode ne "seq") {
	open(ALL, ">", $allout);
	print ALL "Kmer\t$allhead\n";
	close ALL;
}

my %all_khash;

### main
foreach my $input (@ARGV) {
	open(IN, "<", $input) || die;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			if (defined $seq_name) {
				&seq2k($seq_name, $seq, $ksize, $sepout);
			}
  		  	$seq_name = $1;
			$seq = "";
	 	 } else {
			$seq .= $_;
			$seq = uc($seq);
		}
	}
	# last element:
	&seq2k($seq_name, $seq, $ksize, $sepout);
	close IN;
	undef $seq_name;
	undef $seq;
#close SEP;
}
### output all set
if ($mode ne "seq") {
	open(ALL, ">>", $allout);
	if ($sort) {
		foreach (sort {$a cmp $b} keys %all_khash) {
			if ($all_khash{$_} >= $min_count) {
				print ALL "$_\t$all_khash{$_}\n";
			}
		}
	} else {
		foreach (keys %all_khash) {
			if ($all_khash{$_} >= $min_count) {
				print ALL "$_\t$all_khash{$_}\n";
			}
		}
	}
	close ALL;
}

############################################################
### modules
############################################################
sub seq2k {
# input: seqname, sequence and k-mer size
# output: add kmers as the key and counts as the hash value
#         to the global hash variable khash
	my ($inseqname, $inseq, $inksize, $insepout) = @_;
	my %seqkhash = ();
	while (my $kmer_1 = substr($inseq, 0, 1, "")) {
		my $kmer_rest = substr($inseq, 0, $inksize - 1);
		my $kmer = $kmer_1.$kmer_rest;
		if (($kmer =~ /^[ATGCatgc]+$/) and (length($kmer) == $inksize)) {
		# no N and length = $inksize
			my $kmer_rc = &revcom($kmer);
			my @kmer_sort = sort {$a cmp $b} ($kmer, $kmer_rc);
			my $out_kmer = $kmer_sort[0];
			$seqkhash{$out_kmer}++;
			$all_khash{$out_kmer}++;
		}
	}
	# output
	open(SEPOUT, ">>", $insepout) || die;
	if (defined $sort) {
		foreach (sort {$a cmp $b} keys %seqkhash) {
			print SEPOUT "$inseqname\t$_\t$seqkhash{$_}\n";
		}
	} else {
		foreach (keys %seqkhash) {
			print SEPOUT "$inseqname\t$_\t$seqkhash{$_}\n";
		}
	}
	close SEPOUT;
}

sub revcom {
### reverse and complement an input sequence
	my $inseq = shift @_; 
	my $revcom = reverse($inseq);
	$revcom =~ tr/AGCTagct/TCGAtcga/;
	return $revcom;
}

