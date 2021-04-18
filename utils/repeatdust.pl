#!/usr/bin/perl -w

# ===============================================================
# repeatdust.pl
# Sanzhen Liu
# 4/17/2021
# 
# usage: perl repeatdust.pl <blast_output> [options]
#
# ===============================================================

use strict;
use warnings;
use Getopt::Long;

my $match = 100;
my $identity = 80;
my $coverage = 80;
my $nrepeat = 5;
my $identity_col = 6;
my $start_col = 8;
my $end_col = 9;
my $help;

sub errINF {
  print <<EOF

  Usage: perl repeatdust.pl <blast_output> [options]
  - based on BLASTN alignments, compare query segments each other to identify repeats
  - identity and coverage are used to define overlapped query segments in different alignments
  Options:
  --match|m <num>    : minimum bp match ($match)
  --identity|i <num> : percentage of identity (0-100) ($identity)
  --coverage|c <num> : coverage (0-100) ($coverage)
  --nrepeat|n <num>  : number of repeated alignments ($nrepeat)
  --identcol|d <num> : column number of identity ($identity_col)
  --startcol|s <num> : column number of query start coordinates ($start_col)
  --endcol|e <num>   : column number of query end coordinates ($end_col)
  --help             : Reminding informatioin
	
EOF
}

GetOptions("match|m=i"    => \$match,
           "identity|i=i" => \$identity,
           "coverage|c=i" => \$coverage,
           "nrepeat|n=i"  => \$nrepeat,
           "identcol|d=i" => \$identity_col,
		   "startcol|s=i" => \$start_col,
           "endcol|e=i"   => \$end_col,
           "help"         => \$help);

# judge parameters input
if ($help or @ARGV<1) {
	&errINF;
	exit;
}

# blastn output:
my (@starts, @ends);
open (IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	push(@starts, $line[$start_col - 1]);
	push(@ends, $line[$end_col - 1]);
}
close IN;

# open Query file:
open (IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	my $start = $line[$start_col - 1];
	my $end = $line[$end_col - 1];
	my $repeat_count = -1;
	for (my $i=0; $i<=$#starts; $i++) {
		my $is_overlap_alignment = &interval_overlap($start, $end, $starts[$i], $ends[$i], $match, $coverage);
		$repeat_count += $is_overlap_alignment;
		if ($repeat_count > $nrepeat) {
			last;
		}
	}
	
	if ($line[$identity_col - 1] >= $identity and $repeat_count <= $nrepeat) {
		print "$_\n";
	}
}
close IN;

sub interval_overlap {
# input includes 4 numbers
# first two are increasingly ordered numbers (first interval)
# second two are increasingly ordered numbers (2nd interval)
# identity
# coverage
	my $is_repeat = 0;
	my @input = @_; # size numbers
	my $match_cutoff = $input[4];
	my $coverage_cutoff = $input[5];

	my $overlap = 0;
	if (($input[0]>=$input[2]) & ($input[0]<=$input[3])) {
		$overlap = 1;
	} elsif (($input[1]>=$input[2]) & ($input[1]<=$input[3])) {
		$overlap = 1;
	} elsif (($input[0]<=$input[2]) & ($input[1]>=$input[3])) {
		$overlap = 1;
	}   	

# overlapping length and percentage
	my $overlap_length = 0;
	my $overlap_perc = 0;
	if ($overlap) {
		my @sort_input = sort {$a <=> $b} @input[0..3];
		$overlap_length = ($sort_input[2] - $sort_input[1] + 1);
		$overlap_perc = 100 * $overlap_length / ($input[1] - $input[0] + 1);
	}

# output
	if ($overlap_length >= $match_cutoff and $overlap_perc >= $coverage_cutoff) {
		$is_repeat = 1;
	}
	return($is_repeat);
}

