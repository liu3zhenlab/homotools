#!/usr/bin/perl -w

# ===============================================================
# redundant.aln.rm.pl
# Sanzhen Liu
# 2/19/2023
# 
# usage: perl redundant.aln.rm.pl <blast_output> [options]
#
# ===============================================================

use strict;
use warnings;
use Getopt::Long;

# qseqid qlen sseqid slen length pident qcovs qstart qend sstart send evalue bitscore

my $min_match = 50;
my $min_coverage = 80;
my $start_col = 8;
my $end_col = 9;
my $bitscore_col = 13;
my $help;

sub errINF {
  print <<EOF

  Usage: perl repeatdust.pl <blast_output> [options]
  - based on BLASTN alignments, compare query segments each other to identify repeats
  - identity and coverage are used to define overlapped query segments in different alignments
  Options:
  --minmatch|m <num> : minimum bp match ($min_match)
  --mincov|c <num>   : minimum coverage (0-100) of the overlapping region out of the smaller interval ($min_coverage)
  --startcol|s <num> : column number of query start coordinates ($start_col)
  --endcol|e <num>   : column number of query end coordinates ($end_col)
  --bscol|b <num>    : column number of bitscores ($bitscore_col)
  --help             : Reminding informatioin
EOF
}

GetOptions("minmatch|m=i" => \$min_match,
           "mincov|c=i"   => \$min_coverage,
		   "startcol|s=i" => \$start_col,
           "endcol|e=i"   => \$end_col,
		   "bscol|b=i"    => \$bitscore_col,
		   "help"         => \$help);

# judge parameters input
if ($help or @ARGV<1) {
	&errINF;
	exit;
}

if (-z $ARGV[0]) { # if an input is an empty file
	exit;
}

# blastn output:
my (@starts, @ends, @bitscore);
open (IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	push(@starts, $line[$start_col - 1]);
	push(@ends, $line[$end_col - 1]);
	push(@bitscore, $line[$bitscore_col -1]);
}
close IN;

# open Query file:
my $row = 0;
my %flag_row;
open (IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	my @line = split(/\t/, $_);
	my $start = $line[$start_col - 1];
	my $end = $line[$end_col - 1];
	my $repeat_count = -1;
	for (my $i=0; $i<=$#starts; $i++) {
		if (!exists $flag_row{$i}) {
			if ($row != $i) {
				my $is_overlap = &interval_overlap($start, $end, $starts[$i], $ends[$i], $min_match, $min_coverage);
				if ($is_overlap) {
					if ($bitscore[$row] >= $bitscore[$i]) {
						$flag_row{$i}++;
					} else {
						$flag_row{$row}++;
					}
				}
			}
		}
	}
	$row++;
}
close IN;

$row = 0;
open (IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	if (!exists $flag_row{$row}) {
		print "$_\n";
	}
	$row++;
}
close IN;

sub interval_overlap {
# input includes 4 numbers and 2 cutoffs
# first two are increasingly ordered numbers (first interval)
# second two are increasingly ordered numbers (2nd interval)
# two cutoffs are identity and coverage
#
# return overlapping codes:
# 0: no overlap
# 1: overlap
	my $overlap_code = 0; # no overlap
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

# overlapping length and percentage (coverage)
	my $overlap_length = 0;
	my $overlap_perc = 0;
	if ($overlap) {
		my @sort_input = sort {$a <=> $b} @input[0..3];
		$overlap_length = ($sort_input[2] - $sort_input[1] + 1);
		my $len_1 = $input[1] - $input[0] + 1;
		my $len_2 = $input[3] - $input[2] + 1;
		if ($len_1 >= $len_2) {
			$overlap_perc = 100 * $overlap_length / $len_2;
		} else {
			$overlap_perc = 100 * $overlap_length / $len_1;
		}
	}

# output
	if ($overlap_length >= $match_cutoff and $overlap_perc >= $coverage_cutoff) {
		$overlap_code = 1;
	}
	return($overlap_code);
}

