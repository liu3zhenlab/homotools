#!/usr/bin/perl -w
#
# ====================================================================================================
# File: fastaSize.pl
# Author: Sanzhen Liu
#
# This file count the size of each sequence in the input file if not specify the option --all. 
# If the option --all is added, count the total size of all sequences.
# If the option --average is specified, the average size will be informed.
# 
# Usage: perl fastaSize.pl <Input Fasta File> [--all] [--average]
# --all is an option, if it isn't specified, size of all sequences will be counted seperately.
# --average is another option, if it is specified, it count the average size for all sequences.
# --help
# ====================================================================================================

use strict;
use warnings;
use Getopt::Long;

my %seq_size;
my $seq_name;
my $size;
my $all; # default value of $all;
my $average; 
my $total;
my $count=0;
my $help;

sub prompt {
	print <<EOF;
	Usage: perl fastaSize.pl <Input Fasta File> [--all] [--average] [--help]
	Count the size for each fasta sequence or count the total size for all sequences in the file.
	--all is an option, if it isn't specified, size of each sequence will be counted seperately.
	--average is another option, if it is specified, it count the average size for all sequences.
	--help 
EOF
exit;
}
# read the parameters:
&GetOptions("all" => \$all, "average" => \$average, "help" => \$help) || &prompt;

if ($help) {
	&prompt;
}

open(IN, $ARGV[0]) || die "The input fasta file cannot be opened.";

# Read all sequence (name and size) into hash;
while (<IN>) {
   $_ =~ s/\R//g;
   chomp;
   if (/^>(\S+)/) {
      if (defined $seq_name) {
         $seq_size{$seq_name} = $size;
      }
      $seq_name = $1;
      $size = 0;
	  $count++;
   }
   else {
      $_ =~ s/\s//g;
	  $size += length($_);
   }
}
# last element:
$seq_size{$seq_name} = $size;
close IN;

# output according to the option input:
if ($all) {
	foreach (keys %seq_size) {
		$total += $seq_size{$_};
	}
	print "The TOTAL size of all sequences in \"$ARGV[0]\" is: $total\n";
} 

if ($average) {
	$total = 0;
	foreach (keys %seq_size) {
		$total += $seq_size{$_};
	}
	my $mean = $total/$count;
	print "The AVERAGE size of all sequences in \"$ARGV[0]\" is: $mean\n";
}

# if both --all and --average are not specified, print size of each sequence:
unless ($all || $average) {
	foreach (keys %seq_size) {
   		print "$_\t$seq_size{$_}\n";
	}
}
