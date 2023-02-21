#!/usr/bin/perl -w
#======================================================================
# cdhit2pggb.pl 
#
# Author: Sanzhen Liu <liu3zhen@ksu.edu>
# 11/27/2021
#
# to reformat fasta data based on cd-hit clustering result
#======================================================================

use strict;
use warnings;
use Getopt::Long;

my $gene = "seq"; 
my $prefix = "postcdhit";

sub prompt {
    print <<EOF;
    Usage: perl $0 --fasta <fasta> --clust <cd-hit cluster> [options]
    [Options]
    --fasta <file>  fasta file containing all sequence in cd-hit clusters; required
    --clust <file>  cluster output from cd-hit-est; required
    --prefix <str>  output prefix ($prefix)
    --gene <str>    gene name ($gene)
    --help          help information
EOF
exit;
}

###############################################
# parameters:
###############################################
my (@fasta, $cluster);

my %opts = ();
&GetOptions(\%opts, "fasta=s@", "clust=s", "prefix",
                    "gene=s", "help");

&prompt if exists $opts{help} or !%opts;

if (!exists $opts{fasta} or !exists $opts{clust}) {
	print STDERR "Both --fasta and --clust are required\n";
	&prompt;
} else {
	@fasta = @{$opts{fasta}} if exists $opts{fasta};
	$cluster = $opts{clust} if exists $opts{clust};
}

$gene = $opts{gene} if exists $opts{gene};
$prefix = $opts{prefix} if exists $opts{prefix};

###############################################
# fasta sequence
###############################################
my (%fasta_seq, $seq, $name);
foreach my $efasta (@fasta) {
	open(IN, $efasta) || die;
	while (<IN>) {
		chomp;
		if (/^>(.+)/) {
			if (defined $name) {
				$fasta_seq{$name} = $seq;
			}	
			$name = $1;
			$seq = '';
		} else {
			$seq .= $_;
		}
	}
	$fasta_seq{$name} = $seq;
	close IN;
}

###############################################
# output files
###############################################
#my $ncluster_out = $prefix.".cdhit.ncluster";
#my $seqlist_out = $prefix."cdhit.seqlist";
#my $seq_out = $prefix.".cdhit.fasta";

#open(NCLUSTER, ">", $ncluster_out) || die;
#open(SEQLIST, ">", $seqlist_out) || die;
#open(SEQ, ">", $seq_out) || die;

###############################################
# cd-hit cluster
###############################################
#>Cluster 0
#0	60001nt, >Ki3... *
#>Cluster 1
#0	30804nt, >CML103... at +/89.42%
#1	30891nt, >CML277... at +/96.62%

my $ncluster = 0;
my (%cluster, $cluster_id, $seqname);
open(CLUSTER, $cluster) || die;
while (<CLUSTER>) {
	chomp;
	if (/^>Cluster (\d+)/) {
		$ncluster++;
		$cluster_id = $1;
	} elsif (/>([^\.]+)\./) {
		$seqname = $1;
		if (exists $fasta_seq{$seqname}) {
			print ">$seqname\#cluster$cluster_id\#$gene\n";
			#print "$fasta_seq{$seqname}\n";
			&format_print($fasta_seq{$seqname}, 80);
		} else {
			print STDERR "$seqname does not exist in fasta input\n";
		}
	} else {
		print STDERR "Unexpected line(s) in $cluster: \n";
		print STDERR "    $_\n";
	}
}
close CLUSTER;

#print NCLUSTER "$ncluster\n";
#close NCLUSTER;
#close SEQLIST;
#close SEQ;

###############################################
### function for formatted output:
###############################################
sub format_print {
	my ($inseq, $formatlen) = @_;
	while (my $chunk = substr($inseq, 0, $formatlen, "")) {
		print "$chunk\n";
	}
}

