#!/usr/bin/perl -w
# Author: Sanzhen Liu
# Date: 3/28/2021

package maxbed;

# select a bed from one or more BED files
# compare a set of BED files and select the one with the maximal total length
#
# input:
# 1. bed files

# output
# 1. the filename with the path for the selected BED file

sub maxbed {
	my $max_total_len = 0;
	my $select_bed;
	my $bedlist = shift;
	my @bed = split(/[ \n]/, $bedlist);
	foreach my $bed (@bed) {
		#chomp $bed;
		my $total_len = 0;
		open(BED, "<", $bed) || die;
		while (<BED>) {
			chomp;
			my @line = split("\t", $_);
			$total_len += $line[2] - $line[1];
		}
		close BED;

		if ($total_len > $max_total_len) {
			$max_total_len = $total_len;
			$select_bed = $bed;
		}
	}
	return $select_bed;
}

1;

