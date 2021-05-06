#!/usr/bin/perl -w
# Author: Sanzhen Liu
# Date: 3/28/2021

package gffadjust;
# adjust position to a newly defined start for gene information in a GTF
# and output a position-adjusted GTF and BED file for each transcript
#
# input:
# 1. gff
# 2. start site
# 3. seqname
#
# output
# 1. pos-adjusted gff (one file per transcript)
# 2. pos-adjusted bed (one file per transcript)
#

sub gffadjust {
	use strict;
	use warnings;
	
	my $feature_height = 0.03;
	my $feature_color = "gray20";

	my ($outpath, $gff, $start, $seqname) = @_;
	my %entries;

	my $gff_prefix = $gff;
	$gff_prefix =~ s/.*\//$outpath\//g;
	$gff_prefix =~ s/.gff3?$//g;
	
	my $adj_gff_out = $gff_prefix.".adjusted.gff";
	open(ADJGFF, ">", $adj_gff_out) || die;
	my $adj_bed_out = $gff_prefix.".adjusted.bed";
	open(ADJBED, ">", $adj_bed_out) || die;

	my $row = 0;
	open(GFF, "<", $gff) || die;
	while (<GFF>) {
		chomp;
		$row++;
		my @line = split(/\t/, $_);
		my @adjust_line = @line;
		#my $info = $line[8];
		my $feature = $line[2];
		my $strand = $line[6];
		$adjust_line[0] = $seqname;
		$adjust_line[3] = abs($adjust_line[3] - $start) + 1;
		$adjust_line[4] = abs($adjust_line[4] - $start) + 1;

		if ($adjust_line[3] > $adjust_line[4]) {
			my $tmp_pos = $adjust_line[3];
			$adjust_line[3] = $adjust_line[4];
			$adjust_line[4] = $tmp_pos;
			if ($strand eq "+") {
				$adjust_line[6] = "-";
			} else {
				$adjust_line[6] = "+";
			}
		}
		
		my $adjust_line = join("\t", @adjust_line);
		$entries{$adjust_line[3]}{$row} = $adjust_line;
	}
	close GFF;
	
	# output entries
	foreach my $pos (sort {$a <=> $b} keys %entries) {
		my %entry_row = %{$entries{$pos}};
		foreach my $rowid (sort {$a <=> $b} keys %entry_row) {
			print ADJGFF "$entry_row{$rowid}\n";
			my $entrybed = &gff2bed($entry_row{$rowid}, $feature_height, $feature_color);
			print ADJBED "$entrybed\n";
		}
	}
	
	# close files
	close ADJGFF;
	close ADJBED;
}


###########################################
# modules
###########################################
### GFF row to BED
sub gff2bed {
# output a row of gtf data to a bed format
# chr start end name height strand color
	my ($gff_row, $height, $color) = @_;
	my @gff = split(/\t/, $gff_row);
	#10	gramene	gene	1001	3460	.
	my $bname = $gff[0];
	my $bstart = $gff[3] - 1;
	my $bend = $gff[4];
	my $group = $gff[2];
	my @bed_row = ($bname, $bstart, $bend, $group, $height, "+", $color);
	my $bed_row = join("\t", @bed_row);
	return $bed_row;
}

1;
