#!/usr/bin/perl -w
# Author: Sanzhen Liu
# Date: 3/28/2021

package genegtfadjust;
#
# adjust position to a newly defined start for gene information in a GTF
# and output a position-adjusted GTF and BED file for each transcript
#
# input:
# 1. gtf
# 2. start site
# 3. outdir
# 4. seqname
#
# output
# 1. pos-adjusted gtf (one file per transcript)
# 2. pos-adjusted bed (one file per transcript)
#

sub genegtfadjust {
	use strict;
	use warnings;
	my %transcripts;	
	my ($genegtf, $start, $outdir, $seqname) = @_;
	open(GENEGTF, "<", $genegtf) || die;
	while (<GENEGTF>) {
		chomp;
		my @line = split(/\t/, $_);
		my @adjust_line = @line;
		my $info = $line[8];
		my $feature = $line[2];
		my $strand = $line[6];
		
		if ($info =/transcript_id \"(.+?)\"/) {
			my $transcript_name = $1;
			if ($feature eq "exon" or $feature eq "CDS") {
				my @exon_CDS_info = @adjust_line;
				$exon_CDS_info[0] = $seqname;
				$adjust_line[3] = abs($adjust_line[3] - $start) + 1;
				$adjust_line[4] = abs($adjust_line[4] - $start) + 1;
				@exon_CDS_info[3..4] = sort {$a <=> $b} @adjust_line[3..4];
				$exon_CDS_info[6] = "+";
				my $exon_CDS_info = join("\t", @exon_CDS_info);
				$transcripts{$transcript_name}{$feature}{$exon_CDS_info[3]} = $exon_CDS_info;
			}
		}
	}
	close GENEGTF;


	###########################################
	### output transcripts
	###########################################
	my $gene_height = 0.015;
	my $gene_color = "gray60";
	my $exon_height = 0.05;
	my $exon_color = "steelblue1";
	my $cds_height = 0.05;
	my $cds_color = "steelblue3";
	foreach my $transcript (keys %transcripts) {
		my $adj_gtf_out = $outdir."/".$transcript.".adjusted.gtf";
		my $adj_bed_out = $outdir."/".$transcript.".adjusted.bed";

		open(ADJGTF, ">", $adj_gtf_out) || die;
		open(ADJBED, ">", $adj_bed_out) || die;

		my %transcript_exon_CDS_info = %{$transcripts{$transcript}};
	
		# exon
		if (exists $transcript_exon_CDS_info{exon}) {
			my %pos_info = %{$transcript_exon_CDS_info{exon}};
			foreach my $pos (sort {$a <=> $b} keys %pos_info) {
				print ADJGTF "$pos_info{$pos}\n";
				my $exonbed = &gtf2bed($pos_info{$pos}, $exon_height, $exon_color);
				print ADJBED "$exonbed\n";
			}
		}

		# CDS
		if (exists $transcript_exon_CDS_info{CDS}) {
			my %pos_info = %{$transcript_exon_CDS_info{CDS}};
			foreach my $pos (sort {$a <=> $b} keys %pos_info) {
				print ADJGTF "$pos_info{$pos}\n";
				my $cdsbed = &gtf2bed($pos_info{$pos}, $cds_height, $cds_color);
				print ADJBED "$cdsbed\n";
			}
		}
	
		# close files
		close ADJGTF;
		close ADJBED;
	}
}


###########################################
# modules
###########################################
### GTF row to BED
sub gtf2bed {
# output a row of gtf data to a bed format
# chr start end name height strand color
	my ($gtf_row, $height, $color) = @_;
	my @gtf = split(/\t/, $gtf_row);
	#10	gramene	gene	1001	3460	.
	my $bname = $gtf[0];
	my $bstart = $gtf[3] - 1;
	my $bend = $gtf[4];
	my $group = $gtf[2];
	my @bed_row = ($bname, $bstart, $bend, $group, $height, "+", $color);
	my $bed_row = join("\t", @bed_row);
	return $bed_row;
}

1;
