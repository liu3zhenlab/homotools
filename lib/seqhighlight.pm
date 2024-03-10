#!/usr/bin/perl -w
#Author: Sanzhen Liu
#Date: 5/6/2021
#
package seqhighlight;

sub seqhighlight {
# print a genomic sequence and highlight cdna region
# cdna with capitalized; others in little cases
# start and stop codons are marked
	my ($infas, $inbed, $outfas, $transcript_name) = @_;
	my $inseq_name;
	my $inseq = "";
	
	# fasta (only one sequence)
	open(INFAS, "<", $infas) || die;
	while (<INFAS>) {
		chomp;
		if (/^>(\S+)/) {
			$inseq_name = $1;
		} else {
			$inseq .= $_;
		}
	}
	
	# bed file
	my (%exon, %cds);
	open(BED, "<", $inbed) || die;
	while (<BED>) {
		chomp;
		if (!/^#/) {
			my @line = split(/\t/, $_);
			my $start = $line[1];
			my $end = $line[2];
			if ($line[3] eq "exon") {
				$exon{$start} = $end;	
			} elsif ($line[3] eq "CDS") {
				$cds{$start} = $end;
			}
		}
	}
	close INFAS;

	# partition
	my @exon_starts = sort {$a <=> $b} keys %exon;
	my @cds_starts = sort {$a <=> $b} keys %cds; 

	my $cds_start_point = $cds_starts[0]; # 0-based start point
	my $exon_start_point = $exon_starts[0]; # 0-based start point
	
	my $cds_end_point = $cds{$cds_starts[$#cds_starts]}; # 1-based end point
	my $exon_end_point = $exon{$exon_starts[$#cds_starts]}; # 1-based end point

	my (%utr5, %utr3);
	
	# 5' upstream
	my $up = $exon_start_point;
	
	# 5' UTR
	foreach my $estart (@exon_starts) {
		if ($estart < $cds_start_point) {
			if ($cds_start_point < $exon{$estart}) {		
				$utr5{$estart} = $cds_start_point;
			} else {
				$utr5{$estart} = $exon{$estart};
			}
		}
	}

	# 3' UTR
	foreach my $estart (@exon_starts) {
		if ($exon{$estart} > $cds_end_point) {
			if ($estart < $cds_end_point) {
				$utr3{$cds_end_point} = $exon{$estart};
			} else {
				$utr3{$estart} = $exon{$estart};
			}
		}
	}

	# 3' downstream
	my $down = $exon_end_point;
	
	# output
	my (%finalseq);
	my $utr5seq = "";
	my $utr3seq = "";
	# upstream seq
	if ($up > 0) {
		my $upseq = &coordinates2str($inseq, 0, $up, "");
		$finalseq{0} = lc $upseq;
	}
	
	# 5' UTR
	if (%utr5) { # not empty
		foreach my $utr5start (sort {$a <=> $b} keys %utr5) {
			my $each_utr5seq = &coordinates2str($inseq, $utr5start, $utr5{$utr5start}, "");
			$finalseq{$utr5start} = uc $each_utr5seq;
		}
	}
	
	# all introns
	for (my $i=0; $i<$#exon_starts; $i++){
		my $intron_start = $exon{$exon_starts[$i]};
		my $intron_end = $exon_starts[$i+1];
		my $intron_seq = &coordinates2str($inseq, $intron_start, $intron_end, "");
		$finalseq{$intron_start} = lc $intron_seq;
	}

	# all cds
	foreach my $cds_start (@cds_starts) {
		my $cds_end = $cds{$cds_start};
		my $cds_seq = &coordinates2str($inseq, $cds_start, $cds_end, "**");
		$finalseq{$cds_start} = uc $cds_seq;
	}
	
	# 3' UTR
	if (%utr3) {
		foreach my $utr3start (sort {$a <=> $b} keys %utr3) {
			my $each_utr3seq = &coordinates2str($inseq, $utr3start, $utr3{$utr3start}, "");
			$finalseq{$utr3start} = uc $each_utr3seq;
		}
	}

	# 3' downstream seq
	if ($down < length($inseq)) {
		my $downseq = &coordinates2str($inseq, $down, length($inseq), "");
		$finalseq{$down} = lc $downseq;
	}

	open(OUTFAS, ">", $outfas) || die;
	# join all strings
	my $finalseq = "";
	foreach (sort {$a <=> $b} keys %finalseq) {
		$finalseq .= $finalseq{$_};
	}
	print OUTFAS ">$transcript_name\n";
	#while (my $chunk = substr($finalseq, 0, 80, "")) {
	#	print OUTFAS "$chunk\n";
	#}
	print OUTFAS "$finalseq\n";
	close OUTFAS;
	
	
	# module
	# extract sequences based on coordinates:
	sub coordinates2str {
		my $extracted_seq = "";
		my ($in_seq, $in_start, $in_end, $flank_character) = @_;
		# $in_start: 0-based
		# $in_end: 1-based
		
		my $in_seq_len = length($in_seq);
		my $extract_length = $in_end - $in_start;
		if ($in_end <= $in_seq_len) {
			$extracted_seq = substr($in_seq, $in_start, $extract_length);
			$extracted_seq = $flank_character.$extracted_seq.$flank_character;
		} elsif ($in_start < $in_seq_len) {
			$in_end = $in_seq_len;
			$extracted_seq = substr($in_seq, $in_start, $extract_length);
		}
		return $extracted_seq;
	}
}

1;
