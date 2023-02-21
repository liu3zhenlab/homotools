#!/usr/bin/perl
# Sanzhen Liu
# 7/6/2019

my $current_version = "v0.1";

use strict;
use warnings;
use Getopt::Long;

# fasta 
# show-snps -T output (test on mummer: 4.0.0beta2)

sub errINF {
	print <<EOF;
	Usage: perl showsnps2vcf.pl --fasta <fasta> --var <show-snps output> [options]
	- to parse show-snps output to vcf format

	Options:
	--fasta|f <fasta>:             fasta sequence file; required
	--refcol|r <num>:              column number for reference sequence names; default=11
	--var|v <show-snps -T output>: show-snps -T output; required
	--geno|g <geno name>:          genotype name; use -var file name if not specified.
	--overlapAllow:                allow multiple variant types on a position if specified.
	--version:                     version information
	--help|h:                      help information
EOF
	exit;
}

my ($fasta, $var, $refcol, $geno, $overlap_allow, $version, $help);
GetOptions("fasta|f=s" => \$fasta,
           "var|v=s" => \$var,
		   "refcol|r=i" => \$refcol,
		   "geno|g=s" => \$geno,
		   "overlapAllow|o" => \$overlap_allow,
		   "version" => \$version,
		   "help|h" => \$help);

# help info:
&errINF if $help;

if ($version) {
	print "$current_version\n";
	exit;
}

if (!defined $geno) {
	$geno = $var;
	$geno =~ s/.*\///g;
}

&errINF if (!defined $fasta or !defined $var);
$refcol = 11 if (!defined $refcol);

#fasta
print STDERR "====\n";
print STDERR "To convert show-snps result to vcf format\n";
print STDERR "o input fasta file: $fasta\n";
print STDERR "o input show-snps output file: $var\n";
&timereport;
print STDERR "Start to read fasta\n";
my (%chr_seq, $seq, $chr);
open(IN, $fasta) || die;
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		if (defined $seq) { 
			$chr_seq{$chr} = $seq;
		}       
		$chr = $1;
		$seq = "";
	} else {
		$seq .= $_; 
	}
}
$chr_seq{$chr} = $seq;
close IN;
&timereport;
print STDERR "Finish reading fasta\n";

#show-snps -T output file
my $nlines = &read_nrows($var);
print STDERR "$nlines lines in $var\n";
my $count = 0;
my $report_point = 10;

open(IN, $var) || die;
while (<IN>) {
	$count++;
	last if /^\[/; # skip heads
}

my (%var, %del, %check, $delpos);
while (<IN>) {
	chomp;
	my @line = split;
	my $ref = $line[$refcol - 1];
	if (!exists $chr_seq{$ref}) {
		print STDERR "$ref does not exist in fasta\n";
		exit;
	}

	my ($pos, $rseq, $qseq) = @line[0..2];
	if ($rseq ne "." and $qseq ne ".") { # snp
		$var{$ref}{$pos}{R} = $rseq; # ref snp
		$var{$ref}{$pos}{Q} = $qseq; # query snp
		$check{$ref.":".$pos}{"snp"} = 1;
	} elsif ($rseq eq "." and $qseq ne ".") { # ins
		if (exists $var{$ref}{$pos}{Q}) { 
			my $querybase = $var{$ref}{$pos}{Q};
			$querybase .= $qseq;
			$var{$ref}{$pos}{Q} = $querybase;
		} else {
			my $refbase = substr($chr_seq{$ref}, $pos - 1, 1);
			$var{$ref}{$pos}{R} = $refbase;
			$var{$ref}{$pos}{Q} = $refbase.$qseq;
			$check{$ref.":".$pos}{"ins"} = 1;
		}
	} elsif ($rseq ne "." and $qseq eq ".") { # del
		if (exists $del{$ref}{$pos - 1}) { # continous deletion?
			my $delbase = $var{$ref}{$delpos}{R};
			$delbase .= $rseq;
			$var{$ref}{$delpos}{R} = $delbase;
			$del{$ref}{$pos} = 1; # record the last del position
		} else {
			$delpos = $pos - 1; # position before del
			my $refbase = substr($chr_seq{$ref}, $delpos - 1, 1); # base before del
			$var{$ref}{$delpos}{R} = $refbase.$rseq;
			$var{$ref}{$delpos}{Q} = $refbase;
			$del{$ref}{$pos} = 1;
			$check{$ref.":".$pos}{"del"} = 1;
		}
	}
	
	# status report:
	if ($count / $nlines * 100 > $report_point) {
		print STDERR "$report_point% have been read\n";
		$report_point += 10;
	}
	$count++;
}
close IN;

# status report
print STDERR "$report_point% have been read\n";
&timereport;

# check
my $pos_overlap_count = 0;
foreach (keys %check) {
	my %erefpos = %{$check{$_}};
	my @types = keys %erefpos;
	if ($#types > 0) {
		my $type = join(",", @types);
		print STDERR "$_ has multiple variant types: $type.\n";
		$pos_overlap_count++;
	}
}

# if overlaps found, report:
if ($pos_overlap_count > 0) {
	print STDERR "$pos_overlap_count positions have multiple variant types (e.g., snps, indel)\n";
	if ($overlap_allow) {
		print STDERR "The result is output because this is allowed\n";
	} else {
		exit;
	}
}

# status report:
if ($pos_overlap_count == 0) {
	print STDERR "Check completed. Hooray! No positions with multiple variant types (e.g., snps, indel)\n";
} else {
	print STDERR "Check completed. $pos_overlap_count positions have multiple variant types (e.g., snps, indel)\n";
}

# output
print "##fileformat=VCFv4.0\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$geno\n";

foreach my $eseq (sort {$a cmp $b} keys %var) { 
	my %allpos = %{$var{$eseq}};
	foreach my $epos (sort {$a <=> $b} keys %allpos) {
		print "$eseq\t$epos\t.\t$allpos{$epos}{R}\t$allpos{$epos}{Q}\t.\t.\t.\tGT\t1\|1\n";
	}
}

# status report:
&timereport;
print STDERR "All completed!\n";
print STDERR "====\n";

######################
# modules
sub timereport {
# print current time
	my $curtime = localtime();
	print STDERR "$curtime\n";
}

sub read_nrows {
# determine number of rows in the input file
	my $file_to_check = shift;
	my $wc_cmd = sprintf("%s%s", "wc -l ", $file_to_check);
	chomp(my $num_rows = `$wc_cmd`);
	$num_rows =~ s/ .*//g;
	return($num_rows);
}

