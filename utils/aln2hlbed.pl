#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
# Sanzhen Liu
# 7/30/2022

sub prompt {
	print <<EOF;
    Usage: perl $0 <Input Fasta File> [--all] [--average] [--help]
    --fmt   : input format: blast or nucmer (blast)
    --height: height for highlighted regions (0.05)
    --col   : color for highlighted regions (darkolivegreen4)
    --help  : help information
EOF
exit;
}


my ($format, $color, $height, $help);
# read the parameters:
&GetOptions("fmt=s"    => \$format,
            "col=s"    => \$color,
			"height=f" => \$height,
			"help"     => \$help) || &prompt;

&prompt if $help || ! defined $ARGV[0];

if (!defined $format) {
	$format = "blast";
} elsif (($format ne "blast") or ($format ne "nucmer")) {
	print STDERR "--fmt must be \"blast\" or \"nucmer\"\n";
	exit;
}

$color="darkolivegreen4" if !defined $color;

if (!defined $height) {
	$height = 0.05
} elsif ($height < 0.01 or $height > 0.1) {
	print STDERR "--height must be from 0.01 to 0.1\n";
	exit;
}

# format conversion

#NODE_656_length_10962_cov_26.916766	B71_pwl2	100.00	438	0	0	5883	6320	1	438	0.0	  809
#chr start(0-based) end(1-based) label height(0.01-0.1) strand(+/-) color(R compatible)

open(IN, "<", $ARGV[0]) || die;
while (<IN>) {
	chomp;
	if ($format eq "blast") {
		&blast2out($_, $height, $color);
	} elsif ($format eq "nucmer") {
		#&nucmer2out;
	}
}
close IN;

sub blast2out {
	my ($in_row, $in_height, $in_color) = @_;
	my @line = split(/\t/, $in_row);
	if ($#line > 10) {
		my $seq = $line[0];
		my $start = $line[6] - 1;
		my $end = $line[7];
		my $strand = "+";
		if ($line[9] < $line[8]) {
			$strand = "-";
		}
		my $label=$line[1];
		print "$seq\t$start\t$end\t$label\t$in_height\t$strand\t$in_color\n";
	}
}

