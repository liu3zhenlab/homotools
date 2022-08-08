#!/usr/bin/perl -w
#======================================================================
# homostack
#
# Author: Sanzhen Liu <liu3zhen@ksu.edu>
# 8/6/2022
#
# to make input sequences the same strand and adjust corresponding BED 
#======================================================================

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $version = "0.1.0";

my $reverse_complement_text = "_rc";

my $identity = 90; 
my $match = 100;
my $prefix = "seqbed";
my $threads = 1;

sub prompt {
    print <<EOF;
    Usage: perl $0 --seq <fasta> --bed <bed_file> [options]
    [Options]
    --seq <file>     fasta file containing a sequence as the query; required
                     multiple sequences can be input by using --seq multiple times
                     the first sequence will be used as the reference sequence
    --bed <file>   bed file; if specified, the same number of --bed needed to be specified as --seq;
                     --seq and --bed are paired by their order, i.e., 1st --seq paired with 1st --bed;
                     if some --bed has no data, input "none".
    --identity <num> minimal percentage of identity from 0 to 100 ($identity)
    --match <num>    minimal bp match of an alignment ($match)
    --prefix <str>   the output directory and the prefix for output files ($prefix)
    --threads <num>  number of cpus ($threads)
    --version        version information
    --help           help information.
EOF
exit;
}

###############################################
# parameters:
###############################################
my %opts = ();
my ($query, $db);

&GetOptions(\%opts, "seq=s@", "bed=s@",
                    "identity=i", "match=i",
                    "prefix=s", "threads=i",
                    "version", "help");

if (exists $opts{version}) {
	print "$version\n";
	exit;
}

&prompt if exists $opts{help} or !%opts;
my (@seq, @bed);

if (!exists $opts{seq}) {
	print STDERR RED, "--seq is required\n", RESET;
	&prompt;
} else {
	@seq = @{$opts{seq}} if exists $opts{seq};
	if (exists $opts{bed}) {
		@bed = @{$opts{bed}};
		if ($#seq != $#bed) {
			print STDERR RED, "Numbers of --seq and --bed MUST be equal\n", RESET;
			&prompt;
		}
	}
}

### check if input files exist
my $error_files_exist = 0;
foreach (@seq) {
	my $check_error = &file_not_exist($_);
	$error_files_exist += $check_error;
}

if (exists $opts{bed}) {
	foreach (@bed) {
		if ($_ ne "none") {
			my $check_error = &file_not_exist($_);
			$error_files_exist += $check_error;
		}
	}
}
exit if ($error_files_exist); # quit if any file is missing

### other parameters:
$identity = $opts{identity} if exists $opts{identity};
$match = $opts{match} if exists $opts{match};
$prefix = $opts{prefix} if exists $opts{prefix};
$threads = $opts{threads} if exists $opts{threads};

###############################################
# preparation
###############################################
# create a directory for outputs
if (-d $prefix) {
	print STDERR RED, "Warning: the directory $prefix exists.\n", RESET;
} else {
	`mkdir $prefix`;
}

# script path:
# my $scriptPath = $FindBin::Bin;
# my $utilsPath = $scriptPath."/utils/";

###############################################
# check requirments
###############################################
&cmd_check("nucmer");
&cmd_check("realpath");
&cmd_check("show-coords");

&runreport("step 01: Requirement checked");

###############################################
# intermediate output
###############################################
my $nucmer_prefix = $prefix."/".$prefix.".2.nucmer.";

###############################################
# compare via nucmer
###############################################

# first sequence and BED:
my $ref_seq = &re2abs($seq[0]);
`ln -s $ref_seq $prefix/`;
if (exists $opts{bed}) {
	if (-f $bed[0]) {
		`cp $bed[0] $prefix/ `;
	}
}
# other sequences and BED:
for (my $i=1; $i<=$#seq; $i++) {
	my $current_seq = &re2abs($seq[$i]);
	my $current_bed;
	if (exists $opts{bed}) {
		$current_bed = &re2abs($bed[$i]);
	}
	# nucmer alignment
	my $delta_out = $nucmer_prefix.$i.".delta";
	my $nucmer_show = $nucmer_prefix.$i.".txt";
	`nucmer --mum --delta $delta_out --threads $threads $ref_seq $current_seq`;
	`show-coords -HTl -I $identity -L $match $delta_out > $nucmer_show`;
	
	# assess alignment output
	#3255    10962   1       7639    7708    7639    99.01   10962   7639    ref	subject
	my $is_rc_needed = 0;
	if (-s $nucmer_show) { # non-zero
		$is_rc_needed = &reverse_gt_forward($nucmer_show);	
	} else {
		print STDERR "WARNING: no alignment output between $ref_seq and $current_seq\n";
	}
	
	if ($is_rc_needed) {
		my $seq_len = &reverse_complement($current_seq, $prefix);
		if ((exists $opts{bed}) and ($bed[$i] ne "none")) {
			&bed_complement($current_bed, $seq_len, $prefix);
		}
	} else {
		`ln -s $current_seq $prefix/`;
	}

	# cleanup
	`rm $delta_out`;
	`rm $nucmer_show`;

}
#close OUT;

&runreport("step 02: sequences aligned");

###############################################
# module 1: check command availability
###############################################
sub cmd_check {
	my $cmd = shift;
	my $cmdPath=`which $cmd 2>/dev/null`;
	if (!$cmdPath) {# not founded
		print STDERR RED, "  $cmd is not found\n", RESET;
		print STDERR RED, "Exit\n", RESET;
		exit;
	}
}

###############################################
# module 2: report running result
###############################################
# funtion to report running return
sub runreport {
	my $injob = shift;
    my $dateinfo = `date +'o %Y-%m-%d %H:%M:%S'`;
	print STDERR MAGENTA, "$dateinfo", RESET;
	print STDERR "  $injob.\n";
}

###############################################
## module: convert relative path to absolute path
###############################################
sub re2abs {
	my $inpath = shift;
	my $outpath = `realpath $inpath`;
	chomp $outpath;
	return $outpath;
}

############################################### 
## module for reverse complement sequences
###############################################
sub reverse_complement {
	my $fas_len = 0;
	my ($in_seqfile, $outdir) = @_; 
	my $original_filename = $in_seqfile;
	$original_filename =~ s/.*\///;
	open(INFAS, "<", $in_seqfile) || die;
	my $outfas_file = $outdir."/".$original_filename;
	open(OUTFAS, ">", $outfas_file) || die;
	while (<INFAS>) {
		chomp;
		if (/^>/) {
			print OUTFAS "$_";
			print OUTFAS "$reverse_complement_text\n";
		} else {
			my $revcomseq = reverse($_);
			$revcomseq =~ tr/AGCTagct/TCGAtcga/;
			$fas_len += length($_);
			print OUTFAS "$revcomseq\n";
		}
	}
	close INFAS;
	close OUTFAS;
	return $fas_len;
}

############################################### 
# module to check if a file exists
############################################### 
sub file_not_exist {
	my $file_not_exist = 0;
	my $file_to_check = shift @_;
	if ( !-f $file_to_check) {
		$file_not_exist = 1;
		print STDERR "ERROR: $file_to_check does not exist\n";
	}
	return $file_not_exist;
}

###############################################
# module to compare forward and reverse alignments
###############################################
#3255    10962   1       7639    7708    7639    99.01   10962   7639    A_12.0.321      O_ARG-60
sub reverse_gt_forward {
	my $nucmer_outfile = shift @_;
	my $is_reverse_gt_forward = 0;
	my $forward_sum_bp = 0;
	my $reverse_sum_bp = 0;
	open(NUCALN, "<", $nucmer_outfile) || die;
	while (<NUCALN>) {
		chomp;
		my @line = split("\t", $_);
		my $aln_len = $line[3] - $line[2];
		if ($aln_len > 0) { # forward
			$forward_sum_bp += $aln_len + 1;
		} else {
			$reverse_sum_bp += abs($aln_len) + 1;
		}
	}
	close NUCALN;

	if ($reverse_sum_bp > $forward_sum_bp) {
		$is_reverse_gt_forward = 1;	
	}
	return $is_reverse_gt_forward;
}

###############################################
# complement BED coordinates
###############################################
sub bed_complement {
	my ($inbed_file, $infas_len, $in_dir) = @_;
	my $bed_filename = $inbed_file;
	$bed_filename =~ s/.*\///;
	my $outbed_file = $in_dir."/".$bed_filename;
	
	open(INBED, "<", $inbed_file) || die;
	open(OUTBED, ">", $outbed_file) || die;
	while(<INBED>) {
		chomp;
		if (/^\#/) {
			print OUTBED  "$_\n";
		} else {
			my @line = split(/\t/, $_);
			my $seqname = $line[0].$reverse_complement_text;
			$line[0] = $seqname;
			my $end = $infas_len - $line[1];
			my $start = $infas_len - $line[2];
			$line[1] = $start;
			$line[2] = $end;
			if ($#line>=5) {
				my $strand = $line[5];
				if ($strand eq "+") {
					$strand = "_";
				} elsif ($strand eq "-") {
					$strand = "+";
				}
				$line[5] = $strand;
			}		
			my $newline = join("\t", @line);
			print OUTBED "$newline\n";
		}
	}
	close INBED;
	close OUTBED;
}
