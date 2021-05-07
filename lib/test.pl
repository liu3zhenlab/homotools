my $files=`ls -1`;
@files = split(/\n/, $files);
foreach (@files) {
	print "$_--\n";
}

