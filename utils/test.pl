my $string = "ANN=C|||missense_variant&splice_region_variant|MODERATE|EXON_Zm00001eb320870_1001_1200|Zm00001eb320870|transcript|Zm00001eb320870_T002|protein_coding|4/8|c.172G>C|p.Val58Leu|380/942|172/360|58/119|||||";
my @string = split(/\|/, $string);
print "$string\n";
print join(";", @string);
print "\n";
