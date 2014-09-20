#!/usr/bin/perl
use strict;
use warnings;
use IDMAP;
my $is_demo_on = 0;

my $file = "./static/String/9606.protein.links.v9.1-GS-dist.txt";
$file = "./static/String/9606.protein.links.v9.1-GS-dist_demo.txt" if $is_demo_on == 1;

my $outfile = "./static/String/9606.protein.links.v9.1-GN-dist.txt";
$outfile = "./static/String/9606.protein.links.v9.1-GN-dist_demo.txt" if $is_demo_on == 1;

my $ug_ref = UniProtKB_genename();	#a hash reference containing UniProtKB => genename conversions
open my $FI, '<', $file or die "Could not open string file $file: $!\n";
open my $FO, '>', $outfile or die "Could not open output file $outfile: $!\n";
print $FO "GeneName1\tGeneName2\tDistance\n";
while(my $line = <$FI>){
	next if $. == 1;
	my @words = split(/\t/, $line);
	my $uniprot1 = shift @words;
	my $uniprot2 = shift @words;
	my $distance = shift @words;
	chomp($distance);

	my $gene1;
	$gene1 = $ug_ref->{$uniprot1} or $gene1 = 0;
	chomp($gene1);
	my $gene2;
	$gene2 = $ug_ref->{$uniprot2} or $gene2 = 0;
	chomp($gene2);

	if ($gene1){
		if ($gene2){	#both conversions exist
			print $FO $gene1, "\t", $gene2, "\t", $distance, "\n";
		} else {	#only first gene converted
			print $FO $gene1, "\t", $uniprot2, "\t", $distance, "\n";
		}
	} else {
		if ($gene2){	#only second gene converted
			print $FO $uniprot1, "\t", $gene2, "\t", $distance, "\n";
		} else{		#no conversion
			print $FO $line;	#no new line needed
		}
	}
}
close $FO;
close $FI;
exit;
