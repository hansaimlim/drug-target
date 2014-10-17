#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;

die "Usage: $0 <zinc id file>\n" unless @ARGV == 1;
my $idfile = shift @ARGV;
open my $ZINC, '<', $idfile or die "Could not open zinc id file $idfile: $!\n";
my %ids;
while (my $line = <$ZINC>){
	chomp($line);
	$ids{$line} = 1;
}
close $ZINC;
foreach my $id (keys %ids){
	chomp($id);
	my $url_front = "http://zinc.docking.org/results/annotation?annotation.name=".$id;
	my $url_smile = $url_front . "&annotation.type=B10&page.format=smiles&page.size=all";
	my $url_mol2 = $url_front . "&annotation.type=B10&page.format=mol2&page.size=all";
	my $smile_file = "./static/zinc/smile/$id".".smi";
	my $mol2_file = "./static/zinc/mol2/$id".".mol2";
	getstore($url_smile, $smile_file);
	getstore($url_mol2, $mol2_file);
}
exit;
