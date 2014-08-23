#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;
use List::MoreUtils qw/ uniq /;
use Data::Dumper;

#----------------------------------------------------------------------------
# Purpose	: To convert list of Drugnames to InChI keys
# Input		: drug list file (one column)
# Output	: drugname-InChI Key; two-column (or more) file (tab separated)
# Author	: Hansaim Lim
# Date		: 22 Aug, 2014
#----------------------------------------------------------------------------

die "Usage: $0 <drug list file>\n" unless @ARGV == 1;
my $druglist = shift @ARGV;
my $outfile = $druglist;
$outfile =~ s/^(.+)(\..+)$/$1InChIKey$2/;

my $Prestwick = "./cMap_Prestwick_map.tsv";	#Prestwick ID with chemical name map
open my $PRESTWICK, '<', $Prestwick or die "Could not open file $Prestwick: $!\n";
my %prestwick_chemicalname = ();
while (<$PRESTWICK>){
	my @ids = split(/\t/, $_);
	my $prestwick = shift @ids;
	my $chemicalname = shift @ids;
	chomp($chemicalname);
	$prestwick_chemicalname{$prestwick} = $chemicalname;
}
close $PRESTWICK;

my $drug_inchikey_manualmap = "./cMap_drug_inchikeys_manualmap.tsv";
open my $MANUALMAP, '<', $drug_inchikey_manualmap or die "Could not open file $drug_inchikey_manualmap: $!\n";
my %manual_drug_inchikey = ();
while (<$MANUALMAP>){
	my @words = split(/\t/, $_);
	my $drug = shift @words;
	my $inchikey = shift @words;
	chomp($inchikey);
	$manual_drug_inchikey{$drug} = $inchikey;
}
close $MANUALMAP;

my $dn_col = 0;	#column index for drugname.
my %dn_ikeys = ();	#contains drugname-InChI key pair
my $line = 1;
open my $DRUGLIST, '<', $druglist or die "Could not open file $druglist: $!";
LINE: while (<$DRUGLIST>) {
	my $drug = $_;
	chomp($drug);
	if ( $line == 1) {
		$line++;
		next LINE;	#skip the first line
	}
	
InChIKey: {
	if ($drug =~ m/^Prestwick/i){	#this drug is named by Prestwick ID
		$drug = $prestwick_chemicalname{$drug};
	}
	my $inchikey = pubchem_inchikey_by_drug($drug);
	if ( not defined $inchikey ) {
		my @cids = pubchem_cids_by_substance($drug);
		unless (@cids){	#the drug name is not searchable in pubchem
			$dn_ikeys{$drug}{0} = $manual_drug_inchikey{$drug};
			next LINE;
		}
		my @unique_cids = uniq @cids;
		my @inchikeys = ();
		my $i = 0;
		foreach my $cid ( @unique_cids ){
			my $inchi = pubchem_inchikey_by_cid($cid);
			next unless defined $inchi;
			$dn_ikeys{$drug}{$i} = $inchi;
			$i++;
		}
	} else {
		my @ikeys = split(/\n/, $inchikey);
		my $i = 0;
		foreach my $ikey ( @ikeys ) {
			chomp($ikey); 
			$dn_ikeys{$drug}{$i} = $ikey;
			$i++;
		}
	 }
}
}
close $DRUGLIST;

open my $OUTPUT, '>', $outfile or die "Could not open file $outfile: $!\n";
foreach my $drug ( keys %dn_ikeys ) {
	print $OUTPUT $drug;
	foreach my $i ( keys $dn_ikeys{$drug} ) {
		print $OUTPUT "\t",$dn_ikeys{$drug}{$i};
	}
	print $OUTPUT "\n";
}
close $OUTPUT;

sub pubchem_inchikey_by_drug {
	my $drug = shift @_;
	my $url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$drug/property/InChIKey/TXT";
	my $inchikey = get $url;
	chomp($inchikey);
	return $inchikey;
}
sub pubchem_cids_by_substance {
	my $substance = shift @_;
	my $url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/$substance/cids/TXT";
	my $cid = get $url;
	chomp($cid);
	my @cids = split(/\n/, $cid);
	return @cids;
}
sub pubchem_inchikey_by_cid {
	my $cid = shift @_;
	my $url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$cid/property/InChIKey/TXT";
	my $inchikey = get $url;
	chomp($inchikey);
	return $inchikey;
}
exit;
