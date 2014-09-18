#!/usr/bin/perl

package PUGREST;
use LWP::Simple;
use IDMAP;
use Data::Dumper;
use DrugTargetBase;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(get_InChIKey_by_name get_InChIKey_by_compound get_CID_by_substance get_InChIKey_by_CID);

#--------------------------TEST AREA------------------------

#--------------------------TEST AREA------------------------

sub get_InChIKey_by_name
{
	#input PubChem compound or substance name then output an array of InChIKeys if exist
	my $name = shift @_;
	chomp($name);
	my $manualInChIKey = manual_get_InChIKey_by_chemicalID($name);
	return $manualInChIKey if ($manualInChIKey);	#stop if found in manual map--save time for accessing PubChem

	my @inchikeys = get_InChIKey_by_compound($name);
	my @cids = get_CID_by_substance($name);
	foreach my $cid ( @cids ){
		my $inchikey = get_InChIKey_by_CID($cid);
		push @inchikeys, $inchikey;
	}
	my @unique_ikeys = unique(\@inchikeys);
	if (@unique_ikeys){
		return @unique_ikeys;
	} else {
		return;
	}
}

sub get_InChIKey_by_compound
{
	#input PubChem compound name then output an array of InChIKeys if exist
	my $compound = shift @_;
	my $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$compound/property/InChIKey/TXT";
	my $inchikey = get $url;
	unless($inchikey){
		$compound =~ s/_/ /g;	#change _ to a space
		$url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$compound/property/InChIKey/TXT";
		$inchikey = get $url;
	}
	return unless $inchikey;

	my @is = split(/\n/, $inchikey);
	my @ikeys = chomp_array(\@is);
	my @unique_ikeys = unique(\@ikeys);
	return @unique_ikeys;
}

sub get_CID_by_substance
{
	#input PubChem substance name then output an array of CIDs if exist
	my $substance = shift @_;
	my $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/$substance/cids/TXT";
	my $cid = get $url;
	unless($cid){
		$substance =~ s/_/ /g;	#change _ to a space
		$url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/$substance/cids/TXT";
		$cid = get $url;
	}
	return unless $cid;
	my @cs = split(/\n/, $cid);
	my @cids = chomp_array(\@cs);
	my @unique_cids = unique(\@cids);
	return @unique_cids;
}
sub get_InChIKey_by_CID
{
        #input PubChem CID then output InChIKey
        #one CID gives only one InChIKey
	my $cid = shift @_;
        my $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$cid/property/InChIKey/TXT";
        my $inchikey = get $url;
        chomp($inchikey);
	if (defined($inchikey)){
		return $inchikey;
	} else {
		return;
	}
}

1;
