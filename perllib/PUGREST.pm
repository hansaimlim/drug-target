#!/usr/bin/perl

package PUGREST;
use LWP::Simple;
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
	my @inchikeys = get_InChIKey_by_compound($name);
	my @cids = get_CID_by_substance($name);
	foreach my $cid ( @cids ){
		my $inchikey = get_InChIKey_by_CID($cid);
		push @inchikeys, $inchikey;
	}
	my @unique_ikeys = unique(\@inchikeys);
	return @unique_ikeys;
}

sub get_InChIKey_by_compound
{
	#input PubChem compound name then output an array of InChIKeys if exist
	my $compound = shift @_;
	my $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$compound/property/InChIKey/TXT";
	my $inchikey = get $url;
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
        return $inchikey;
}

1;
