#!/usr/bin/perl
use strict;
use warnings;

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose	: To extract drug-targets (genesymbols) for common drugs (drugs appearing in both (DrugBank+STITCH) and cMap)
#		  Also prints out the names of common drugs in cMap for destination finding
# Input		: DrugBank (name-target symbol; name-inchikey), STITCH (inchikey-genesymbol-score), cMap (drugname-inchikey)
# Output	: STDB sources (genesymbols in STITCH and/or DrugBank); cMap drug names (shared only)
# Author	: Hansaim Lim
# Date		: 23 Aug, 2014
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
die "Usage: $0 <DrugBank name-targetsymbol file> <DrugBank name-inchikey file> <STITCH inchikey-genesymbol-score file> <cMap drugname-inchikey file>\n" unless @ARGV == 4;
my $DrugBank_name_target_file = shift @ARGV;
my $DrugBank_name_inchikey_file = shift @ARGV;
my $STITCH_inchikey_symbol_file = shift @ARGV;
my $cMap_name_inchikey_file = shift @ARGV;

my $STDB_source_file = "./STDB_source_genesymbol.txt";
my $cMap_shared_drug_file = "./cMap_common_drugnames.txt";

my %DrugBank_inchikey_drugnames = inchikeys_to_hash($DrugBank_name_inchikey_file, 1);	#1st column is drugname
my %STITCH_inchikey_genesymbols = inchikeys_to_hash($STITCH_inchikey_symbol_file, 2);	#2nd column is Gene Symbol

my %DrugBank_shared_drugs = ();	#drugnames => 1 for shared drugs
my %STDB_targets_shared = ();	#genesymbols => 1 for targets for shared drugs
open my $CMAP, '<', $cMap_name_inchikey_file or die "Could not open file $cMap_name_inchikey_file: $!\n";
open my $CMAP_SHARED_DRUGS, '>', $cMap_shared_drug_file or die "Could not open file $cMap_shared_drug_file: $!\n";
while(<$CMAP>){
	my @words = split(/\t/, $_);
	my $drug = shift @words;
	foreach my $ikey (@words){
		chomp($ikey);
		next if (!$DrugBank_inchikey_drugnames{$ikey} and !$STITCH_inchikey_genesymbols{$ikey});	#skip in inchikey does not appear in either DrugBank or STITCH
		print $CMAP_SHARED_DRUGS $drug, "\t", $ikey,  "\n";	#print output for cmap
		$DrugBank_shared_drugs{$DrugBank_inchikey_drugnames{$ikey}} = 1 if $DrugBank_inchikey_drugnames{$ikey};
		$STDB_targets_shared{$STITCH_inchikey_genesymbols{$ikey}} = 1 if $STITCH_inchikey_genesymbols{$ikey};	#get shared targets for STITCH first
	}
}
close $CMAP_SHARED_DRUGS;
close $CMAP;

open my $DrugBank_nametarget, '<', $DrugBank_name_target_file or die "Could not open file $DrugBank_name_target_file: $!\n";
my %DrugBank_targets_shared = ();
while(<$DrugBank_nametarget>){
	my @words = split(/\t/, $_);
	my $drug = shift @words;
	next unless $DrugBank_shared_drugs{$drug};	#skip lines not shared
	foreach my $genesymbol (@words){
		chomp($genesymbol);
		$STDB_targets_shared{$genesymbol} = 1;	#get shared targets for DrugBank here
	}
}	#now %STDB_targets_shared contains genesymbol => 1 pair for drug targets for shared drugs
close $DrugBank_nametarget;

open my $STDBoutput, '>', $STDB_source_file or die "Could not open file $STDB_source_file: $!\n";
foreach my $gs (keys %STDB_targets_shared){
	chomp($gs);
	print $STDBoutput $gs, "\n";	#print unique gene symbols (hash keys are not duplicated)
}
close $STDBoutput;


#----------------------------------------Subroutines-----------------------------------------------------------------------------------------------------------------------
sub inchikeys_to_hash {	#takes two arguments, filename and column number for hash values (NOT the index)
	my $file = shift @_;
	my $value_col = shift @_;
	$value_col--;

	my %inchikey_values = ();
	open my $FILE, '<', $file or die "Could not open file $file: $!\n";
	while (<$FILE>){
		my @words = split(/\t/, $_);
		my $value = $words[$value_col];	#the drugname
		chomp($value);
		foreach my $word ( @words ){
			chomp($word);
			next unless $word =~ m/^[A-Z]{14}\-[A-Z]{10}\-[A-Z]$/; #skip if not InChIKey format
			$inchikey_values{$word} = $value;
		}
	}
	close $FILE;
	return %inchikey_values;
}

exit;
