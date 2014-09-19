#!/usr/bin/perl

use DrugTargetBase;
use Data::Dumper;
use strict;
use warnings;
#this script output contains redundant drug names (no redundancies for InChIKeys)
#such redundancies are due to the one to many relationships between drug and InChIKeys,
#especially ones categorized as substances under PubChem
#redundancies may be easily removed by linux commands, sort, uniq.
my $CM_ref = read_InChIKey_drugname("./static/idmap/cMap_InChIKey_drugname.tsv");
my $DB_ref = read_InChIKey_drugname("./static/idmap/DrugBank_InChIKey_drugname.txt");
my $ST_ref = read_InChIKey_drugname("./static/STITCH/9606.protein_chemical.links.v4.0InChIKey_GS_min900.tsv");	#format different
my @union_ikey = union_DrugBank_STITCH($CM_ref, $DB_ref, $ST_ref);
my @intersect_ikey = intersect_DrugBank_STITCH($CM_ref, $DB_ref, $ST_ref);

my $union_drugbank_file = "./static/common_drugs/union/DrugBank_drugnames.txt";	#union means (DrugBank Union STITCH)
my $union_cmap_file = "./static/common_drugs/union/cMap_drugnames.txt";
my $union_stitch_file = "./static/common_drugs/union/STITCH_InChIKeys.txt";
print_commondrug($CM_ref, $DB_ref, $union_cmap_file, $union_drugbank_file, $union_stitch_file, \@union_ikey);

my $inter_drugbank_file = "./static/common_drugs/intersect/DrugBank_drugnames.txt";
my $inter_cmap_file = "./static/common_drugs/intersect/cMap_drugnames.txt";
my $inter_stitch_file = "./static/common_drugs/intersect/STITCH_InChIKeys.txt";
print_commondrug($CM_ref, $DB_ref, $inter_cmap_file, $inter_drugbank_file, $inter_stitch_file, \@intersect_ikey);
#foreach my $ikey (@$common_keys_ref){
#	print $cMap_ref->{$ikey}," is for cMap\t";
#	print $DB_ref->{$ikey}," is for DrugBank\n";
#}
sub intersect_DrugBank_STITCH
{
	#input references to InChIKey -> drugname hashes
	#output InChIKeys commonly found in cMap and (DrugBank intersect STITCH)
	my ($cMap_ref, $DrugBank_ref, $STITCH_ref) = @_;	
	my $CD_common_ref = intersection_keys($cMap_ref, $DrugBank_ref);	#CD stands for cMap and DrugBank
	my @CD_arr = @$CD_common_ref;
	my $CS_common_ref = intersection_keys($cMap_ref, $STITCH_ref);	#CS stands for cMap and STITCH
	my @CS_arr = @$CS_common_ref;
	my %seen;
	$seen{$_}++ for @CD_arr;
	my @intersect;
	foreach my $ikey (@CS_arr){
		push (@intersect, $ikey) if $seen{$ikey};
	}
	my @combined_uniq = unique(\@intersect);
	return @combined_uniq;
}
sub union_DrugBank_STITCH
{
	#input references to InChIKey -> drugname hashes
	#output InChIKeys commonly found in cMap and (DrugBank Union STITCH)
	my ($cMap_ref, $DrugBank_ref, $STITCH_ref) = @_;	
	my $CD_common_ref = intersection_keys($cMap_ref, $DrugBank_ref);	#CD stands for cMap and DrugBank
	my @CD_arr = @$CD_common_ref;
	my @CD_uniq = unique(\@CD_arr);
	my $CS_common_ref = intersection_keys($cMap_ref, $STITCH_ref);	#CS stands for cMap and STITCH
	my @CS_arr = @$CS_common_ref;
	my @CS_uniq = unique(\@CS_arr);
	push (@CD_uniq, @CS_uniq);
	my @combined_uniq = unique(\@CD_uniq);
	return @combined_uniq;
}
sub print_commondrug
{
	my ($CM_ref, $DB_ref, $cmap_file, $drugbank_file, $stitch_file, $combined_uniq_ikey_ref) = @_;
	open my $DB, '>', $drugbank_file or die "Could not open file $drugbank_file: $!\n";
	open my $CM, '>', $cmap_file or die "Could not open file $cmap_file: $!\n";
	open my $ST, '>', $stitch_file or die "Could not open file $stitch_file: $!\n";
	foreach my $ikey (@$combined_uniq_ikey_ref){
		print $DB $DB_ref->{$ikey}, "\n";
		print $CM $CM_ref->{$ikey}, "\n";
		print $ST $ikey, "\n";	#STITCH only needs the InChIKey, not a drugname
	}
	close $ST;
	close $CM;
	close $DB;
	return;
}
sub read_InChIKey_drugname
{
	#input file formatted InChIKey (tab) drugname
	#output hash reference
	my $file = shift @_;
	my %hash;
	open my $FH, '<', $file or die "Could not open cmap-inchikey file $file: $!\n";
	while (my $line = <$FH>){
		my @words = split(/\t/, $line);
		my $ikey = shift @words;
		my $drug = shift @words;
		chomp($drug);
		$hash{$ikey} = $drug;
	}
	close $FH;
	return \%hash;
}
sub intersection_keys
{
	#input two hashes (same formatted)
	#output array reference containing commonly appearing keys (InChIKeys)
	my( $hash1_ref, $hash2_ref) = @_;
	my @common_keys;
	foreach my $key1 (%$hash1_ref){
		push (@common_keys, $key1) if defined($hash2_ref->{$key1});
	}
	return \@common_keys;
}
1;
