#!/usr/bin/perl

package NETWORK;
use cMap;
use DrugBank;
use STITCH;
use String;
use DrugTargetBase;
use Data::Dumper;
use strict;
use warnings;

#-----------------------------------------------------TEST AREA-------------------------------------------------------------------
my $cmapobj = cMap->new("50up");
my $dbobj = DrugBank->new();
my $STobj = STITCH->new();
my $presrcref = get_pre_source($cmapobj, $dbobj, $STobj);
print Dumper($presrcref);
#-----------------------------------------------------TEST AREA-------------------------------------------------------------------

sub get_pre_source
{
	my ($cmap_ref, $drugbank_ref, $stitch_ref) = @_;	#the references to databases
	my @cMap_InChIKeys = $cmap_ref->get_cMap_InChIKey();
	my %pre_sources;
	foreach my $ikey (@cMap_InChIKeys){
		my $drugname = $cmap_ref->get_cMap_drugname_by_InChIKey($ikey);	#each drugname on cMap
		my @genes;
		my $dbtarget_ref = $drugbank_ref->get_DrugBank_targets_by_InChIKey($ikey);
		my %dbtargets = %$dbtarget_ref;
		foreach my $target (keys %dbtargets){
			next if $dbtargets{$target} =~ m/unknown/i;
			next if $dbtargets{$target} =~ m/other/i;
			next if $dbtargets{$target} =~ m/product\s*of/i;
			push @genes, $target;
		}
		my $sttarget_ref = $stitch_ref->get_STITCH_targets_by_InChIKey($ikey);
		my @sttargets = @$sttarget_ref;
		foreach my $target2 (@sttargets){
			next if $stitch_ref->get_STITCH_score($ikey, $target2) < 900;	#at least 900
			push @genes, $target2;
		}
		my @uniq_pre_srcs = unique(@genes);
		$pre_sources{$drugname} = [@uniq_pre_srcs];
	}
	return \%pre_sources;
}
sub get_destination
{

}
sub get_node
{

}
sub get_edge
{
	my $string = String->new();
	my $G1 = "ARF5";
	my $G2 = "VAMP3";
	my @nodes = ('PSD1', 'VAMP3', 'ACAP1', 'QCR1', 'MLF2');
	print $string->get_String_distance($G1, $G2);   #0.673
	my $edgeref = $string->get_String_edges_by_nodes(\@nodes);
	my %edges = %$edgeref;
	foreach my $node1 (keys %edges){
		foreach my $node2 (keys $edges{$node1}){
			print "$node1\t$node2\t", $edges{$node1}{$node2}, "\n";
		}
	}
}

sub printout_network
{

}

1;
