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
#my $cmapobj = cMap->new("50up");
#my $dbobj = DrugBank->new();
#my $STobj = STITCH->new();
#my $presrcref = get_pre_source($cmapobj, $dbobj, $STobj);
#print Dumper($presrcref);
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
		my @uniq_pre_srcs = unique(\@genes);
		$pre_sources{$drugname} = [@uniq_pre_srcs];
	}
	return \%pre_sources;
}
sub get_pre_destination
{
	my ($cmap_ref) = @_;
	my @cMap_InChIKeys = $cmap_ref->get_cMap_InChIKey();
	my %pre_destinations;
	foreach my $ikey (@cMap_InChIKeys){
		my $drugname = $cmap_ref->get_cMap_drugname_by_InChIKey($ikey);
		my $targets_ref = $cmap_ref->get_cMap_targets_by_InChIKey($ikey);
		my @genes = @$targets_ref;
		$pre_destinations{$drugname} = [@genes];
	}
	return  \%pre_destinations;
}
sub make_node
{
	
}
sub get_network
{
	my ($string_ref, $pre_sources_ref, $pre_destinations_ref) = @_;
	my %pre_sources = %$pre_sources_ref;
	my %pre_destinations = %$pre_destinations_ref;
	my (%sources, %destinations, %nodes, %edges);	#containers for final networks

	foreach my $drug (keys %pre_destinations){
		next unless defined($pre_sources{$drug});	#skip if not found in source (STITCH + DrugBank)
		my @pre_dests = $pre_destinations{$drug};
		my @destinations;	#for final destinations
		#get destinations if corresponding edge exists
		foreach my $pre_dest (@pre_dests){
			my $edge_ref = $string_ref->get_String_edges_by_single_node($pre_dest);
			if ($edge_ref){
				push @destinations, $pre_dest;
			}
		}
		#get sources if corresponding edge exists
		my @pre_srcs = $pre_sources{$drug};
		my @sources;
		foreach my $pre_src (@pre_srcs){
			my $edge_ref = $string_ref->get_String_edges_by_single_node($pre_src);
			if ($edge_ref){
				push @sources, $pre_src;
			}
		}
		#get nodes from sources and destinations

		#get edges from nodes

		#insert into the bigger container for each type of networks by drugnames
	}
	
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
