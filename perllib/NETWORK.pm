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
my $is_intersect = 1;	#0 if using union; union means drugs appearing in both cMap and (DrugBank 'union' STITCH)

die "Please specify the parent folder for network files\n" unless @ARGV == 1;
my $outdir = shift @ARGV;
chomp($outdir);
$outdir = dirname_add_slash($outdir);
make_dir($outdir);
my $cMap_obj = cMap->new("100up");
my $DrugBank_obj = DrugBank->new();
my $STITCH_obj = STITCH->new();
my $String_obj = String->new();
get_network($cMap_obj, $DrugBank_obj, $STITCH_obj, $String_obj);
#-----------------------------------------------------TEST AREA-------------------------------------------------------------------
#my $cmapobj = cMap->new("50up");
#my $dbobj = DrugBank->new();
#my $STobj = STITCH->new();
#my $presrcref = get_pre_source($cmapobj, $dbobj, $STobj);
#print Dumper($presrcref);
#-----------------------------------------------------TEST AREA-------------------------------------------------------------------
sub get_network
{
	my ($cmap_ref, $drugbank_ref, $stitch_ref, $string_ref) = @_;	#the references to databases
	my @cMap_InChIKeys = $cmap_ref->get_cMap_InChIKey();
	my $type = "union";
	$type = "intersect" if $is_intersect == 1;
	my $path_commondrug = "./static/common_drugs/";
	my $cmap_drug = $path_commondrug . $type . "/cMap_drugnames.txt";
	my $drugbank_drug = $path_commondrug . $type . "./DrugBank_drugnames.txt";
	my $stitch_drug = $path_commondrug . $type . "./STITCH_InChIKeys.txt";

	my $drugbank_exist = one_column_file_switch($drugbank_drug);
	my $stitch_exist = one_column_file_switch($stitch_drug);
	foreach my $ikey (@cMap_InChIKeys){	#at this point it exists in cmap data; no need to check if it is in cmap
		my $db_drug = $drugbank_ref->get_DrugBank_drugname_by_InChIKey($ikey);	#may or may not exist
		next if ($db_drug == 0 && !$stitch_exist->{$ikey});	#skip if not found in either DrugBank or STITCH
	
		#found at least in one database (DrugBank, STITCH)
		my $pre_network = get_pre_source_dest_by_InChIKey($cmap_ref, $drugbank_ref, $stitch_ref, $ikey);
		next unless $pre_network;
		my $pre_src_ref = $pre_network->{"pre_source"};
		my @pre_src = @$pre_src_ref;
		my $pre_dest_ref = $pre_network->{"pre_dest"};
		my @pre_dest = @$pre_dest_ref;
		
		#now need to check if source and dest genes appear in PPI
		my @sources;	#final sources
		my @destinations;	#destinations for the given drug
		my @edges;
		foreach my $source (@pre_src){
			chomp($source);
			my $edge_ref = $string_ref->get_String_edges_by_single_node($source);
			next unless $edge_ref;	#skip if not found in String (0 is returned)
			push (@sources, $source);
			push (@edges, @$edge_ref);
		}
		foreach my $dest (@pre_dest){
			chomp($dest);
			my $edge_ref = $string_ref->get_String_edges_by_single_node($dest);
			next unless $edge_ref;	#skip if not found, 0 is returned
			push (@destinations, $dest);
			push (@edges, @$edge_ref);
		}
		my @pre_nodes;
		push (@pre_nodes, @sources);
		push (@pre_nodes, @destinations);
		my @uniq_pre_nodes = unique(\@pre_nodes);
		my @nodes;
		foreach my $node (@uniq_pre_nodes){
			chomp($node);
			my $node_line = $node . "\t1";	#1 for the fold change
			push (@nodes, $node_line);
		}
		#now ready for output networks
		my $num_node = scalar(@nodes);	#Number of Nodes: $num_node
		my $num_edge = scalar(@edges);	#Number of Edges: $num_edge
		my $current_drug = $cmap_ref->get_cMap_drugname_by_InChIKey($ikey);
		$current_drug = rm_special_char_in_drugname($current_drug);
		chomp($current_drug);
		my $current_dir = $outdir . $current_drug;
		$current_dir = dirname_add_slash($current_dir);
		make_dir($current_dir);	#create directory for individual drug

		my $source_file = $outdir . $current_drug . "/source_" . $current_drug . ".txt";
		my $dest_file = $outdir . $current_drug . "/dest_" . $current_drug . ".txt";
		my $node_file = $outdir . $current_drug . "/node_" . $current_drug . ".txt";
		my $edge_file = $outdir . $current_drug . "/edge_" . $current_drug . ".txt";
		open my $SRC, '>', $source_file or die "Could not open source file $source_file: $!\n";
		foreach my $src (@sources){
			print $SRC $src, "\n";
		}
		close $SRC;
		open my $DST, '>', $dest_file or die "Could not open dest file $dest_file: $!\n";
		foreach my $dst (@destinations){
			print $DST $dst, "\n";
		}
		close $DST;
		open my $NOD, '>', $node_file or die "Could not open node file $node_file: $!\n";
		print $NOD "Number of Nodes: $num_node", "\n";
		foreach my $node (@nodes){
			print $NOD $node, "\n";
		}
		close $NOD;
		open my $EDG, '>', $edge_file or die "Could not open edge file $edge_file: $!\n";
		print $EDG "Number of Edges: $num_edge", "\n";
		foreach my $edge (@edges){
			print $EDG $edge, "\n";
		}
		close $EDG;
	}
	return;
}
sub get_pre_source_dest_by_InChIKey
{
	#usage: my $prenetwork = get_pre_source_dest_by_InChIKey($....);
	#	my $pre_src_ref = $prenetwork->{"pre_source"};
	#	my @pre_src = @$pre_src_ref;	#source genes for given InChIKey (drug)
	#
	#	my $pre_dest_ref = $prenetwork->{"pre_dest"};
	#	my @pre_dest = @$pre_dest_ref;	#destination genes for given InChIKey (drug)
	my ($cmap_ref, $drugbank_ref, $stitch_ref, $ikey) = shift @_;
	
	my $dest_ref = $cmap_ref->get_cMap_targets_by_InChIKey($ikey);	#ref to an array of targets in cmap
	my @pre_dest = @$dest_ref;	#targets in cmap
	my @pre_sources;

	#may produce warnings or errors here when: 1)no targets for the given drug, 2)only one database contains the drug--in case of union
	my $drugbank_targets_ref = $drugbank_ref->get_DrugBank_targets_by_InChIKey($ikey);	#hash ref
	my $stitch_targets_ref = get_STITCH_targets_by_InChIKey($ikey);	#array ref
	foreach my $drug (keys %$drugbank_targets_ref){
		chomp($drug);
		push (@pre_sources, $drug);
	}
	push (@pre_sources, @$stitch_targets_ref);	#join array
	my @temp = @pre_sources;
	@pre_sources = unique(\@temp);	#remove redundancies
	my %hash;
	$hash{"pre_source"} = \@pre_sources;
	$hash{"pre_dest"} = \@pre_dest;
	return \%hash;
}
1;
