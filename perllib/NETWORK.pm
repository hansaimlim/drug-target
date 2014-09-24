#!/usr/bin/perl

package NETWORK;
use cMap;
use DrugBank;
use STITCH;
use String;
use DrugTargetBase;
use strict;
use warnings;
my $is_intersect = 1;	#0 if using union; union means drugs appearing in both cMap and (DrugBank 'union' STITCH)

die "Usage: $0 <parent directory for network files> <cMap range; 100up 100down..>\n" unless @ARGV == 2;
my $outdir = shift @ARGV;
my $range = shift @ARGV;
chomp($outdir);
chomp($range);
$outdir = dirname_add_slash($outdir);
make_dir($outdir);

#need to be modified---------use json to load hash structure instead of primitive objects
my $cMap_obj = cMap->new($range, "off");	#matching range file must exist and defined in cMap.pm . check if error occurs
my $DrugBank_obj = DrugBank->new("off");	#off means PUGREST off
#need to be modified---------use json to load hash structure instead of primitive objects

my $STITCH_obj = new STITCH();	#STITCH obj does not need PUGREST--json NOT needed
my $String_obj = new String();	#String obj does not need PUGREST--json NOT needed
get_network($cMap_obj, $DrugBank_obj, $STITCH_obj, $String_obj);


sub get_network
{
	my ($cmap_ref, $drugbank_ref, $stitch_ref, $string_ref) = @_;	#the references to databases
	my @cMap_InChIKeys = keys (%$cmap_ref);
	my $type = "union";
	$type = "intersect" if $is_intersect == 1;
	my $path_commondrug = "./static/common_drugs/";
	my $cmap_drug = $path_commondrug . $type . "/cMap_drugnames.txt";
	my $drugbank_drug = $path_commondrug . $type . "/DrugBank_drugnames.txt";
	my $stitch_drug = $path_commondrug . $type . "/STITCH_InChIKeys.txt";

	my $drugbank_exist = one_column_file_switch($drugbank_drug);
	my $stitch_exist = one_column_file_switch($stitch_drug);
	foreach my $ikey (@cMap_InChIKeys){	#at this point it exists in cmap data; no need to check if it is in cmap
		my $db_drug = $drugbank_ref->get_DrugBank_drugname_by_InChIKey($ikey);	#may or may not exist
		next if (!$db_drug && !$stitch_exist->{$ikey});	#skip if not found in either DrugBank or STITCH
	
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
		foreach my $source (@pre_src){
			chomp($source);
			my $edge_ref = $string_ref->get_String_edges_by_single_node($source);
			next unless $edge_ref;	#skip if not found in String (0 is returned)
			push (@sources, $source);
		}
		foreach my $dest (@pre_dest){
			chomp($dest);
			my $edge_ref = $string_ref->get_String_edges_by_single_node($dest);
			next unless $edge_ref;	#skip if not found, 0 is returned
			push (@destinations, $dest);
		}
		my $num_source = scalar(@sources);
		my $num_dest = scalar(@destinations);

		if ($num_source == 0 or $num_dest == 0){	#skip if no source of destination
			print "No sources or destinations found for $db_drug\n";
			next;
		}
	
		my $node_edge_ref = get_node_and_edges(\@sources, \@destinations, $string_ref);
		my $node_ref = $node_edge_ref->{"nodes"};
		my $edge_ref = $node_edge_ref->{"edges"};
		my @edges = @$edge_ref;
		my @nodes;
		foreach my $node (@$node_ref){
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
sub get_node_and_edges
{
	#input: source_ref, dest_ref, and String_obj
	#output: hash ref containing ref to nodes and edges
	my ($source_ref, $dest_ref, $string_obj) = @_;
	my @sources = @$source_ref;
	my @destinations = @$dest_ref;
	my %node_seen;	#node switch
	my %edge_seen;	#edge switch
	my @current_nodes = (@sources, @destinations);
	@current_nodes = unique(\@current_nodes);
	my @current_edges;
	my @new_nodes;
	my @new_edges;
	foreach my $node (@current_nodes){
		chomp($node);
		$node_seen{$node} = 1;	#node switch on
		my $edge_ref = $string_obj->get_String_edges_by_single_node($node);
		my @edge = @$edge_ref;
		if (scalar(@edge) == 0){	#skip if no edge found
			next;
		}
		foreach my $edge (@edge){
			chomp($edge);
			push (@current_edges, $edge);
			$edge_seen{$edge} = 1;	#edge switch on
		}
	}
	foreach my $ed (@current_edges){
		chomp($ed);
		my @words = split(/\t/, $ed);
		my $n1 = shift @words;	#node 1
		my $n2 = shift @words;	#node 2
		#my $ds = shift @words;	#distance ---unnecessary
		if (!$node_seen{$n1}){
			push (@new_nodes, $n1);
			$node_seen{$n1} = 1;
		}
		if (!$node_seen{$n2}){
			push (@new_nodes, $n2);
			$node_seen{$n2} = 1;
		}
	}
	my $num_new_node = scalar(@new_nodes);	#0 if no more new nodes needed
	while ($num_new_node){
		foreach my $node (@new_nodes){
			chomp($node);
			$node_seen{$node} = 1;	#node switch on
			my $edge_ref = $string_obj->get_String_edges_by_single_node($node);
			my @edge = @$edge_ref;
			if (scalar(@edge) == 0){	#skip if no edge found
				next;
			}
			foreach my $edge (@edge){
				chomp($edge);
				push (@new_edges, $edge);
				$edge_seen{$edge} = 1;	#edge switch on
			}
		}
		push (@current_nodes, @new_nodes);
		undef @new_nodes;
		foreach my $ed (@new_edges){
			chomp($ed);
			my @words = split(/\t/, $ed);
			my $n1 = shift @words;	#node 1
			my $n2 = shift @words;	#node 2
			#my $ds = shift @words;	#distance unnecessary
			if (!$node_seen{$n1}){
				push (@new_nodes, $n1);
				$node_seen{$n1} = 1;
			}
			if (!$node_seen{$n2}){
				push (@new_nodes, $n2);
				$node_seen{$n2} = 1;
			}
		}
		push (@current_edges, @new_edges);
		undef @new_edges;
		$num_new_node = scalar(@new_nodes);
	}
	my @final_edges = unique(\@current_edges);
	my @final_nodes = unique(\@current_nodes);
	my %hash;
	$hash{"edges"} = \@final_edges;
	$hash{"nodes"} = \@final_nodes;
	return \%hash;
}
sub get_pre_source_dest_by_InChIKey
{
	#usage: my $prenetwork = get_pre_source_dest_by_InChIKey($....);
	#	my $pre_src_ref = $prenetwork->{"pre_source"};
	#	my @pre_src = @$pre_src_ref;	#source genes for given InChIKey (drug)
	#
	#	my $pre_dest_ref = $prenetwork->{"pre_dest"};
	#	my @pre_dest = @$pre_dest_ref;	#destination genes for given InChIKey (drug)
	my ($cmap_ref, $drugbank_ref, $stitch_ref, $ikey) = @_;

	my $dest_ref = $cmap_ref->get_cMap_targets_by_InChIKey($ikey);	#ref to an array of targets in cmap
	my @pre_dest = @$dest_ref;	#targets in cmap
	my @pre_sources;

	#may produce warnings or errors here when: 1)no targets for the given drug, 2)only one database contains the drug--in case of union
	my $drugbank_targets_ref = $drugbank_ref->get_DrugBank_targets_by_InChIKey($ikey);	#hash ref
	my $stitch_targets_ref = $stitch_ref->get_STITCH_targets_by_InChIKey($ikey);	#array ref
	return 0 if ($drugbank_targets_ref == 0 && $stitch_targets_ref == 0);	#if no targets found for the given inchikey
	if ($drugbank_targets_ref){	#skip if no drugbank targets found
		foreach my $drug (keys %$drugbank_targets_ref){
			chomp($drug);
			push (@pre_sources, $drug);
		}
	}
	push (@pre_sources, @$stitch_targets_ref) if $stitch_targets_ref;	#push stitch targets if found
	return 0 if scalar(@pre_sources) == 0;	#skip if no targets found

	my @temp = @pre_sources;
	@pre_sources = unique(\@temp);	#remove redundancies
	my %hash;
	$hash{"pre_source"} = \@pre_sources;
	$hash{"pre_dest"} = \@pre_dest;
	return \%hash;
}
1;
