#!/usr/bin/perl

package String;
use DrugTargetBase;
use PUGREST;
use IDMAP;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;
use strict;
use warnings;
#----------------------------------------------------------TEST AREA-----------------------------
#my $input = shift @ARGV;
#my $string = new String();
#my $G1 = "ARF5";
#my $G2 = "VAMP3";
#my @nodes = ('PSD1', 'VAMP3', 'ACAP1', 'QCR1', 'MLF2');
#print $string->get_String_distance($G1, $G2);	#0.673
#my $edgeref = $string->get_String_edges_by_nodes(\@nodes);
#my %edges = %$edgeref;
#my $num_edges = 0;
#foreach my $node1 (keys %edges){
#	foreach my $node2 (keys $edges{$node1}){
#		print "$node1\t$node2\t", $edges{$node1}{$node2}, "\n";
#		$num_edges++;
#	}
#}
#print $num_edges,"\n";
#----------------------------------------------------------TEST AREA-----------------------------
my $is_demo_on = 1;	#read demo data (shorter list) if 1
sub new
{
        my $class = shift @_;
        my $self = StringData();
        bless $self, $class;
        return $self;
}
sub get_String_edges_by_single_node
{
	#input : a single node reference (genesymbol)
	#output: reference to array of lines (edges)
	my( $self, $node ) = @_;
	my $g1 = $node;
	chomp($g1);
	my @edges;	#will contain sorted edges as lines (one element contains g1 g2 and distance, but no new line character)
	if ( !( defined($self->{$g1})) ){
		return 0;	#not found in the PPI, return 0
	}
	if (defined($self->{$g1})){
		my $ref = $self->{$g1};
		my %n = %$ref;
		my @temp;
		foreach my $g2 (keys %n){
			chomp($g2);
			push (@temp, $g1);
			push (@temp, $g2);
			my @sorted_nodes = sort(@temp);	#sort to remove redundant edges
			my $first = shift @sorted_nodes;
			chomp($first);
			my $second = shift @sorted_nodes;
			chomp($second);
			my $distance = $n{$g2};
			my $edge = $first . "\t" . $second . "\t" . $distance;
			push (@edges, $edge);
			undef @sorted_nodes;
			undef @temp;
		}
	} 
	my @uniq_edges = unique(\@edges);
        return \@uniq_edges;
}
sub get_String_edges_by_nodes
{
	#preserving redundant edges is useful here
	#final set of edges will not include redundant edges
        #input : reference to array of nodes (GeneSymbol)
        #output: reference to hash of edges
        my( $self, $noderef ) = @_;
        my @nodes = @$noderef;	#nodes
	my %edges;
	foreach my $node (@nodes){
		my $g1 = $node;
		if (defined($self->{$node})){
			my $ref = $self->{$node};
			my %n = %$ref;
			foreach my $g2 (keys %n){
				$edges{$g1}{$g2} = $n{$g2};
			}
		} 
	}
        return \%edges;
}
sub get_String_distance
{
        #input : Two GeneSymbols
        #output: Distance or 0
        my( $self, $gene1, $gene2 ) = @_;
        if (defined($self->{$gene1}->{$gene2})){
                my $distance = $self->{$gene1}->{$gene2};
                return $distance;
        } else {
                return 0;
        }
}
sub Randomize_String
{
	#randomize PPI network
	#follows the algorithm introduced in Random Network Plugin for Cytoscape - reference below:
	#https://sites.google.com/site/randomnetworkplugin/Home/randomization-of-existing-networks
	my( $self, $n ) = @_;	#will shuffle n times
	my %ppi = %$self;	#do not use reference below; to preserve original ppi
        my $file = "./static/String/9606.protein.links.v9.1-GN-dist.txt";
        $file = "./static/demo/nonrandomppi.tsv" if $is_demo_on;
	my @edges_original;
	my @distances;	#collection of distances; duplicate values are allowed
	my %edge_count;
	open my $String, '<', $file or die "Could not open PPI file, $file: $!\n";
	while (my $line = <$String>){
		next if $. == 1;	#skip first line
		chomp($line);
		push (@edges_original, $line);
		$edge_count{$line}++;

		my @words = split(/\t/, $line);
		my $dist = $words[2];	#the third elements are distances (with a possible new line)
		chomp($dist);
		push (@distances, $dist);
	}
	close $String;
	SHUFFLE: for (my $i=1; $i <= $n; $i++){
		my @edges = fisher_yates_shuffle(\@edges_original);
		my @distances_shuffled = fisher_yates_shuffle(\@distances);
		my ($edge1, $edge2);
		my ($new_edge1, $new_edge2, $dist1, $dist2);
		my ($u, $v, $s, $t);	#nodes
		EDGE: while(@edges){
			my $first_edge = shift @edges;
			my @words1 = split(/\t/, $first_edge);
			my $u_temp = shift @words1;
			my $v_temp = shift @words1;
			chomp($u_temp);
			chomp($v_temp);
			if ($u_temp ne $v_temp){
				$edge1 = $first_edge;
				$u = $u_temp;
				$v = $v_temp;
			} else {
				next EDGE;
			}
			EDGE2: foreach my $e (@edges){
				my $second_edge = shift @edges;
				my @words2 = split(/\t/, $second_edge);
				my $s_temp = shift @words2;
				my $t_temp = shift @words2;
				chomp($s_temp);
				chomp($t_temp);
				if ($s_temp ne $t_temp){
					if ( ($v ne $s_temp) && ($u ne $t_temp) && ($u ne $s_temp) && ($v ne $t_temp) ){
						next if (($ppi{$u}{$t_temp} or $ppi{$t_temp}{$u} or $ppi{$v}{$s_temp} or $ppi{$s_temp}{$v}));	#should not exist in current ppi
						$edge2 = $second_edge;
						$s = $s_temp;
						$t = $t_temp;
						last EDGE2;
					}
				}
			}
		}
		unless ($u && $v && $s && $t){	#failed to get shuffled edge; perform $i th iteration again
			$i--;
			next SHUFFLE;
		}
		$dist1 = shift @distances_shuffled;
		$dist2 = shift @distances_shuffled;
		chomp($dist1);
		chomp($dist2);
		$ppi{$u}{$t} = $dist1;	
		$ppi{$s}{$v} = $dist2;
		$new_edge1 = "$u\t$t\t$dist1";
		$new_edge2 = "$s\t$v\t$dist2";
		push (@edges, $new_edge1);
		push (@edges, $new_edge2);
		delete($ppi{$u}{$v});	#remove original edge
		delete($ppi{$s}{$t});
	}
	return \%ppi;	#this ppi is separate from the original string data; it is shuffled, but may contain some identical edges.
}
sub StringData
{
        #create String object from pre-converted (Gene Gene Distance)
	#Distances are calculated based on the confidence score: D = (1000-S)/1000
	#At this point, the redundant networks are preserved intentionally
	#The redundancies will be removed when getting edges from node
        my $file = "./static/String/9606.protein.links.v9.1-GN-dist.txt";
        $file = "./static/demo/nonrandomppi.tsv" if $is_demo_on;
        my %Data;
        open my $String, '<', $file or die "Could not open PPI file, $file: $!\n";
        while (my $line = <$String>){
		next if $. == 1;	#skip first line
                my @words = split(/\t/, $line);
                my $gene1 = shift @words;
                my $gene2 = shift @words;
                my $distance = shift @words;
                chomp($distance);
                $Data{$gene1}{$gene2} = $distance;
		$Data{$gene2}{$gene1} = $distance;	#make duplicated data for convenience; Assuming that the network is UNdirected
        }
        close $String;
        return \%Data;
}
1;
