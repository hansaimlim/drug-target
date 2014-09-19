#!/usr/bin/perl

package String;
use DrugTargetBase;
use PUGREST;
use IDMAP;
use Data::Dumper;
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
my $is_demo_on = 0;	#read demo data (shorter list) if 1
sub new
{
        my $class = shift;
        my $self = StringData();
        bless $self, $class;
        return $self;
}
sub reverse_edges
{
	#take String PPI hash reference and reverse the order
	#output: reference to hash of reverse edges
	my( $self ) = @_;
	my %ppi = %$self;
	my %rev;
	foreach my $g1 (keys %ppi){
		foreach my $g2 (keys $ppi{$g1}){
			$rev{$g2}{$g1} = $ppi{$g1}{$g2};
		}
	}
	return \%rev;
}
sub get_String_edges_by_single_node
{
	#input : a single node reference (genesymbol)
	#output: reference to array of lines (edges)
	my( $self, $node ) = @_;
	my $g1 = $node;
	chomp($g1);
	my @edges;	#will contain sorted edges as lines (one element contains g1 g2 and distance, but no new line character)
	my $rev = reverse_edges($self);
	if ( !( defined($self->{$g1}) or defined($rev->{$g1})) ){
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
			my @sorted = sort(@temp);	#sort to remove redundant edges
			my $first = shift @sorted;
			chomp($first);
			my $second = shift @sorted;
			chomp($second);
			my $distance = $n{$g2};
			my $edge = $first . "\t" . $second . "\t" . $distance;
			push (@edges, $edge);
		}
	} 
	if (defined($rev->{$g1})) {
		my $ref = $rev->{$g1};
		my %n = %$ref;
		my @temp;
		foreach my $g2 (keys %n){
			chomp($g2);
			push (@temp, $g1);
			push (@temp, $g2);
			my @sorted = sort(@temp);	#sort to remove redundant edges
			my $first = shift @sorted;
			chomp($first);
			my $second = shift @sorted;
			chomp($second);
			my $distance = $n{$g2};
			my $edge = $first . "\t" . $second . "\t" . $distance;
			push (@edges, $edge);
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
	my $rev = reverse_edges($self);	#reverse PPI, a hash reference
	foreach my $node (@nodes){
		my $g1 = $node;
		if (defined($self->{$node})){
			my $ref = $self->{$node};
			my %n = %$ref;
			foreach my $g2 (keys %n){
				$edges{$g1}{$g2} = $n{$g2};
			}
		} 
		if (defined($rev->{$node})) {
			my $ref = $rev->{$node};
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
sub StringData
{
        #create String object from pre-converted (Gene Gene Distance)
	#Distances are calculated based on the confidence score: D = (1000-S)/1000
	#At this point, the redundant networks are preserved intentionally
	#The redundancies will be removed when getting edges from node
        my $file = "./static/String/9606.protein.links.v9.1-GS-dist.txt";
        $file = "./static/String/9606.protein.links.v9.1-GS-dist.txt" if $is_demo_on;
        my %Data;
        open my $String, '<', $file or die "Could not open DrugBank file, $file: $!\n";
        while (my $line = <$String>){
		next if $. == 1;	#skip first line
                my @words = split(/\t/, $line);
                my $gene1 = shift @words;
                my $gene2 = shift @words;
		$gene1 = get_genename_by_UniProtKB($gene1) if get_genename_by_UniProtKB($gene1);	#target id format conversion to genename
		$gene2 = get_genename_by_UniProtKB($gene2) if get_genename_by_UniProtKB($gene2);	#target id format conversion to genename
                my $distance = shift @words;
                chomp($distance);
                $Data{$gene1}{$gene2} = $distance;
        }
        close $String;
        return \%Data;
}
1;
