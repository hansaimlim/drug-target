#!/usr/bin/perl
use strict;
use warnings;
use String;
use DrugTargetBase;

my $string = new String();
my %hash = %$string;
my @genes;
foreach my $gene1 (keys %hash){
	chomp($gene1);
	push (@genes, $gene1);
	foreach my $gene2 (keys $hash{$gene1}){
		chomp($gene2);
		push (@genes, $gene2);
	}
}
my @nodes = unique(\@genes);
my $num_nodes = scalar(@nodes);
my $nodefile = "./node.txt";

open my $NODE, '>', $nodefile or die "Could not open node file $nodefile: $!\n";
print $NODE "Number of Nodes: $num_nodes\n";
foreach my $node (@nodes){
	print $NODE "$node\t1\n";	#1 for default fold change
}
close $NODE;

exit;
