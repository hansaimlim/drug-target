#!/usr/bin/perl

package IDMAP;
use DrugTargetBase;
use strict;
use warnings;
use Data::Dumper;
#------------------------------------------------TEST AREA------------------------------------
#my $b = genename_UniProtKB();
my $a = UniProtKB_genename();
print $a->{"GGACT"};	#A2LD1
print $a->{"SYAC"};	#AARS
#------------------------------------------------TEST AREA------------------------------------

sub genename_UniProtKB
{
	#returns a reference to hash containing (genename => UniProtKB) pairs
	my $file = "./static/idmap/genename_UniProtsymbol.tsv";
	my %genename_UniProtKB;
	open my $GU, '<', $file or die "Could not open Idmap file $file: $!\n";
	while(my $line = <$GU>){
		next if $. == 1;
		my @words = split(/\t/, $line);
		my $genename = shift @words;
		my $UniProtKB = shift @words;
		chomp($UniProtKB);
		if ($UniProtKB =~ m/_HUMAN/i){
			$UniProtKB =~ s/_HUMAN//ig;	#remove the suffix if exist
		}
		next if $genename eq $UniProtKB;	#skip meaningless conversion (i.e. ABL1 => ABL1)
		$genename_UniProtKB{$genename} = $UniProtKB;
	}
	close $GU;
	return \%genename_UniProtKB;
}
sub UniProtKB_genename
{
	#returns a reference to hash containing (UniProtKB => genename) pairs
	my $gu_ref = genename_UniProtKB();
	my $ug_ref = reverse_simple_hash($gu_ref);
	return $ug_ref;
}
1;
