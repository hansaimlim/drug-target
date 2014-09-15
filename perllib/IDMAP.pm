#!/usr/bin/perl

package IDMAP;
use DrugTargetBase;
use warnings;
use Data::Dumper;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(UniProtKB_genename genename_UniProtKB);

#------------------------------------------------TEST AREA------------------------------------
#my $u = "GGACT_HUMAN";
#my $g = get_genename_by_UniProtKB($u);	#A2LD1
#print $g, "\n";
#------------------------------------------------TEST AREA------------------------------------
sub get_genename_by_UniProtKB
{
	#input  : UniProtKB (one)
	#output : genename (one)
	my $uniprot = shift @_;
	$uniprot = remove_HUMAN_suffix($uniprot);
	my $ug_ref = UniProtKB_genename();
	my $genename = $ug_ref->{$uniprot};
	return $genename;	
}
sub UniProtKB_genename
{	#han@naz:~/drug-target/perllib/static/idmap$ cut -f1 genename_UniProtsymbol.tsv | sort | uniq | wc -l
	#2137 (num of genenames)
	#han@naz:~/drug-target/perllib/static/idmap$ cut -f2 genename_UniProtsymbol.tsv | sort | uniq | wc -l
	#2152 (num of UniProtKBs)
	#returns a reference to hash containing (genename => UniProtKB) pairs
	my $file = "./static/idmap/genename_UniProtsymbol.tsv";
	my %UniProtKB_genename;
	open my $GU, '<', $file or die "Could not open Idmap file $file: $!\n";
	while(my $line = <$GU>){
		next if $. == 1;
		my @words = split(/\t/, $line);
		my $genename = shift @words;
		my $UniProtKB = shift @words;
		chomp($UniProtKB);
		$UniProtKB = remove_HUMAN_suffix($UniProtKB);
		next if $genename eq $UniProtKB;	#skip meaningless conversion (i.e. ABL1 => ABL1)
		$UniProtKB_genename{$UniProtKB} = $genename;
	}
	close $GU;
	return \%UniProtKB_genename;
}
sub genename_UniProtKB
{
	#returns a reference to hash containing (genename => UniProtKB) pairs
	#this action results in loss of some information due to repeated genenames
	my $ug_ref = genename_UniProtKB();
	my $gu_ref = reverse_simple_hash($ug_ref);
	return $gu_ref;
}
sub remove_HUMAN_suffix
{
	my $input = shift @_;
	if ($input =~ m/_HUMAN/i){
		$input =~ s/_HUMAN//ig;
	}
	return $input;
}
1;
