#!/usr/bin/perl

package IDMAP;
use DrugTargetBase;
use warnings;
use Data::Dumper;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(manual_get_genename_by_UniProtKB manual_get_InChIKey_by_chemicalID);

#------------------------------------------------TEST AREA------------------------------------
#my $u = "GGACT_HUMAN";
#my $g = get_genename_by_UniProtKB($u);	#A2LD1
#print $g, "\n";
#print get_genename_by_UniProtKB("nothing"), "\n";	# 0
#print get_genename_by_UniProtKB("GGACT"), "\n";	#A2LD1
#------------------------------------------------TEST AREA------------------------------------
sub manual_get_InChIKey_by_chemicalID
{
	#use manually found IDmap file in idmap static folder
	#input  : chemical ID (e.g. Chicago_Sky_Blue_6B)
	#output : InChIKey    (e.g. UPKAWFACSJWKND-MAQYYZNPSA-J)
	my $chemical = shift @_;
	my $hash_ref = InChIKey_by_chemicalID();	#hash reference (chemical ID => InChIKey)
	my $ikey = $hash_ref->{$chemical};
	return unless defined($ikey);
	return if $ikey =~ m/unknown/i;
	return $ikey;
}
sub InChIKey_by_chemicalID
{
	my $file = "./static/idmap/cMap_manuallyfound_maps.tsv";
	my %hash;
	open my $ID, '<', $file or die "Could not open Idmap file $file: $!\n";
	while (my $line = <$ID>){
		my @words = split(/\t/, $line);
		my $chemical = shift @words;
		my $ikey = shift @words;
		chomp($ikey);
		$hash{$chemical} = $ikey;
	}
	close $ID;
	return \%hash;
}
sub manual_get_genename_by_UniProtKB
{
	#input  : UniProtKB (one)
	#output : genename (one) or 0
	my $uniprot = shift @_;
	$uniprot = remove_HUMAN_suffix($uniprot);
	my $ug_ref = UniProtKB_genename();
	my $genename = $ug_ref->{$uniprot};
	if (defined($genename)){
		return $genename;
	} else {
		return 0;
	}	
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
