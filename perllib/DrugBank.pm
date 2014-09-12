#!/usr/bin/perl

package DrugBank;
use DrugTargetBase;
use PUGREST;
use Data::Dumper;
use strict;
use warnings;
#----------------------------------------------------------TEST AREA-----------------------------
#my $input = shift @ARGV;
#my $drugbank = new DrugBank();
#my $ikey = "WGWPRVFKDLAUQJ-UHFFFAOYSA-N";
#my $ikey2 = "WGWPRVFKDLAUQJ-MITYVQBRSA-N";
#print $drugbank->get_DrugBank_drugname_by_InChIKey($ikey);
#print Dumper($drugbank->get_DrugBank_targets_by_InChIKey($ikey));
#print $drugbank->get_DrugBank_drugname_by_InChIKey($ikey2);
#print Dumper($drugbank->get_DrugBank_targets_by_InChIKey($ikey2));
#----------------------------------------------------------TEST AREA-----------------------------
sub new
{
        my $class = shift;
        my $range = shift;
        my $self = DrugBankData();
        bless $self, $class;
        return $self;
}
sub get_DrugBank_drugname_by_InChIKey
{
        #input InChIKey then output the drugname
        my( $self, $ikey ) = @_;
        my $drugname = $self->{$ikey}->{drugname};
        return $drugname;
}
sub get_DrugBank_targets_by_InChIKey
{
        #input InChIKey then output a reference to targets in a hash
        my( $self, $ikey ) = @_;
        my $targetref = $self->{$ikey}->{targets};
        return $targetref;
}
sub DrugBankData
{
	my $file = "./static/DrugBank/DrugBank_name_target_action.tsv";
	my %DrugBankData;
	open my $DrugBank, '<', $file or die "Could not open DrugBank file, $file: $!\n";
	while (my $line = <$DrugBank>){
		my @words = split(/\t/, $line);
		my $drugname = shift @words;
		my %target_action;
		while(@words){
			my $target = shift @words;
			chomp($target);
			my $action = shift @words;
			chomp($action);
			$target_action{$target} = $action;
		}
		my @ikeys = get_InChIKey_by_name($drugname);
		foreach my $ikey (@ikeys){
			chomp($ikey);
			$DrugBankData{$ikey}{drugname} = $drugname;
			foreach my $target (keys %target_action){
				$DrugBankData{$ikey}{targets}{$target} = $target_action{$target};
			}
		}
	
	}
	close $DrugBank;
	my $ref = \%DrugBankData;
	return $ref;
}

1;
