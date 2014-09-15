#!/usr/bin/perl

package STITCH;
use DrugTargetBase;
use PUGREST;
use Data::Dumper;
use strict;
use warnings;
#----------------------------------------------------------TEST AREA-----------------------------
#my $input = shift @ARGV;
#my $stitch = new STITCH();
#my $ikey = "GGAPOKLXWIOMRA-UHFFFAOYSA-N";
#my $target = "OPRK";
#my $ikey2 = "XBSWVSBTGPVAJQ-UHFFFAOYSA-N";
#my $target2 = "CAH12";
#print Dumper $stitch->get_STITCH_targets_by_InChIKey($ikey);	#should print OPRK
#print $stitch->get_STITCH_score($ikey, $target);	#955
#print Dumper $stitch->get_STITCH_targets_by_InChIKey($ikey2);	#should print OPRD
#print $stitch->get_STITCH_score($ikey2, $target2);	#900
#----------------------------------------------------------TEST AREA-----------------------------
sub new
{
        my $class = shift;
        my $self = STITCHData();
        bless $self, $class;
        return $self;
}
sub get_STITCH_targets_by_InChIKey
{
	#input : InChIKey
	#output: reference to an array of targets
        my( $self, $ikey ) = @_;
        my $targetsref = $self->{$ikey};
	my @targets;
	foreach my $target (keys %$targetsref){
		push @targets, $target;
	}
        return \@targets;
}
sub get_STITCH_score
{
	#input : InChIKey and GeneSymbol (a target)
	#output: Prediction score or 0
        my( $self, $ikey, $target ) = @_;
	if (defined($self->{$ikey}->{$target})){
		my $score = $self->{$ikey}->{$target};
		return $score;
	} else {
		return 0;
	}
}
sub STITCHData
{
	#create STITCH object from pre-converted (CIDs to InChIKeys; ENSPs to UniProtKB) with minimum prediction score 900
	#the targets are in genename format; NOT in UniProtKB
        my $file = "./static/STITCH/9606.protein_chemical.links.v4.0InChIKey_GS_min900.tsv";
        my %Data;
        open my $STITCH, '<', $file or die "Could not open DrugBank file, $file: $!\n";
        while (my $line = <$STITCH>){
                my @words = split(/\t/, $line);
                my $ikey = shift @words;
		my $target = shift @words;
		$target = get_genename_by_UniProtKB($target);	#targets are in genenames
		my $score = shift @words;
		chomp($score);
		$Data{$ikey}{$target} = $score;
        }
        close $STITCH;
        return \%Data;
}

1;
