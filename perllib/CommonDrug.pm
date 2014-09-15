#!/usr/bin/perl

package CommonDrug;
use cMap;
use DrugBank;
use STITCH;
use String;
use DrugTargetBase;
use Data::Dumper;
use strict;
use warnings;

my $cmap50up = new cMap("50up");	#the range (50up) does not matter for extraction of common drugs
my $drugbank = new DrugBank();
my $stitch = new STITCH();

get_common_drugs_by_union($cmap50up, $drugbank, $stitch);
get_common_drugs_by_intersection($cmap50up, $drugbank, $stitch);

sub get_common_drugs_by_union
{
        #input: objects (references) to DrugBank, STITCH and cMap
        #output: list of drug names for each database of drugs in both (cMap) and (DrugBank U STITCH)
        my ($cmap_ref, $drugbank_ref, $stitch_ref) = @_;        #the references to databases
        my @cMap_InChIKeys = $cmap_ref->get_cMap_InChIKey();
        my (%cMap_drugnames, %DrugBank_drugnames, %STITCH_drugikeys);
        foreach my $ikey (@cMap_InChIKeys){
                my $cmapdrug = $cmap_ref->get_cMap_drugname_by_InChIKey($ikey);
                my $DBdrug = $drugbank_ref->get_DrugBank_drugname_by_InChIKey($ikey);
                my $STdrugtarget_ref = $stitch_ref->get_STITCH_targets_by_InChIKey($ikey);
                if ( defined($DBdrug) ){
                        $DrugBank_drugnames{$DBdrug} = 1;
                        $cMap_drugnames{$cmapdrug} = 1;
                }
                if ( defined($STdrugtarget_ref) ){
                        $STITCH_drugikeys{$ikey} = 1;
                        $cMap_drugnames{$cmapdrug} = 1;
                }
        }
	my $dir = "./static/common_drugs/union/";
	print_common_drugs(\%cMap_drugnames, \%DrugBank_drugnames, \%STITCH_drugikeys, $dir);
	return;
}
sub get_common_drugs_by_intersection
{
        #input: objects (references) to DrugBank, STITCH and cMap
        #output: list of drug names for each database of drugs in both (cMap) and (DrugBank && STITCH)
        my ($cmap_ref, $drugbank_ref, $stitch_ref) = @_;        #the references to databases
        my @cMap_InChIKeys = $cmap_ref->get_cMap_InChIKey();
        my (%cMap_drugnames, %DrugBank_drugnames, %STITCH_drugnames);
        foreach my $ikey (@cMap_InChIKeys){
                my $cmapdrug = $cmap_ref->get_cMap_drugname_by_InChIKey($ikey);
                my $DBdrug = $drugbank_ref->get_DrugBank_drugname_by_InChIKey($ikey);
                my $STdrug = $stitch_ref->get_STITCH_drugname_by_InChIKey($ikey);
                if ( defined($DBdrug) && defined($STdrug) ){
                        $DrugBank_drugnames{$DBdrug} = 1;
                        $STITCH_drugnames{$STdrug} = 1;
                        $cMap_drugnames{$cmapdrug} = 1;
                }
        }
	my $dir = "./static/common_drugs/intersection/";
	print_common_drugs(\%cMap_drugnames, \%DrugBank_drugnames, \%STITCH_drugnames, $dir);
	return;
}
sub print_common_drugs
{
        my ($cMap_drugnames_ref, $DrugBank_drugnames_ref, $STITCH_drugnames_ref, $dir);
	make_dir($dir);
	my $cmapfile = $dir . "cMap_commondrugs.txt";
	my $dbfile = $dir . "DrugBank_commondrugs.txt";
	my $stfile = $dir . "STITCH_commondrugs.txt";
	my %cMap_drugnames = %$cMap_drugnames_ref;
	my %DrugBank_drugnames = %$DrugBank_drugnames_ref;
	my %STITCH_drugnames = %$STITCH_drugnames_ref;
	open my $CMAP_FH, '>', $cmapfile or die "Could not open cmap commondrug file $cmapfile: $!\n";
	foreach my $cdrug (keys %cMap_drugnames){
		print $CMAP_FH $cdrug, "\n";
	}
	close $CMAP_FH;
	open my $DB_FH, '>', $dbfile or die "Could not open drugbank commondrug file $dbfile: $!\n";
	foreach my $dbdrug (keys %DrugBank_drugnames){
		print $DB_FH $dbdrug, "\n";
	}
	close $DB_FH;
	open my $ST_FH, '>', $stfile or die "Could not open stitch commondrug file $stfile: $!\n";
	foreach my $stdrug (keys %STITCH_drugnames){
		print $ST_FH $stdrug, "\n";
	}
	close $ST_FH;
	return;
}
1;
