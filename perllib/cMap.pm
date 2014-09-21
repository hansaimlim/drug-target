#!/usr/bin/perl

package cMap;
use DrugTargetBase;
use PUGREST;
use Data::Dumper;
use IDMAP;
use strict;
use warnings;
#----------------------------------------------------------TEST AREA-----------------------------
#my $input = shift @ARGV;
#my $cmap100up = new cMap("100up");
#my $ikey = "PIRYBCJVLMHZOK-UHFFFAOYSA-N";
#print $cmap100up->get_cMap_drugname_by_InChIKey($ikey);
#print Dumper($cmap100up->get_cMap_targets_by_InChIKey($ikey));
#----------------------------------------------------------TEST AREA-----------------------------
my $is_demo_on = 0;	#use demo file (shorter) if 1
my $is_PUGREST_needed = 1;	#0 to turn of PUGREST and speed up reading
sub new
{
	my $class = shift;
	my $range = shift;
	my $self;
	if ($is_PUGREST_needed == 1){
		$self = cMapdata($range);	#performs PUGREST InChIKey converting
	} else {
		$self = cMapSimple($range);	#does not perform PUGREST; use json to load hash
	}
	bless $self, $class;
	return $self;
}
sub get_cMap_InChIKey
{
	#get all the available InChIKeys in cMap DrugRL data
	#returns a reference array
	my( $self ) = @_;
	my @ikeys = keys( %$self );
	return \@ikeys;
}
sub get_cMap_drugname_by_InChIKey
{
	#input InChIKey then output the drugname
	my( $self, $ikey ) = @_;
	my $drugname = $self->{$ikey}->{drugname};
	return $drugname;
}
sub get_cMap_targets_by_InChIKey
{
	#input InChIKey then output a reference to targets in array
	my( $self, $ikey ) = @_;
	my $targetref = $self->{$ikey}->{targets};
	return $targetref;
}
sub cMapSimple
{	
	my $range = shift @_;
	my $file = "unknown";
	#decide which file to open based on the range input
	$file = "./static/json/cMap/100down.json" if $range eq "100down";
	$file = "./static/json/cMap/100up.json" if $range eq "100up";
	$file = "./static/json/cMap/50down.json" if $range eq "50down";
	$file = "./static/json/cMap/50up.json" if $range eq "50up";
	die "cMap range unspecified or file does not exist" if ($file =~ m/unknown/i);
	my $cMap_ref = load_hash($file);
	return $cMap_ref;
}
sub cMapdata
{	
	#to create the $self body
	#this returns a reference to the hash containing cMap information
	my $range = shift @_;
	my $file = "unknown";
	#decide which file to open based on the range input
	$file = "./static/cMap/cMap_drugRL_top50.txt" if $range eq "50up";
	$file = "./static/cMap/cMap_drugRL_bot50.txt" if $range eq "50down";
	$file = "./static/cMap/cMap_drugRL_rand50.txt" if $range eq "50rand";
	$file = "./static/cMap/cMap_drugRL_top100.txt" if $range eq "100up";
	$file = "./static/cMap/cMap_drugRL_bot100.txt" if $range eq "100down";
	$file = "./static/cMap/cMap_drugRL_rand100.txt" if $range eq "100rand";
	$file = "./static/cMap/cMap_drugRL_top1000.txt" if $range eq "1000up";
	$file = "./static/cMap/cMap_drugRL_bot1000.txt" if $range eq "1000down";
	$file = "./static/cMap/cMap_drugRL_rand1000.txt" if $range eq "1000rand";
	print "Using cMap demo file\n" if ($file eq "unknown" or $is_demo_on == 1);
	$file = "./static/cMap/cMap_drugRL_top50_demo.txt" if $is_demo_on;

	my ($is_drugname, $drug) = (0, 0);
	my @ikeys;
	my @targets;
	my %cMap;
	open my $CMAP, '<', $file or die "Could not open cmap file $file: $!\n";
	while(my $line = <$CMAP>){
		if ($. == 1){	#1st line, 'x'
			$is_drugname = 1;
			next;
		}
		if ($line =~ m/^x\s*$/i){
			$is_drugname = 1;
			foreach my $ikey (@ikeys) {
				my @uniq_targets = unique(\@targets);	#remove redundant targets
				$cMap{$ikey}{drugname} = $drug;
				$cMap{$ikey}{targets} = [@uniq_targets]; 
			}
			undef @ikeys;	#should take a new set of InChIKeys
			undef $drug;	#should take a new drug name
			undef @targets;
			next;
		}
		my @words = split(/\t/, $line);
		if ($is_drugname){
			$drug = $words[1];
			chomp($drug);
			@ikeys = get_InChIKey_by_name($drug);	#get InChIKeys
			$is_drugname = 0;
		} else {
			my $target = shift @words;
			chomp($target);
			my $genename = manual_get_genename_by_UniProtKB($target);
			if ($genename){		#if the target format is UniProtKB and listed in IDMAP file
				push @targets, $genename;
			} else {
				push @targets, $target;
			}
		}
	}
	close $CMAP;
	foreach my $ikey (@ikeys) {#for the last drug
		my @uniq_targets = unique(\@targets);	#remove redundant targets
		$cMap{$ikey}{drugname} = $drug;
		$cMap{$ikey}{targets} = [@uniq_targets]; 
	}
	my $cMapref = \%cMap;
	return $cMapref;
}

1;
