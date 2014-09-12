#!/usr/bin/perl

package cMap;
use DrugTargetBase;
use PUGREST;
use Data::Dumper;
use strict;
use warnings;
#----------------------------------------------------------TEST AREA-----------------------------
#my $input = shift @ARGV;
#my $cmap100up = new cMap("100up");
#my $ikey = "PIRYBCJVLMHZOK-UHFFFAOYSA-N";
#print $cmap100up->get_cMap_drugname_by_InChIKey($ikey);
#print Dumper($cmap100up->get_cMap_targets_by_InChIKey($ikey));
#----------------------------------------------------------TEST AREA-----------------------------
sub new
{
	my $class = shift;
	my $range = shift;
	my $self = cMapdata($range);
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
sub cMapdata
{	
	#to create the $self body
	#this returns a reference to the hash containing cMap information
	my $range = shift @_;
	die "Please specify cMap range in new cMap(\"range\")\n"unless $range;
	my $file = "unknown";
	#decide which file to open based on the range input
	$file = "./static/cMap/cMap_drugRL_top50.txt" if $range eq "50up";
	$file = "./static/cMap/cMap_drugRL_top100.txt" if $range eq "100up";
	$file = "./static/cMap/cMap_drugRL_bot50.txt" if $range eq "50down";
	$file = "./static/cMap/cMap_drugRL_bot100.txt" if $range eq "100down";
	die "cMap ranges can be one of 50up, 100up, 50down, and 100down\n" if $file eq "unknown";

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
				$cMap{$ikey}{drugname} = $drug;
				$cMap{$ikey}{targets} = [@targets]; 
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
			push @targets, $target;
		}
	}
	close $CMAP;
	foreach my $ikey (@ikeys) {#for the last drug
		$cMap{$ikey}{drugname} = $drug;
		$cMap{$ikey}{targets} = [@targets]; 
	}
	my $cMapref = \%cMap;
	return $cMapref;
}

1;
