#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
#------------------------------------------------------------------------------------------------------
# Purpose	: To run shortestpath algorithm on multiple drugs
# Input		: directory where sub directories for network files are stored (node, edge, source,dest)
#		  shortestpath output file suffix (ex. 100upmaxdist0.5)
# Output	: shortestpath output for each drug
# Author	: Hansaim Lim
# Date		: 16 Aug, 2014
#------------------------------------------------------------------------------------------------------
die "Usage: $0 <network files directory> <suffix for shortestpath output files> <output directory>\n" unless @ARGV == 3;
my $upperdir = shift @ARGV;
my $filesuffix = shift @ARGV;
my $outdir = shift @ARGV;
SLASH: {
	$upperdir .= "/" unless $upperdir =~ m/\/$/;
	$outdir .= "/" unless $outdir =~ m/\/$/;
}

opendir my $dir1, $upperdir or die "Cannot open directory $upperdir: $!\n";
my @lowerdirs = readdir $dir1;
closedir $dir1;

foreach my $lowerdir ( @lowerdirs ){
	my $currentdir = $upperdir . $lowerdir."/";
	my ($node, $edge, $source, $dest, $outfile_shortestpath);
	my $spathcmd = 'shortestpath';
#for one drug------------------------------------
	opendir my $dir_drug, $currentdir or die "Cannot open directory $currentdir: $!\n";
	my @files = readdir $dir_drug;
	foreach my $file ( @files ){
		next if $file =~ m/^\.+$/;
		$node = $currentdir.$file if $file =~ m/^node_/;
		$edge = $currentdir.$file if $file =~ m/^edge_/;
		$source = $currentdir.$file if $file =~ m/^source_/;
		$dest = $currentdir.$file if $file =~ m/^dest_/;
		
		if ($file =~ m/^source_(.*)(\.txt|tsv|dat)$/){	#shortestpath output filename (with path)
			$file =~ m/^source_(.*)(\.txt|tsv|dat)$/;
			$outfile_shortestpath = $outdir . "spath_" . $1 . $filesuffix . ".dat";
		}
	}
	closedir $dir_drug;
	$spathcmd .= " -n $node";
	$spathcmd .= " -e $edge";
	$spathcmd .= " -s $source";
	$spathcmd .= " -d $dest";
	$spathcmd .= " >> $outfile_shortestpath";
	system($spathcmd);
#for one drug-------------------------------------
}
exit;
