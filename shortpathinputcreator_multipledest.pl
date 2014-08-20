#!/usr/bin/perl
use strict;
use warnings;
#---------------------------------------------------------------------------------------------------------------
# Purpose	: To combine edge.pl, edgeclear.pl, tempTonode.pl as well as commands to create temp (node) file
# Input		: same as edge.pl
# Output	: four input files for shortestpath algorithm
# Author	: Hansaim Lim
# Date		: 13 Aug, 2014
#---------------------------------------------------------------------------------------------------------------
die "Usage: $0 <source.txt> <directory for dest files> <String-UniProt-distance file>\n" unless @ARGV == 3;
my $source = shift @ARGV;
my $destdir = shift @ARGV;
my $string = shift @ARGV;


my %sources = ();
open(SRC, "<$source") or die "Could not open file $source: $!\n";
while(<SRC>){
        my $gene = $_;
        chomp($gene);
        $sources{$gene} = 1;    #switches
}
close(SRC);

opendir my $dir, $destdir or die "Could not open directory $destdir: $!\n";
my @destfiles = readdir $dir;
for my $dest ( @destfiles ){
	my %destinations = ();
	next if $dest =~ m/^\.+$/;
	my $outfile_edge_temp = $dest;
	my ($subdirectory, $outfile_source, $outfile_dest, $tempfile1, $tempfile2, $edgeoutput, $nodeoutput);
	OUTFILE1: {
        $outfile_edge_temp = "./.edgetemp.txt";
        $tempfile1 = "./.spathtemp1.txt";
        $tempfile2 = "./.spathtemp2.txt";
	}
	OUTFILE2:{
		$dest =~ m/dest(.*)(\.txt|tsv|dat)$/;
		$subdirectory = "dir" . $1;
		$edgeoutput = "./" .$subdirectory. "/edge" . $1 . $2;
		$nodeoutput = "./" . $subdirectory. "/node" . $1 . $2;
		$outfile_source = "./" .$subdirectory. "/source" . $1. $2;
		$outfile_dest = "./" . $subdirectory. "/dest" . $1 . $2;
	}
	unless(-e $subdirectory or mkdir $subdirectory){
		die "Unable to create $subdirectory\n";
	}
	my $dest_with_path = $destdir . $dest;
	open(DEST, "<$dest_with_path") or die "Could not open file $dest: $!\n";
	while(<DEST>){
		my $gene = $_;
		chomp($gene);
		$destinations{$gene} = 1;       #switches
	}
	close(DEST);
	open(STRING, "<$string") or die "Could not open file $string: $!\n";
	open(OUTFILE, ">$outfile_edge_temp") or die "Could not open file $outfile_edge_temp: $!\n";
	my $line = 1;   #to skip first line
	STRING: while(<STRING>){
		if ($line == 1){
			$line++;
			next STRING;
		}
		my $line_intact = $_;   #the whole line
		my @words = split(/\t/, $_);
		my $gene1 = shift @words;
		my $gene2 = shift @words;

		if ($sources{$gene1}) { $sources{$gene1} = "found"; }
		if ($sources{$gene2}) { $sources{$gene2} = "found"; }
		if ($destinations{$gene1}) { $destinations{$gene1} = "found"; }
		if ($destinations{$gene2}) { $destinations{$gene2} = "found"; }
		if ($sources{$gene1} || $sources{$gene2} || $destinations{$gene1} || $destinations{$gene2}){
			print OUTFILE $line_intact;
		}
	}
	close(OUTFILE);
	close(STRING);
	open(SOURCE, ">$outfile_source") or die "Could not open file $outfile_source: $!\n";
	foreach my $source ( keys %sources ){
		if ($sources{$source} eq "found"){
			print SOURCE $source, "\n";
		}
	}
	close(SOURCE);
	open(DEST, ">$outfile_dest") or die "Could not open file $outfile_dest: $!\n";
	foreach my $dest ( keys %destinations ){
		if ($destinations{$dest} eq "found"){
			print DEST $dest, "\n";
		}
	}
	close(DEST);

	# Source and destination files are done at this point
	# Only edge files will be re-opened for the steps below

	my %ppi = (); #A->B = c hash container
	my $skipped = 0;
	my $edgenum = 0;	#to count Number of Edges:
	open my $edge_temp, '<', $outfile_edge_temp  or die "Could not open file $outfile_edge_temp: $!\n";
	open my $temp1, '>', $tempfile1 or die "Could not open file $tempfile1: $!\n";
	EDGE: while(<$edge_temp>){
		my @words = split(/\t/, $_);
		if ($words[0] =~ m/^Number/i){
			next EDGE;
		}
		my $p1 = shift @words;  #protein 1
		my $p2 = shift @words;  #protein 2
		my $d = shift @words;   #distance
		chomp($d);
		unless ($ppi{$p1}{$p2}){ #if the relationship is undefined yet
			unless ($ppi{$p2}{$p1}) {       #if the reverse relationship is undefined as well
				print $temp1 $p1,"\t",$p2,"\t",$d,"\n";
				$ppi{$p1}{$p2} = $d;
				$edgenum++;
				next EDGE;
			}
			#when the reverse is already defined
			$skipped++;
			next EDGE;
		}
	}
	close $edge_temp;
	close $temp1;
	print $edgenum, "edges were found for $dest\n";
	print $skipped, "edges were skipped for $dest\n";

	#--------------------create node file from temporary edge file-----------------------
	my %nodes = ();
	open my $edgetemp, '<', $tempfile1 or die "Could not open file $tempfile1: $!\n";
	while(<$edgetemp>){
		my @words = split(/\t/, $_);
		my $g1 = shift @words;
		my $g2 = shift @words;
		$nodes{$g1} = 1;
		$nodes{$g2} = 1;
	}
	close $edgetemp;

	open my $nodetemp, '>', $nodeoutput or die "Could not open file $nodeoutput: $!\n";
	my $nodenum = keys %nodes;
	print $nodetemp "Number of Nodes: ", $nodenum, "\n";
	foreach my $node ( keys %nodes ){
		print $nodetemp $node, "\t1\n";	#1 is for fold change
	}
	close $nodetemp;
	#--------------------create node file from temporary edge file-----------------------

	#--------------------add the header, Number of Edges: xxxx------------------
	open my $temp2, '<', $tempfile1 or die "Could not open file $tempfile1: $!\n";
	open my $edge_final, '>', $edgeoutput or die "Could not open file $edgeoutput: $!\n";
	print $edge_final "Number of Edges: ", $edgenum, "\n";
	while(<$temp2>){
		print $edge_final $_;
	}
	close $edge_final;
	close $temp2;
	#--------------------add the header, Number of Edges: xxxx------------------
}
closedir $dir;
exit;
