#!/usr/bin/perl
use warnings;
use strict;

#-----------------------------------------------------------------------------------------------------------------
# Purpose	: To convert ENSP id to Gene Symbols (or orf ID if not found) and score to distance
# Input		: StringDB file, UniProt ID map file, orfmap file
# Output	: Text, same formatted, Gene Symbol (or orf ID) and distance
#		  Distance = (1000-score)/1000
# Author	: Hansaim Lim
# Date		: 21 Aug, 2014
#-----------------------------------------------------------------------------------------------------------------
die "Usage: $0 <String-db file> <UniProt map file> <orf map file>\n" unless @ARGV==3;
my $string = shift @ARGV;
my $uniprot_idmap = shift @ARGV;
my $orfmap = shift @ARGV;
my $outfile = $string;
$outfile =~ s/^(.*)(\..*)$/$1_GS_dist$2/;

open my $UNIPROT, '<', $uniprot_idmap or die "Could not open file $uniprot_idmap: $!\n";
my %ensp_uniprots = ();
my %ensp_symbols = ();
while(<$UNIPROT>){
	next unless $_ =~ m/ENSP\d+/;
	my @ids = split(/\t/, $_);
	my $u = shift @ids;	#UniProtKB
	my $gs = shift @ids;	#gene symbol (with _HUMAN suffix)
	$gs =~ s/(.*)_HUMAN/$1/;#remove suffix

	foreach my $id ( @ids ){
		my @ensps = split(/; /, $id) if $id =~ m/ENSP\d+/;
		foreach my $ensp ( @ensps ){
			$ensp_symbols{$ensp} = $gs;
			$ensp_uniprots{$ensp} = $u;
		}
	}
}
close $UNIPROT;


open my $ORFMAP, '<', $orfmap or die "Could not open file $orfmap: $!\n";
my %ensp_orfs = ();
my $orfmap_line = 1;	#to skip the first line
ORFMAP: while(<$ORFMAP>){
	if ($orfmap_line == 1){
		$orfmap_line++;
		next ORFMAP;
	}
	my @ids = split(/\t/, $_);
	my $orf = shift @ids;
	my @ensps = shift @ids;
	my @uniprots = shift @ids;

	if ($_ =~ m/ENSP\d+/){ #if ENSP id found
		foreach my $ensp (@ensps){
			$ensp_orfs{$ensp} = $orf;
		}
	} else {	#try to match using UniProtKB when the orfmap line does not contain ENSP information
		my @ENSPs;
		foreach my $uniprot (@uniprots){
			foreach my $ENSP (keys %ensp_uniprots){
				push @ENSPs, $ENSP if $ensp_uniprots{$ENSP} eq $uniprot;	#push the ensp ids for given uniprot
			}
		}
		foreach my $ensp_additional (@ENSPs){
			$ensp_orfs{$ensp_additional} = $orf;	#ensp-orf pairs for orfs without ensp listed on orfmap file
		}
	}
}
close $ORFMAP;

open my $STRING, '<', $string or die "Could not open file $string: $!\n";
open my $OUTPUT, '>', $outfile or die "Could not open file $outfile: $!\n";
my $line_num = 1;	#to skip the first line
my $current_ensp1;
my $current_uniprot1;
my $lineskip = 0;	#count the number of skipped lines
String: while(<$STRING>){
	if ($line_num == 1){	#the first line is skipped
		print $OUTPUT "GeneSymbol1\tGeneSymbol22\tDistance";
		$line_num++;
		next String;
	}
	my $line = $_;
	my $line_intact = $_;	#this variable will hold the line as it is
	my @words = split(/ /, $line);	#[0] = protein1; [1] = protein2; [2] = score
	my $ensp1 = $words[0];
	my $ensp2 = $words[1];
	my $score = $words[2];

	$ensp1 =~ s/9606\.(ENSP\d{11})\D*/$1/;
	$ensp2 =~ s/9606\.(ENSP\d{11})\D*/$1/;
;
	my ($gs1, $gs2, $distance);
	$gs1 = $ensp_symbols{$ensp1} or $gs1 = $ensp_orfs{$ensp1};	#put gene symbol if gene symbol found, otherwise put orf
	$gs2 = $ensp_symbols{$ensp2} or $gs2 = $ensp_orfs{$ensp2};
	if (!$gs1 or !$gs2) {	#if anyone missing even after considering orfmaps
		$lineskip++;
		print $line_intact;
		next String;
	}
	chomp($score);
	$distance = ( (1000-$score)/1000 );	#score conversion to distance; higher score, closer distance
	print $OUTPUT "\n", $gs1, "\t", $gs2, "\t", $distance;	
}
close $OUTPUT;
close $STRING;
print $lineskip, " lines were skipped because any one of two proteins were not found in either UniProt map or Orfmap\n";
exit;
