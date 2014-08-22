#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;

#----------------------------------------------------------------------------------------------------------
# Purpose	: To process extreme-N-genes.txt file from cMap:drugRL data
#		  The extreme-N-genes.txt contains drugnames and gene symbols that are either highly upregulated
#		  or downregulated. (top 100 ranked OR bottom 100 ranked)
# Requirement	: The extreme-N-genes.txt file should start with an 'x' before drugname;
#		  After 'x' each column shows gene symbol and its rank
# Input		: extreme-N-genes.txt file; source-InChIKey file
# Output	: one text file per one drug; containing geme symbols and rank
#		  One text file reporting the number of shared targets for each drug
# Author	: Hansaim Lim
# Date		: 14 Jul, 2014
#----------------------------------------------------------------------------------------------------------
die "Usage: $0 <extreme-N-genes.txt> <source-InChIKey>\n" unless @ARGV == 2;
my $cmap = shift @ARGV;
my $source = shift @ARGV;

my %shared = ();
open my $src, '<', $source or die "Could not open file $source: $!\n";
while(<$src>){
	return if $_ =~ m/^\s*$/;	#skip empty line
	my $ikey = $_;
	chomp($ikey);
	$shared{$ikey} = 1;	#like switches
}
close $src;

open my $CMAP, '<', $cmap or die "Could not open file $cmap: $!\n";
my $is_drugname = 0;
my $is_drug_shared = 0;
my $drug;
my $line = 1;	#to skip first line ('x')
my $outfilehandle;
LINE: while(<$CMAP>){
	if ($line == 1){ 
		$line++;
		$is_drugname = 1;
		next LINE;
	}
	my $line_intact = $_;	#unmodified line
	my @words = split(/\t/, $_);
	my $symbol = $words[0];
	my $rank;
	my $drugname;
	if ($_ =~ m/^x\s*$/i){# the previous filehandler should be closed here
		close $outfilehandle;
		$is_drugname = 1;	#the next line contains drugname
		$is_drug_shared = 0;	#switch on when the drug is shared
		next LINE;
	}
	if ($is_drugname){
		$drug = $words[1];
		chomp($drug);
		#get InChIKey of the drug and compare if shared
		my $inchikey = pubchem_inchikey_by_drug($drug);
		if ($shared{$inchikey}){
			$is_drug_shared = 1;	#valid until find next 'x'
		}
		#modifying drug names
		$drugname = $drug;
		$drugname =~ s/(\-)(\d|\w)/_$2/;	#dash to underscore (not minus sign for stereochemistry)
		$drugname =~ tr/()/__/;			#remove parenthesis
		$drugname =~ tr/+-/pm/;			#+ -> p; - -> m (stereochemistry, plus or minus)
		$drugname =~ s/\\/___/;			#backslash -> ___ (three underscores)
		$drugname =~ s/\//__/;			#slash -> __ (two underscores)
		my $fileout = "dest_" . $drugname . ".txt";
		open $outfilehandle, '>', $fileout or die "Could not open file $fileout: $!\n";
		$is_drugname = 0;	#the next line is NOT drugname
	} elsif ($is_drug_shared){
		print $outfilehandle $line_intact;
	}
}
close $CMAP;

sub pubchem_inchikey_by_drug {
	my $drug = shift @_;
	my $url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$drug/property/InChIKey/TXT";
	my $inchikey = get $url;
	chomp($inchikey);
	return $inchikey;
}
exit;
