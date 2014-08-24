#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Purpose	: To get destination genes from  extreme-N-genes.txt file from cMap:drugRL data
#		  The extreme-N-genes.txt contains drugnames and gene symbols that are either highly upregulated
#		  or downregulated. (top 100 ranked OR bottom 100 ranked)
# Requirement	: The extreme-N-genes.txt file should start with an 'x' before drugname;
#		  After 'x' each column shows gene symbol and its rank
# Input		: extreme-N-genes.txt file; source-InChIKey file
# Output	: one text file per one drug; containing geme symbols and rank
#		  One text file reporting the number of shared targets for each drug
# Author	: Hansaim Lim
# Date		: 23 Aug, 2014
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
die "Usage: $0 <extreme-N-genes.txt> <cMap common drugname list> <output folder>\n" unless @ARGV == 3;
my $extreme_N_genes_file = shift @ARGV;
my $cMap_common_drugs_file = shift @ARGV;
my $output_folder = shift @ARGV;
$output_folder = dirname_endslash($output_folder);
make_dir($output_folder);

my %common_drugnames = ();
open my $cMap_commondrugs, '<', $cMap_common_drugs_file or die "Could not open file $cMap_common_drugs_file: $!\n";
while(<$cMap_commondrugs>){
	return if $_ =~ m/^\s*$/;	#skip empty line
	my @words = split(/\t/,$_);
	my $drug = shift @words;
	my $ikey = shift @words;
	$common_drugnames{$drug} = 1;	#like switches
}
close $cMap_commondrugs;

open my $CMAP, '<', $extreme_N_genes_file or die "Could not open file $extreme_N_genes_file: $!\n";
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
	if ($_ =~ m/^x\s*$/i){# the previous filehandler should be closed here
		close $outfilehandle if $is_drug_shared;	#close output filehandler if it was opened
		$is_drugname = 1;	#the next line contains drugname
		$is_drug_shared = 0;	#switch on when the drug is shared
		next LINE;
	}
	my $line_intact = $_;	#unmodified line
	my @words = split(/\t/, $_);
	if ($is_drugname){
		$drug = $words[1];
		chomp($drug);
		if ($common_drugnames{$drug}){
			$is_drug_shared = 1;	#valid until find next 'x'
			#modifying drug names
			my $drugname = $drug;
			$drugname =~ s/(\-)(\d|\w)/_$2/;	#dash to underscore (not minus sign for stereochemistry)
			$drugname =~ tr/()/__/;			#remove parenthesis
			$drugname =~ tr/+-/pm/;			#+ -> p; - -> m (stereochemistry, plus or minus)
			$drugname =~ s/\\/___/;			#backslash -> ___ (three underscores)
			$drugname =~ s/\//__/;			#slash -> __ (two underscores)
			my $fileout = $output_folder. "dest_" . $drugname . ".txt";
			open $outfilehandle, '>', $fileout or die "Could not open file $fileout: $!\n";
		}
		$is_drugname = 0;	#the next line is NOT drugname
	} elsif ($is_drug_shared){
		print $outfilehandle $line_intact;
	}
}
close $CMAP;

sub make_dir {
        my $dir_to_create = shift @_;
        unless(-e $dir_to_create or mkdir $dir_to_create){
                die "Unable to create $dir_to_create\n";
        }
        return;
}
sub dirname_endslash {
        my $dir_to_slash = shift @_;
        $dir_to_slash .= '/' unless $dir_to_slash =~ m/\/$/;    #put the end slash if missing
        return $dir_to_slash;
}

exit;
