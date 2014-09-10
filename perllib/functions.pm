#!/usr/bin/perl

#--------------------------------------------------
# Collection of functions that are used for
# drug-target network analysis
# A preparation for creating a set of modules
#--------------------------------------------------

package FunctionsUnsorted;
use strict;

#------------TEST AREA-----------------

#------------TEST AREA-----------------

sub get_GeneSymbol_by_ENSP
{
	
}


sub rm_special_char_in_drugname
{
	my $drug = shift @_;
	$drug =~ s/\///g;	#remove slash
	$drug =~ s/(\-)(\d|\w)/_$2/;	#a dash to an underscore (but not the stereochemical minus sign)
	$drug =~ tr/()/__/;	#parenthesis to underscores
	$drug =~ tr/+-/pm/;	#stereochemical signs to letters
	$drug =~ tr/\\//;	#remove backslash -- just in case
	return $drug;
}
