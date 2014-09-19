#!/usr/bin/perl
use DrugBank;
use DrugTargetBase;
use strict;
use warnings;

my $drugbank = new DrugBank();
my $file = "./static/idmap/DrugBank_InChIKey_drugname.txt";
open my $FH, '>', $file or die "Could not open drugbank outfile $file: $!\n";
foreach my $ikey (keys %$drugbank){
	my $drugname = $drugbank->get_DrugBank_drugname_by_InChIKey($ikey);
	print $FH $ikey, "\t", $drugname, "\n";
}
close $FH;
exit;
