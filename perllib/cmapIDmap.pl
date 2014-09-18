#!/usr/bin/perl
use cMap;
use IDMAP;
use DrugTargetBase;
use strict;
use warnings;

my $cmap50up = new cMap("50up");      #the range (50up) does not matter for extraction of common drugs
my $ref = $cmap50up->get_cMap_InChIKey();
my @arr = @$ref;
my $file = "./static/idmap/cMap_InChIKey_drugname.tsv";
open my $FH, '>', $file or die "Could not open cmapfile $file: $!\n";
foreach my $ikey (@arr){
        print $FH $ikey, "\t";
        print $FH $cmap50up->get_cMap_drugname_by_InChIKey($ikey), "\n";
}
close $FH;
exit;
