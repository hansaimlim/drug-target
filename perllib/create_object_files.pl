#!/usr/bin/perl
use DrugBank;
use cMap;
use DrugTargetBase;
use strict;
use warnings;

my $cmap_obj = new cMap("100up");
store_hash($cmap_obj, "./static/json/cMap/100up.json");	#store cmap 100 up-regulated gene data in a file

my $db_obj = new DrugBank();
store_hash($db_obj, "./static/json/DrugBank/DrugBank.json");	#store DrugBank data in a file

exit;
