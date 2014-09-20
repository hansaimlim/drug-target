#!/usr/bin/perl
use DrugBank;
use cMap;
use DrugTargetBase;
use strict;
use warnings;

create_cmap_json("100down");

sub create_cmap_json
{
	#input : range
	my $range = shift @_;
	chomp($range);
	my $file = "./static/json/cMap/" . $range . ".json";
	my $cmap = new cMap($range);
	store_hash($cmap, $file);
	return;
}
sub create_drugbank_json
{
	my $db_obj = new DrugBank();
	store_hash($db_obj, "./static/json/DrugBank/DrugBank.json");
	return;
}
exit;
