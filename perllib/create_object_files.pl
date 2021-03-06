#!/usr/bin/perl
#use DrugBank;
#use cMap;
use String;
use DrugTargetBase;
use strict;
use warnings;

#PUGREST must be on in cMap.pm script
#create_cmap_json("100up");	#done
#create_cmap_json("100down");	#done
#create_cmap_json("100rand");	#done
#create_cmap_json("200up");	#done
#create_cmap_json("200down");	#done
#create_cmap_json("200rand");	#done
#create_cmap_json("1000up");	#done
#create_cmap_json("1000down");	#done
#create_cmap_json("1000rand");	#done
#print_random_ppi(20);	#testing the randomization mechanism
#print_random_ppi(30);
create_random_ppi_json(10);	#low shuffle numbers are for controlled experiments
create_random_ppi_json(20);
create_random_ppi_json(30);
create_random_ppi_json(100000);
create_random_ppi_json(1000000);
create_random_ppi_json(2000000);
create_random_ppi_json(3000000);
create_random_ppi_json(4000000);
create_random_ppi_json(5000000);
sub print_random_ppi
{
	my $n = shift @_;
	my $string_obj = new String();
	my $string_shuffled_n = $string_obj->Randomize_String($n);
	my $file = "./static/demo/randomppi".$n.".tsv";
	my %hash = %$string_shuffled_n;
	open my $FH, '>', $file or die "Could not open file $file: $!\n";
	print $FH "Gene1\tGene2\tDistance\n";
	foreach my $g1 (keys %hash){
		foreach my $g2 (keys $hash{$g1}){
			print $FH "$g1\t$g2\t$hash{$g1}{$g2}\n";
		}
	}
	close $FH;
}
sub create_cmap_json
{
	#input : range
	my $range = shift @_;
	chomp($range);
	my $file = "./static/json/cMap/" . $range . ".json";
	my $cmap = new cMap($range, "on");	#on means PUGREST on
	store_hash($cmap, $file);
	return;
}
sub create_drugbank_json
{
	my $db_obj = new DrugBank("on");	#on means PUGREST on
	store_hash($db_obj, "./static/json/DrugBank/DrugBank.json");
	return;
}
sub create_random_ppi_json
{
	my $n = shift @_;
	my $string_obj = new String();
	my $string_shuffle_n = $string_obj->Randomize_String($n);
	make_dir("./static/json/String/");
	my $file = "./static/json/String/Random_String" . $n . ".json";
	store_hash($string_shuffle_n, $file);
	return;
}
exit;
