#!/usr/bin/perl
use strict;
use warnings;
use DrugTargetBase;

#----------------------------------------------------------------------------------------------------------------------------------------------------
# Perform Dijkstra's shortestpath for each drug in the input folder
# Output the shortest pathways to output folders
# Output the distances for a whole set of drugs (up, down, random are separated)
#----------------------------------------------------------------------------------------------------------------------------------------------------
die "Usage: $0 <input folder> <output folder> <distance output filename>\n" unless @ARGV == 3;
my $input_folder = dirname_add_slash(shift @ARGV);
my $output_folder = dirname_add_slash(shift @ARGV);
my $distance_file = shift @ARGV;
chomp($distance_file);

my $edge_file = "./static/edge/9606.protein.links.v9.1-GN-dist.txt";
my $node_file = "./static/node/node.txt";
make_dir($output_folder);

open my $Distance_file, '>', $distance_file or die "Could not open distance file $distance_file: $!\n";

my @drug_dirs = read_directory($input_folder);	#sub directories with drug names
foreach my $drug (@drug_dirs){
	next if $drug =~ m/^\.+/;
	my $sub_dir = $input_folder . $drug;
	my @files = read_directory($sub_dir);
	my $source_file;
	my $dest_file;
	foreach my $file (@files){
		$source_file = $file if $file =~ m/source/;
		$dest_file = $file if $file =~ m/dest/;
	}
	my $spath_outfile = $output_folder . "spath_" . $drug . ".txt";
        my $spathcmd = 'shortestpath';
        $spathcmd .= " -n $node_file";
        $spathcmd .= " -e $edge_file";
        $spathcmd .= " -s $source_file";
        $spathcmd .= " -d $dest_file";
        $spathcmd .= " >> $spath_outfile";
        system($spathcmd);      #this command produces shortestpath output

	open my $SPATH, '<', $spath_outfile or die "Could not open shortest path output file $spath_outfile: $!\n";
	while (my $line = <$SPATH>){
		next if $line =~ m/^the\snumber\sof/i;
                next if $line =~ m/^0\.\..*0$/;
                next if $line =~ m/^\->1e\+07/i;
		next if $line =~ m/^No\spath\sfrom/i;
		$line =~ m/\->(.+)$/;
		my $distance = $1;
		print $Distance_file $distance, "\n";
	}
	close $SPATH;
}

close $Distance_file;

sub read_directory
{
	my $dir = shift @_;
	chomp($dir);
	opendir my $Directory, $dir or die "Could not open directory $dir: $!\n";
	my @files = readdir $Directory;
	close $Directory;
	return @files;
}
