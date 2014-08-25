#!/usr/bin/perl
use strict;
use warnings;
#---------------------------------------------------------------------------------------------------------------------------------------------
# Purpose	: Automated work from prepared destination files to shortestpath output and to shortestpath analyzer 
# Input		: source (STDBsource), destfile directory, String proteinlink_maxdistXX file, output dir for network files
# Output	: four input files; shortestpath output files; shortestpath analysis files
# Author	: Hansaim Lim
# Date		: 19 Aug, 2014
#---------------------------------------------------------------------------------------------------------------------------------------------
die "Usage: $0 <source.txt> <directory for dest files> <String-distance file> <foldername: networks_XXupmaxXX (do not add ./)>\n" unless @ARGV == 4;
my $source = shift @ARGV;
my $destdir = shift @ARGV;
my $string = shift @ARGV;
my $network_upper_dir = shift @ARGV;
chomp($network_upper_dir);
$network_upper_dir = dirname_endslash($network_upper_dir);
make_dir($network_upper_dir);
#-------------------------------------Creation of the input files (network files; node, edge, source, dest) started-------------------------------------
my %sources = ();
open my $src, '<', $source or die "Could not open file $source: $!\n";
while(<$src>){
        my $gene = $_;
        chomp($gene);
        $sources{$gene} = 1;    #switches
}
close $src;

opendir my $dir, $destdir or die "Could not open destdir directory $destdir: $!\n";
my @destfiles = readdir $dir;
foreach my $dest ( @destfiles ){
	next if $dest =~ m/^\.+$/;
	my ($network_lower_dir, $sourceoutput, $destoutput, $edgeoutput, $nodeoutput);
	my ($outfile_edge_temp, $tempfile1);
	TEMPFILE:{	#every temp file should be different for multi-task
		$network_upper_dir =~ m/([a-zA-Z0-9.]+)\/$/;
		$outfile_edge_temp = "./.edge$1temp.txt";
		$tempfile1 = "./.spath$1temp.txt";
	}
	OUTFILE2:{
		$dest =~ m/dest_(.*)(\.txt|tsv|dat)$/;
		$network_lower_dir = $network_upper_dir . $1;
		$edgeoutput = "./" .$network_lower_dir. "/edge_" . $1 . $2;
		$nodeoutput = "./" . $network_lower_dir. "/node_" . $1 . $2;
		$sourceoutput = "./" .$network_lower_dir. "/source_" . $1. $2;
		$destoutput = "./" . $network_lower_dir. "/dest_" . $1 . $2;
	}
	make_dir($network_lower_dir);
	my $dest_with_path = $destdir . $dest;
	open my $destination, '<', $dest_with_path or die "Could not open dest file $dest: $!\n";
	my %destinations = ();
	while(<$destination>){
		next if $_ =~ m/^x\s*$/;	#skip the line with only 'x' and/or white spaces
		my @words = split(/\t/, $_);
		my $gene = shift @words;
		chomp($gene);
		$destinations{$gene} = 1;       #switches
	}
	close $destination;

	open my $string, '<', $string or die "Could not open string file $string: $!\n";
	open my $outfile, '>', $outfile_edge_temp or die "Could not open outfile_edge_temp file $outfile_edge_temp: $!\n";
	my $line = 1;   #to skip first line
	STRING: while(<$string>){
		if ($line == 1){
			$line++;
			next STRING;
		}
		my $line_intact = $_;   #the whole line
		my @words = split(/\t/, $_);
		my $gene1 = shift @words;
		my $gene2 = shift @words;

		if ($sources{$gene1}) { $sources{$gene1} = "found"; }
		if ($sources{$gene2}) { $sources{$gene2} = "found"; }
		if ($destinations{$gene1}) { $destinations{$gene1} = "found"; }
		if ($destinations{$gene2}) { $destinations{$gene2} = "found"; }
		if ($sources{$gene1} || $sources{$gene2} || $destinations{$gene1} || $destinations{$gene2}){
			print $outfile $line_intact;
		}
	}
	close $outfile;
	close $string;

	open my $sourceout, '>', $sourceoutput or die "Could not open sourceoutput file $sourceoutput: $!\n";
	foreach my $source ( keys %sources ){
		if ($sources{$source} eq "found"){
			print $sourceout $source, "\n";
		}
	}
	close $sourceout;
	open my $destout, '>', $destoutput or die "Could not open destoutput file $destoutput: $!\n";
	foreach my $dest ( keys %destinations ){
		if ($destinations{$dest} eq "found"){
			print $destout $dest, "\n";
		}
	}
	close $destout;

	# Source and destination files are done at this point
	# Only edge files will be re-opened for the steps below

	my %ppi = (); #A->B = c hash container
	my $skipped = 0;
	my $edgenum = 0;	#to count Number of Edges:
	open my $edge_temp, '<', $outfile_edge_temp  or die "Could not open outfile_edge_temp file $outfile_edge_temp: $!\n";
	open my $temp1, '>', $tempfile1 or die "Could not open tempfile1 file $tempfile1: $!\n";
	EDGE: while(<$edge_temp>){
		my @words = split(/\t/, $_);
		if ($words[0] =~ m/^Number/i){
			next EDGE;
		}
		my $p1 = shift @words;  #protein 1
		my $p2 = shift @words;  #protein 2
		my $d = shift @words;   #distance
		chomp($d);
		unless ($ppi{$p1}{$p2}){ #if the relationship is undefined yet
			unless ($ppi{$p2}{$p1}) {       #if the reverse relationship is undefined as well
				print $temp1 $p1,"\t",$p2,"\t",$d,"\n";
				$ppi{$p1}{$p2} = $d;
				$edgenum++;
				next EDGE;
			}
			#when the reverse is already defined
			$skipped++;
			next EDGE;
		}
	}
	close $edge_temp;
	close $temp1;
	print $edgenum, "edges were found for $dest\n";
	print $skipped, "edges were skipped for $dest\n";

	#--------------------create node file from temporary edge file-----------------------
	my %nodes = ();
	open my $edgetemp, '<', $tempfile1 or die "Could not open file tempfile1 $tempfile1: $!\n";
	while(<$edgetemp>){
		my @words = split(/\t/, $_);
		my $g1 = shift @words;
		my $g2 = shift @words;
		$nodes{$g1} = 1;
		$nodes{$g2} = 1;
	}
	close $edgetemp;

	open my $nodetemp, '>', $nodeoutput or die "Could not open file nodeoutput $nodeoutput: $!\n";
	my $nodenum = keys %nodes;
	print $nodetemp "Number of Nodes: ", $nodenum, "\n";
	foreach my $node ( keys %nodes ){
		print $nodetemp $node, "\t1\n";	#1 is for fold change
	}
	close $nodetemp;
	#--------------------create node file from temporary edge file-----------------------

	#--------------------add the header, Number of Edges: xxxx------------------
	open my $temp2, '<', $tempfile1 or die "Could not open tempfile1 file $tempfile1: $!\n";
	open my $edge_final, '>', $edgeoutput or die "Could not open edgeoutput file $edgeoutput: $!\n";
	print $edge_final "Number of Edges: ", $edgenum, "\n";
	while(<$temp2>){
		print $edge_final $_;
	}
	close $edge_final;
	close $temp2;
	#--------------------add the header, Number of Edges: xxxx------------------
}
closedir $dir;
#-------------------------------------Creation of the input files (network files; node, edge, source, dest) completed-------------------------------------
#-------------------------------------Running shortestpath algorithm starts here--------------------------------------------------------------------------

my $suffix_for_spath;	#file suffix for shortestpath algorithm
my $outdir_for_spath;	#output directory for shortestpath algorithm
SUFFIX: {
	$network_upper_dir =~ m/networks_(.*)\/$/;	#the suffix should NOT have an end slash
	$suffix_for_spath = $1;
}
$network_upper_dir = dirname_endslash($network_upper_dir);
OUTDIR_SPATH: {
	$network_upper_dir =~ m/([a-zA-Z0-9.]+)\/$/;
	$outdir_for_spath = "./spath_" . $1 . "/"; 
}

make_dir($outdir_for_spath);
opendir my $dir1, $network_upper_dir or die "Cannot open network_upper_dir directory $network_upper_dir: $!\n";
my @lowerdirs = readdir $dir1;
closedir $dir1;

LOWERDIR: foreach my $lowerdir ( @lowerdirs ){
	next LOWERDIR if $lowerdir =~ m/^\.+$/;	#skip current and previous dir (. and ..)
	my $currentdir = $network_upper_dir . $lowerdir."/";
	my ($node, $edge, $source, $dest, $outfile_shortestpath);
#---------------------------for one drug-----run shortestpath and analyze the output---------------------------------------------------------
	opendir my $dir_drug, $currentdir or die "Cannot open currentdir directory $currentdir: $!\n";
	my @files = readdir $dir_drug;
	foreach my $file ( @files ){
		next if $file =~ m/^\.+$/;
		$node = $currentdir.$file if $file =~ m/^node_/;
		$edge = $currentdir.$file if $file =~ m/^edge_/;
		$source = $currentdir.$file if $file =~ m/^source_/;
		$dest = $currentdir.$file if $file =~ m/^dest_/;
		
		if ($file =~ m/^source_(.*)(\.txt|tsv|dat)$/){	#shortestpath output filename (with path)
			$file =~ m/^source_(.*)(\.txt|tsv|dat)$/;
			$outfile_shortestpath = $outdir_for_spath . "spath_" . $1 . $suffix_for_spath . ".dat";
		}
	}
	closedir $dir_drug;
	my $spathcmd = 'shortestpath';
	$spathcmd .= " -n $node";
	$spathcmd .= " -e $edge";
	$spathcmd .= " -s $source";
	$spathcmd .= " -d $dest";
	$spathcmd .= " >> $outfile_shortestpath";
	system($spathcmd);	#this command produces shortestpath output
}
#---------------------------for one drug-----run shortestpath and analyze the output---------------------------------------------------------
#-------------------------------------Running shortestpath algorithm completed--------------------------------------------------------------------------
#-------------------------------------Analyzing shortestpath output starts here-------------------------------------------------------------------------
my $spath_dir = $outdir_for_spath;
my $spath_analysis_outdir;
my ($pathwayproperties, $pathwaydisconnect);
SPATHANALYSIS: {
	$network_upper_dir =~ m/networks_(.*)$/;
	$spath_analysis_outdir = "./spathanalysis_" . $1;
}
$spath_analysis_outdir = dirname_endslash($spath_analysis_outdir);
opendir my $dir_spath, $spath_dir or die "Could not open spath_dir directory $spath_dir: $!\n";
my @spathfiles = readdir $dir_spath;
closedir $dir_spath;

make_dir($spath_analysis_outdir);

for my $spathfile ( @spathfiles ){
        next if $spathfile =~ m/^\.+$/; #skip currentdir and upper dir
        OUTPUTFILE: {
                $spathfile =~ m/spath_(.*)\.(txt|tsv|dat)$/;
                $pathwayproperties = $spath_analysis_outdir . $1 . "_property." . $2;
                $pathwaydisconnect = $spath_analysis_outdir . $1 . "_disconnect." . $2;
                $spathfile = $spath_dir.$spathfile;
        }


        open my $spathfileinput, '<', $spathfile or die "Could not open spathfile file $spathfile: $!\n";
        open my $spath_property, '>', $pathwayproperties or die "Could not open spathwayproperties file $pathwayproperties: $!\n";
        open my $spath_disconnect, '>', $pathwaydisconnect or die "Could not open spathwaydisconnect file $pathwaydisconnect: $!\n";
        print $spath_property "source\tdest\tedges\tcum_dist\tavg_dist_per_edge";
        print $spath_disconnect "source\tdestinations_disconnected_to";
        my $current_source_disc = '';   #current source for disconnected pathways

        SPATH_LINE1: while(<$spathfileinput>){
                my $line_intact = $_;
                next SPATH_LINE1 if $line_intact =~ m/^the\snumber\sof/i;
                next SPATH_LINE1 if $line_intact =~ m/^0\.\..*0$/;
                next SPATH_LINE1 if $line_intact =~ m/^\->1e\+07/i;

                if ($line_intact =~ m/^No\spath\sfrom/i){
                        #organize and output for disconnected pathways
                        $line_intact =~ m/^No\spath\sfrom\s(.+)\sto\s(.+)$/;
                        my $source_disc = $1;
                        my $dest_disc = $2;
                        if ($source_disc eq $current_source_disc) {
                                #no need to change the line of output
                                print $spath_disconnect "\t", $dest_disc;
                        } else {
                                $current_source_disc = $source_disc;    #update current source of disconnection
                                print $spath_disconnect "\n";
                                print $spath_disconnect $source_disc, "\t", $dest_disc;
                        }
                        next SPATH_LINE1;
                }
                my @words = split(/\.\./, $_);
                my $last_index = scalar(@words) - 1;    #equals to pathway size (number of nodes in the pathway)
                my $dest_index = $last_index - 1;       #index for destination; equal to the number of edges
                chomp($words[$last_index]);
                my $cum_dist = $words[$last_index];
                DISTANCE: {
                        $cum_dist =~ m/\->(.+)$/;
                        $cum_dist = $1;
                }
                my $avg_dist_per_edge = $cum_dist / $dest_index;
                print $spath_property "\n", $words[0], "\t", $words[$dest_index], "\t", $dest_index, "\t", $cum_dist, "\t", $avg_dist_per_edge;
        }

        close $spath_disconnect;
        close $spath_property;
        close $spathfileinput;
}
#-------------------------------------Analyzing shortestpath output completed---------------------------------------------------------------------------
#-------------------------------------Subroutines-------------------------------------------------------------------------------------------------------
sub make_dir {
	my $dir_to_create = shift @_;	
	unless(-e $dir_to_create or mkdir $dir_to_create){
		die "Unable to create $dir_to_create\n";
	}
	return;
}
sub dirname_endslash {
	my $dir_to_slash = shift @_;
	$dir_to_slash .= '/' unless $dir_to_slash =~ m/\/$/;	#put the end slash if missing
	return $dir_to_slash;
}
#-------------------------------------Subroutines-------------------------------------------------------------------------------------------------------



exit;
