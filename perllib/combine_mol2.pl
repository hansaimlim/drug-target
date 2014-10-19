#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

die "Usage: $0 <mol2 parent directory> <outfile>\n" unless @ARGV == 2;
my $parent_dir = shift @ARGV;
my $outfile = shift @ARGV;
chomp($outfile);
my @files = read_directory($parent_dir);

my %mol2_data;	#keys are atom and bond information; values are molecule information
foreach my $file (@files){
	next if $file =~ m/^\.+/;
	my $mol2file = $parent_dir.$file;
	open my $MOL2, '<', $mol2file or die "Could not open mol2 file $file: $!\n";
	my $datatype = 0;	#1->MOLECULE; 2->ATOM; 3->BOND
	my $atom;	#long string that may contain new line characters
	my $bond;
	my $molecule;
	while (my $line = <$MOL2>){
		$line =~ s/[ \t]+/\t/g;
		if ($line =~ m/^@<TRIPOS>MOLECULE/i){
			$datatype = 1;
			my $atom_bond = $atom.$bond;
			$mol2_data{$atom_bond} = $molecule unless $. == 1;	#avoid empty lines; this line may produce warnings;
			undef $molecule;
			undef $atom;
			undef $bond;
		}
		elsif ($line =~ m/^@<TRIPOS>ATOM/i){
			$datatype = 2;
		}
		elsif ($line =~ m/^@<TRIPOS>BOND/i){
			$datatype = 3;
		}

		$molecule .= $line if $datatype == 1;
		$atom .= $line if $datatype == 2;
		$bond .= $line if $datatype == 3;
		
	}
	my $last_atom_bond = $atom.$bond;
	$mol2_data{$last_atom_bond} = $molecule;	#include the last molecule data;
	close $MOL2;
}

open my $MOL2_OUT, '>', $outfile or die "Could not open output mol2 file $outfile: $!\n";
foreach my $atombond (keys %mol2_data){
	my $molecule = $mol2_data{$atombond};
	print $MOL2_OUT $molecule, $atombond;
}
close $MOL2_OUT;
sub read_directory
{
        my $dir = shift @_;
        chomp($dir);
        opendir my $Directory, $dir or die "Could not open directory $dir: $!\n";
        my @files = readdir $Directory;
        close $Directory;
        return @files;
}
exit;
