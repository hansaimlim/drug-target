#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Data::Dumper;

die "Usage: $0 <username> <password>\n" unless @ARGV == 2;

my $user = shift @ARGV;
my $password = shift @ARGV;
chomp($password);
my $dbh = connect_chembl_mysql($user, $password);
my @multi_domain_targets = get_accession_multidomain($dbh);	#Targets in this array appear to have multiple binding domains; 37 accessions for chembl_19
my $zincfile = "./zinc/idmap/zinc_idmap.tsv";
open my $ZINC, '<', $zincfile or die "Could not open zinc idmap file $zincfile: $!\n";
while (my $line = <$ZINC>){
	next if $. == 1;	#skip first line
	my @words = split(/\t/, $line);
	my $genename = shift @words;
	my $accession = shift @words;
	my $chembl = shift @words;
	my $annotation = shift @words;
	chomp($annotation);
	
	my $start = 71;
	my $end = 380;
	my $subseq = get_substring($dbh, $accession, $start, $end);
}
close $ZINC;
sub get_pdb {

}
sub connect_chembl_mysql {
	my $user = shift @_;
	my $pw = shift @_;
	chomp($user);
	chomp($pw);
	my $driver = "mysql"; 
	my $database = "chembl_19";
	my $dsn = "DBI:$driver:database=$database";
	my $dbh = DBI->connect($dsn, $user, $pw) or die $DBI::errstr;
	return $dbh;
}

sub get_substring {
	my ($dbh, $accession, $start, $end) = @_;
	my $sth = $dbh->prepare("SELECT SUBSTRING(sequence, $start, ($end-$start+1)) FROM component_sequences WHERE accession = '$accession'");
	$sth->execute() or die $DBI::errstr;
	my $num_row = $sth->rows;
	my $sequence = $sth->fetchrow_array();
	$sth->finish();
	return $sequence;
}

sub get_num_site {

}

sub get_accession_multidomain {
	my $dbh = shift @_;
	my $sth = $dbh->prepare("SELECT cs.accession FROM component_sequences cs INNER JOIN component_domains cd ON cs.component_id = cd.component_id INNER JOIN site_components sc ON cs.component_id = sc.component_id AND cd.domain_id = sc.domain_id GROUP BY cs.accession HAVING (COUNT(DISTINCT(cd.compd_id)) > 1 AND COUNT(DISTINCT(cd.domain_id)) > 1)");
	$sth->execute() or die $DBI::errstr;
	my @accessions;
	while(my $accession = $sth->fetchrow_array()){
		chomp($accession);
		push @accessions, $accession;
	}
	$sth->finish();
	return @accessions;
}

sub make_dir
{
        #create a directory unless already created
        my $dir = shift @_;
        unless (-e $dir or mkdir $dir){
                die "Unable to create $dir\n";
        }
        return;
}

sub dirname_add_slash
{
        #add a slash at the end of a directory name unless already attached
        my $dir = shift @_;
        $dir .= '/' unless $dir =~ m/\/$/;
        return $dir;
}

exit;
