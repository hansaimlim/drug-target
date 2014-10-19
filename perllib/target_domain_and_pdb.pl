#!/usr/bin/perl
use strict;
use warnings;
no warnings ('uninitialized', 'substr');
use DBI;
use List::Util qw( min max );
use LWP::Simple qw( $ua getstore );
use Data::Dumper;

die "Usage: $0 <username> <password>\n" unless @ARGV == 2;

my $user = shift @ARGV;
my $password = shift @ARGV;
chomp($password);
my $dbh = connect_chembl_mysql($user, $password);
#my @multi_domain_targets = get_accession_multidomain($dbh);	#Targets in this array appear to have multiple binding domains; 37 accessions for chembl_19
my $zincfile = "./static/idmap/zinc_idmap.tsv";
open my $ZINC, '<', $zincfile or die "Could not open zinc idmap file $zincfile: $!\n";
while (my $line = <$ZINC>){
	next if $. == 1;	#skip first line
	my @words = split(/\t/, $line);
	my $genename = shift @words;
	my $accession = shift @words;
	my $chembl = shift @words;
	my $annotation = shift @words;
	chomp($annotation);
	my $num_site = get_num_site($dbh, $accession);
	if ($num_site == 1){	#Protein with single binding domain
		my $start = get_start_position($dbh, $accession);
		my $end = get_end_position($dbh, $accession);
		my $subseq = get_substring($dbh, $accession, $start, $end);
		output_fasta($genename, $accession, $start, $end, $annotation, $subseq);
		get_pdb($genename, $accession);
		next;
	} elsif ($num_site > 1) {	#Protein with multiple binding domains
		my @starts = get_start_position($dbh, $accession);
		my @ends = get_end_position($dbh, $accession);
		my $start = min (@starts);
		my $end = max (@ends);
		my $subseq = get_substring($dbh, $accession, $start, $end);
		output_fasta($genename, $accession, $start, $end, $annotation, $subseq);
		get_pdb($genename, $accession);
		next;
	} else {	#0 binding domain or no information (NULL)
		print "Protein $accession: num_site is 0 or NULL ($num_site)\n";
		next;
	}
}
close $ZINC;

sub get_pdb {
	my ($genename, $accession) = @_;
	my @PDBs = get_PDB_by_UniProt($accession);
	foreach my $pdb (@PDBs){
		chomp($pdb);
		download_PDB($genename, $pdb);
	}
	return;
}

sub download_PDB {
        my ($genename, $pdb) = @_;
        chomp($pdb);
        my $url = "http://www.rcsb.org/pdb/files/".$pdb.".pdb";
	my $filepath = "../protein/$genename"."/";
	make_dir($filepath);
	my $pdb_file = $filepath . $pdb . ".pdb";
        getstore($url, $pdb_file);
	return;
}
sub get_PDB_by_UniProt {
        my $UniProt = shift @_;
        my $XML_query = qq(
<?xml version="1.0" encoding="UTF-8"?>
<orgPdbQuery>    
    <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
    <description>Simple query for a list of UniprotKB Accession IDs: $UniProt</description>   
    <accessionIdList>$UniProt</accessionIdList>
</orgPdbQuery>
        );
        # Create a request                                                                                  
        my $request = HTTP::Request->new( POST => 'http://www.rcsb.org/pdb/rest/search/');
        $request->content_type( 'application/x-www-form-urlencoded' );
        $request->content( $XML_query );

        # Post the XML query                                                                                
        my $response = $ua->request( $request );

        # Check to see if there is an error
        unless( $response->is_success ) {
            print "\n an error occurred: ", $response->status_line, "\n";
        }
        # Print response content in either case
        my $result = $response->content;
        my @PDBs = split(/\n/, $response->content);
        return @PDBs;
}

sub divide_string {
	my ($string, $cut_length) = @_;
	my @cut_strings;
	while ($string){
		my $cut_string = substr ($string, 0, $cut_length);
		push (@cut_strings, $cut_string);
		$string = substr($string, $cut_length);	#the left string
	}
	return @cut_strings;
}
sub output_fasta {
	my ($genename, $accession, $start, $end, $annotation, $sequence) = @_;
	my $filepath = "../protein/$genename"."/";
	make_dir($filepath);
	my $outfile = $filepath . $genename . ".fas";
	my $length = $end - $start + 1;
	open my $FASTA, '>', $outfile or die "Could not open FASTA file $outfile: $!\n";
	print $FASTA ">emb|$accession|length:$length($start:$end)|$annotation\n";
	my @divided = divide_string($sequence, 70);
	foreach my $seq (@divided){
		print $FASTA $seq, "\n";
	}
	close $FASTA;
	return;
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

sub get_start_position {
	my ($dbh, $accession) = @_;
	my $sth = $dbh->prepare("SELECT cd.start_position FROM component_domains cd 
				INNER JOIN component_sequences cs ON cd.component_id = cs.component_id 
				INNER JOIN site_components sc ON cs.component_id = sc.component_id AND cd.domain_id = sc.domain_id 
				WHERE cs.accession = '$accession'");
        $sth->execute() or die $DBI::errstr;
	my $num_row = $sth->rows;
	if ($num_row == 1){
		my $start = $sth->fetchrow_array();
		$sth->finish();
		return $start;	#single site
	} else {
		my @start_positions;
		while (my $start = $sth->fetchrow_array()){
			chomp($start);
			push (@start_positions, $start);
		}
		$sth->finish();
		return @start_positions;
	}
}
sub get_end_position {
	my ($dbh, $accession) = @_;
	my $sth = $dbh->prepare("SELECT cd.end_position FROM component_domains cd 
				INNER JOIN component_sequences cs ON cd.component_id = cs.component_id 
				INNER JOIN site_components sc ON cs.component_id = sc.component_id AND cd.domain_id = sc.domain_id 
				WHERE cs.accession = '$accession'");
        $sth->execute() or die $DBI::errstr;
	my $num_row = $sth->rows;
	if ($num_row == 1){
		my $end = $sth->fetchrow_array();
		$sth->finish();
		return $end;	#single site
	} else {
		my @end_positions;
		while (my $end = $sth->fetchrow_array()){
			chomp($end);
			push (@end_positions, $end);
		}
		$sth->finish();
		return @end_positions;
	}
}
sub get_substring {
	my ($dbh, $accession, $start, $end) = @_;
	my $length = $end - $start + 1;
	my $sth = $dbh->prepare("SELECT SUBSTRING(sequence, $start, $length) FROM component_sequences 
				WHERE accession = '$accession'");
	$sth->execute() or die $DBI::errstr;
	my $num_row = $sth->rows;
	my $sequence = $sth->fetchrow_array();
	$sth->finish();
	return $sequence;
}

sub get_num_site {
	my ($dbh, $accession) = @_;
	my $sth = $dbh->prepare("SELECT COUNT(DISTINCT(cd.compd_id)) FROM component_domains cd 
				INNER JOIN component_sequences cs ON cd.component_id = cs.component_id 
				INNER JOIN site_components sc ON (sc.component_id = cs.component_id) AND (cd.domain_id = sc.domain_id) 
				WHERE cs.accession = '$accession'"); 
	$sth->execute() or die $DBI::errstr;
	my $num_site = $sth->fetchrow_array();
	$sth->finish();
	return $num_site;
}

sub get_accession_multidomain {
	my $dbh = shift @_;
	my $sth = $dbh->prepare("SELECT cs.accession FROM component_sequences cs 
				INNER JOIN component_domains cd ON cs.component_id = cd.component_id 
				INNER JOIN site_components sc ON cs.component_id = sc.component_id AND cd.domain_id = sc.domain_id 
				GROUP BY cs.accession 
				HAVING (COUNT(DISTINCT(cd.compd_id)) > 1 AND COUNT(DISTINCT(cd.domain_id)) > 1)");
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
