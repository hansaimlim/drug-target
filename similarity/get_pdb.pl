#!/usr/bin/perl
use warnings;
use strict;
use LWP::Simple qw( $ua getstore );

my $uni = 'P41595';
my @PDBs = get_PDB_by_UniProt($uni);
foreach my $pdb (@PDBs){
	chomp($pdb);
	download_PDB($pdb);
}

sub download_PDB {
	my $pdb = shift @_;
	chomp($pdb);
	my $url = "http://www.rcsb.org/pdb/files/".$pdb.".pdb";
	my $pdb_file = "./static/pdb/$pdb".".pdb";
        getstore($url, $pdb_file);
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

	# you can configure a proxy...                                                                          
	#$ua->proxy( http => 'http://yourproxy:8080' );

	# Create a request                                                                                  

	my $request = HTTP::Request->new( POST => 'http://www.rcsb.org/pdb/rest/search/');

	$request->content_type( 'application/x-www-form-urlencoded' );

	$request->content( $XML_query );

	# Post the XML query                                                                                
	print "\n querying PDB...";
	print "\n";
	my $response = $ua->request( $request );

	# Check to see if there is an error
	unless( $response->is_success ) {
	    print "\n an error occurred: ", $response->status_line, "\n";

	}

	# Print response content in either case
#	print $response->content;
	my $result = $response->content;
	my @PDBs = split(/\n/, $response->content);
	return @PDBs;
}
exit;
