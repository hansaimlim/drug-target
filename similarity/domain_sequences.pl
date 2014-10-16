#!/usr/bin/perl
use strict;
use warnings;
use DBI;

die "Usage: $0 <username> <password>\n" unless @ARGV == 2;

my $driver = "MySQL"; 
my $database = "chembl_16";
my $dsn = "DBI:$driver:database=$database";
my $userid = shift @ARGV;
my $password = shift @ARGV;
chomp($password);

my $accession = 'P41595';
my $dbh = DBI->connect($dsn, $userid, $password ) or die $DBI::errstr;
my $sth = $dbh->prepare("SELECT * FROM component_sequences WHERE accession = '$accession'");
$sth->execute() or die $DBI::errstr;
print "Number of rows found :" + $sth->rows;
while (my @row = $sth->fetchrow_array()) {
   my ($first_name, $last_name ) = @row;
   print "First Name = $first_name, Last Name = $last_name\n";
}
$sth->finish();
exit;
