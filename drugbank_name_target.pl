#!/usr/bin/perl
use strict;
use warnings;
use XML::SAX::ParserFactory;
use XML::SAX::Base;
use Data::Dumper;

#------------------------------------------------------------------------------------------------
# Purpose	: To extract drugname-target(UniProtKB) pair (No Action Information)
#		  from drugbank.xml
# Input		: drugbank.xml
# Output	: tsv file with multiple columns. drugname-target1-target2-....
# Author	: Hansaim Lim
# Date		: 25 Jul, 2014; last modified 23 Jul 2014
#-----------------------------------------------------------------------------------------------
die "Usage: $0 <drugbank.xml>\n" unless @ARGV == 1;
my $xmlfile = shift @ARGV;

my %nta = ();	#name-target-action hash
my %switches = ();	#define as they are seen
my ($name, $target, $action);	#buffer variables
my $is_uniprotaccession = 0;	#1 if UniProt Accession is open
my $known_action = 0;	#at the time of writing this code, drugbank.xml contains three answers for known-action for targets; yes, no, unknown.

my $output = "db-name-targetsymbol.tsv";
open my $OUTPUT, '>', $output or die "Could not open file $output: $!\n";

my $factory = new XML::SAX::ParserFactory;
my $handler = new XML::SAX::Base;
my $parser = $factory->parser(
	Handler => $handler,
	Methods => {
		start_element => sub {
			my ($self, $dataca) = @_;
			my $start = $self->{Name};
			{no warnings 'uninitialized';
				$switches{$start} = 1 unless defined;
				$switches{$start} = 1;
			}
			$known_action = 0 if $start =~ m/^known-action$/;	#action unknown unless specified as "yes"
		},
		characters => sub {
		#list of data to skip
			return if $switches{"description"};
			return if $switches{"groups"};
			return if $switches{"general-references"};
			return if $switches{"synthesis-reference"};
			return if $switches{"indication"};
			return if $switches{"pharmacodynamics"};
			return if $switches{"mechanism-of-action"};
			return if $switches{"toxicity"};
			return if $switches{"metabolism"};
			return if $switches{"absorption"};
			return if $switches{"half-life"};
			return if $switches{"protein-binding"};
			return if $switches{"route-of-elimination"};
			return if $switches{"volume-of-distribution"};
			return if $switches{"clearance"};
			return if $switches{"classification"};
			return if $switches{"salts"};
			return if $switches{"synonyms"};
			return if $switches{"brands"};
			return if $switches{"mixtures"};
			return if $switches{"packagers"};
			return if $switches{"manufacturers"};
			return if $switches{"prices"};
			return if $switches{"categories"};
			return if $switches{"affected-organisms"};
			return if $switches{"dosages"};
			return if $switches{"atc-codes"};
			return if $switches{"ahfs-codes"};
			return if $switches{"patents"};
			return if $switches{"food-interactions"};
			return if $switches{"drug-interactions"};
			return if $switches{"sequences"};
			return if $switches{"experimental-properties"};
			return if $switches{"external-links"};
			return if $switches{"pathways"};
			return if $switches{"reactions"};
			return if $switches{"snp-effects"};
			return if $switches{"snp-adverse-drug-reactions"};
			return if $switches{"enzymes"};
			return if $switches{"carriers"};
			return if $switches{"transporters"};
			return if $switches{"calculated-properties"};	# SMILEs are under this tag
		#end list of data to skip
			my ($self, $dataca) = @_;
			my $char = $self->{Data};
			return if $char =~ m/^\s*$/;	#skip empty data
			SWITCH: {
				if ($switches{"name"}) { 
					return if $switches{"targets"}; #no need for target name
					$name = $char;
				#-------$OUTPUT----------------
				print $OUTPUT $name;
				#-------$OUTPUT---------------
				}
				if ($char =~ m/UniProt\s*Accession/i) { $is_uniprotaccession = 1; }
				if ($switches{"identifier"}) {
					return unless $switches{"targets"};
					return unless $is_uniprotaccession == 1;	#skip other identifiers
					$target = $char;
				#-----------$OUTPUT-----------------
				print $OUTPUT "\t", $target;
				#-----------$OUTPUT----------------
				}
			}
		},
		end_element => sub {
			my $end = shift->{Name};
			$switches{$end} = 0;
			$is_uniprotaccession = 0 if $end =~ m/^identifier$/;
		#--------------$OUTPUT-----------
		print $OUTPUT "\n" if $end =~ m/^targets$/;
		#--------------$OUTPUT-----------
		}
	}
);
$parser->parse_uri($xmlfile);
close $OUTPUT;
exit;

