#!/usr/bin/perl
use strict;
use warnings;
use XML::SAX::ParserFactory;
use XML::SAX::Base;
use Data::Dumper;

#------------------------------------------------------------------------------------------------
# Purpose	: To extract genename and UniProtKB genesymbol match map
#		  from drugbank.xml
# Input		: drugbank.xml
# Author	: Hansaim Lim
#-----------------------------------------------------------------------------------------------
die "Usage: $0 <drugbank.xml>\n" unless @ARGV == 1;
my $xmlfile = shift @ARGV;

my %switches = ();	#define as they are seen
my $is_uniprotaccession = 0;	#1 if UniProt Accession is open
my $is_human = 0;		#1 if the target is for human

my $output2 = "genename_UniProtsymbol.tsv";
open my $OUTPUT2, '>', $output2 or die "Could not open file $output2: $!\n";
print $OUTPUT2 "gene-name\tUniProt-symbol\n";
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
			chomp($char);
			SWITCH: {
				if ($switches{"organism"}){
					return unless $switches{"target"};	#should be within a target tag
					if ($char =~ m/^human$/i){ $is_human = 1; }
				}
				if ($switches{"gene-name"}){
					return unless $switches{"target"};
					return unless $is_human;
					#-----------$OUTPUT-----------------
					print $OUTPUT2 $char, "\t";	#for gene-name genesymbol map
					#-----------$OUTPUT----------------
				}
				if ($char =~ m/UniProt\s*Accession/i) {$is_uniprotaccession = 1;}
				if ($switches{"identifier"}) {
					return unless $switches{"targets"};
					return unless $is_human;
					return unless $is_uniprotaccession;	#skip other identifiers
					$char =~ s/^(.*)_HUMAN/$1/ig;	#remove suffix _HUMAN
				#-----------$OUTPUT-----------------
				print $OUTPUT2 $char, "\n";	#for gene-name genesymbol map
				#-----------$OUTPUT----------------
				}
			}
		},
		end_element => sub {
			my $end = shift->{Name};
			$switches{$end} = 0;
			$is_human = 0 if $end =~ m/^target$/;	#should check if human for each target
			$is_uniprotaccession = 0 if $end =~ m/^identifier$/;
		}
	}
);
$parser->parse_uri($xmlfile);
close $OUTPUT2;
exit;

