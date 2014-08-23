#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;
use List::MoreUtils qw/ uniq /;
#--------------------------------------------------------------------------------------------------
# Purpose	: To convert PubChem CID to InChI keys && ENSP to GeneSymbols from STITCH
# Input		: STITCH chemical-protein-score file
# Output	: CID-InChI Key; two-column file (tab separated)
# Author	: Hansaim Lim
# Date		: 22 Aug, 2014
#--------------------------------------------------------------------------------------------------

die "Usage: $0 <STITCH file> <UniProt map file> <orf map file> <minimum score>\n" unless @ARGV == 4;
my $stitch = shift @ARGV;
my $uniprot_idmap = shift @ARGV;
my $orfmap = shift @ARGV;
my $min_score = shift @ARGV;
my $outfile = $stitch;
$outfile =~ s/^(.+)(\..+)$/$1InChIKey_GS_min$min_score$2/;

open my $UNIPROT, '<', $uniprot_idmap or die "Could not open file $uniprot_idmap: $!\n";
my %ensp_uniprots = ();
my %ensp_symbols = ();
while(<$UNIPROT>){
        next unless $_ =~ m/ENSP\d+/;
        my @ids = split(/\t/, $_);
        my $u = shift @ids;     #UniProtKB
        my $gs = shift @ids;    #gene symbol (with _HUMAN suffix)
        $gs =~ s/(.*)_HUMAN/$1/;#remove suffix

        foreach my $id ( @ids ){
                my @ensps = split(/; /, $id) if $id =~ m/ENSP\d+/;
                foreach my $ensp ( @ensps ){
                        $ensp_symbols{$ensp} = $gs;
                        $ensp_uniprots{$ensp} = $u;
                }
        }
}
close $UNIPROT;


open my $ORFMAP, '<', $orfmap or die "Could not open file $orfmap: $!\n";
my %ensp_orfs = ();
my $orfmap_line = 1;    #to skip the first line
ORFMAP: while(<$ORFMAP>){
        if ($orfmap_line == 1){
                $orfmap_line++;
                next ORFMAP;
        }
        my @ids = split(/\t/, $_);
        my $orf = shift @ids;
        my @ensps = shift @ids;
        my @uniprots = shift @ids;

        if ($_ =~ m/ENSP\d+/){ #if ENSP id found
                foreach my $ensp (@ensps){
                        $ensp_orfs{$ensp} = $orf;
                }
        } else {        #try to match using UniProtKB when the orfmap line does not contain ENSP information
                my @ENSPs;
                foreach my $uniprot (@uniprots){
                        foreach my $ENSP (keys %ensp_uniprots){
                                push @ENSPs, $ENSP if $ensp_uniprots{$ENSP} eq $uniprot;        #push the ensp ids for given uniprot
                        }
                }
                foreach my $ensp_additional (@ENSPs){
                        $ensp_orfs{$ensp_additional} = $orf;    #ensp-orf pairs for orfs without ensp listed on orfmap file
                }
        }
}
close $ORFMAP;


my $line = 1;
my $lineskip = 0;	#number of skipped lines
my $ikeyskip = 0;	#number of inchikeys not found
my $gsskip = 0;		#number of genesymbols not matched
open my $STITCH, '<', $stitch or die "Could not open file $stitch: $!";
open my $OUTPUT_STITCH, '>', $outfile or die "Could not open file $outfile: $!\n";
STITCH: while (<$STITCH>) {
	my @words = split(/\t/, $_);
	if ( $line == 1) {
		$line++;
		next STITCH;	#skip the first line
	}
	my $chemical = shift @words;
	my $ensp = shift @words;
	$ensp =~ s/9606\.(ENSP\d{11})/$1/;	#remove NCBI taxonomy ID (9606 for human)

	my $score = shift @words;	#no need to convert this score
	chomp($score);
	next STITCH unless $score ge $min_score;

	my $cid = $chemical;	#the chemical contains CID(1|0) and pubchem CID 
	$cid =~ s/CID\d(\d+)$/$1/;	#only pubchem cid number

	my $inchikey = pubchem_inchikey_by_cid($cid);
	my $genesymbol;
	$genesymbol = $ensp_symbols{$ensp} or $genesymbol = $ensp_orfs{$ensp};
	if (!$inchikey or !$genesymbol){
		$ikeyskip++ unless $inchikey;
		$gsskip++ unless $genesymbol;
		$lineskip++;
		next STITCH;
	}
	print $OUTPUT_STITCH $inchikey, "\t", $genesymbol, "\t", $score, "\n";
}
close $OUTPUT_STITCH;
close $STITCH;

#-----------------------------------SUBROUTINES----------------------------------------------------------------
sub pubchem_inchikey_by_cid {
	my $cid = shift @_;
        my $url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$cid/property/InChIKey/TXT";
        my $inchikey = get $url;
        chomp($inchikey);
        return $inchikey;
}
exit;
