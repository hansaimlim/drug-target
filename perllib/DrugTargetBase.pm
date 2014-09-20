#!/usr/bin/perl
package DrugTargetBase;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(chomp_array unique make_dir dirname_add_slash reverse_simple_hash one_column_file_switch rm_special_char_in_drugname);

sub rm_special_char_in_drugname
{
	#specialized for cmap drugnames, some of which contain slashes, dashes, parenthesis or signs
        my $drug = shift @_;
        $drug =~ s/\///g;       #remove slash
        $drug =~ s/(\-)(\d|\w)/_$2/;    #a dash to an underscore (but not the stereochemical minus sign)
        $drug =~ tr/()/__/;     #parenthesis to underscores
        $drug =~ tr/+-/pm/;     #stereochemical signs to letters
        $drug =~ tr/\\//;       #remove backslash -- just in case
        return $drug;
}
sub one_column_file_switch
{
        #input file must be one element per line
        my ($file) = shift @_;
        my %hash;
        open my $FH, '<', $file or die "Could not open file $file: $!\n";
        while (my $line = <$FH>){
                chomp($line);
                $hash{$line} = 1;
        }
        close $FH;
        return \%hash;
}

sub chomp_array
{
        #input an array then chomp each element
        my $r = shift @_;
        my @array = @$r;
        my @chomparray;
        foreach my $i ( @array ){
                chomp($i);
                push @chomparray, $i;
        }
        return @chomparray;
}

sub unique
{
        #input an array then output an array with duplicated elements removed
        my $r = shift @_;
        my @array = @$r;
        my %seen;
        my @unique = grep { ! $seen{$_}++ } @array;
        return @unique;
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
sub reverse_simple_hash
{
        #input: a hash reference
        #output: a reference to the reverse hash
        #warning: the original hash must be 1 dimensional--$hash{key} = scalar value
        #warning: the key value pairs must be 1 to 1 relationship; possible loss of information otherwise
        my $hash_ref = shift @_;
        my %hash = %$hash_ref;
        my %reverse_hash;
        foreach my $key (keys %hash){
                my $val = $hash{$key};
                $reverse_hash{$val} = $key;
        }
        return \%reverse_hash;
}

1;
