#!/usr/bin/perl
package DrugTargetBase;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(chomp_array unique make_dir dirname_add_slash);

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

1;
