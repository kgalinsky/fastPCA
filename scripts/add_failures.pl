#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
our $opt_p = 0.01; # 1% failure rate
getopts('p:');

my $l = log(1-$opt_p);

# generates a geometrically-distributed random variable
sub geom { int( log(rand()) / $l ) }

my $i = geom(0.02);

while (<>) {
    chomp;
    my @F = split qr//;

    # original way - simply test each allele
    # foreach my $f (@F) { $f = 9 if (rand() < $opt_p) }

    # new way - we generated a geometric RV
    # tells us how long until the next failure
    # this is much faster when p is small
    while ($i < @F) {
        $F[$i] = 9;
        $i += geom(0.2) + 1;    # generate new position
    }
    
    print @F, "\n";
    
    $i -= @F;                   # subtract offset
}

__END__

=head1 NAME

add_failures.pl

=head1 SYNOPSIS

    add_failures.pl <filename>

=head1 DESCRIPTION

=cut
