#!/usr/bin/env perl

use strict;
use warnings;

use Math::Random qw/
  random_uniform
  random_beta
  random_binomial
  /;

my $unif_low  = 0.1;
my $unif_high = 0.9;

my $m = $ARGV[0];
$m =~ s/k/000/;

my $struct;
$ARGV[1] =~ s/k/000/g;
eval "\$struct = $ARGV[1];";
die $@ if ($@);

sub calc_F {
    my ($struct) = @_;
    if ( ref $struct ) {
        my $FST = $struct->[0];
        my $F   = ( 1 - $FST ) / $FST;
        $struct->[0] = $F;

        foreach my $substruct (@$struct[1..$#$struct]) {
	    calc_F($substruct);
        }
    }
}

calc_F($struct);

sub random_pop {
    my ( $p, $struct ) = @_;

    if ( ref $struct ) {
        my $F  = $struct->[0];
        my $aa = $p * $F;
        my $bb = ( 1 - $p ) * $F;

        return (
            [
                $p,
                map {
                    my $p = 0;
                    while ($p < 0.01 || $p > 0.99) { $p = random_beta( 1, $aa, $bb ) }
                    random_pop( $p, $struct->[$_] )
                } ( 1 .. $#$struct )
            ]
        );
    }
    else {
        return ( [ $p, [ random_binomial( $struct, 2, $p ) ] ] );
    }
}

sub p_str {
    my ($out) = @_;
    if ( @$out == 2 ) {
        return ( $out->[0] );
    }
    else {
        return ($out->[0] . '('
              . join( ',', map { p_str($_) } @$out[ 1 .. $#$out ] )
              . ')' );
    }
}

sub pop_str {
    my ($out) = @_;
    if ( @$out == 2 ) {
        return ( join( '', @{ $out->[1] } ) );
    }
    else {
        return ( join( '', map { pop_str($_) } @$out[ 1 .. $#$out ] ) );
    }
}

local $\ = "\n";
while ( $m-- ) {
    my $p = random_uniform( 1, $unif_low, $unif_high );
    my $out = random_pop( $p, $struct );

    print STDERR p_str($out);
    print STDOUT pop_str($out);
}
