#!/usr/bin/env perl

use strict;
use warnings;

use Math::Random qw/
  random_uniform
  random_beta
  random_binomial
  /;

my $unif_low  = 0.05;
my $unif_high = 0.95;

my $m = $ARGV[0];

my $struct;
eval "\$struct = $ARGV[1];";

sub random_pop {
    my ( $p, $struct ) = @_;

    if ( ref $struct ) {
        my $FST = $struct->[0];
        my $F   = ( 1 - $FST ) / $FST;

        my $aa = $p * $F;
        my $bb = ( 1 - $p ) * $F;

        return (
            [
                $p,
                map { random_pop( random_beta( 1, $aa, $bb ), $struct->[$_] ) }
                  ( 1 .. $#$struct )
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
