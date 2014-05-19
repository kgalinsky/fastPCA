#!/usr/bin/env perl

use strict;
use warnings;

open my $fh, '<', $ARGV[0] or die $!;
my $line = <$fh>;
close $fh;

my $inds = length($line) - 1;

local $, = "\t";
local $\ = "\n";

my @gender = qw/ M F /;
my @label  = qw/ Case Control /;
for ( my $i = 0 ; $i < $inds ; $i++ ) {
    print sprintf( 'sample%08d', $i+1 ), $gender[rand(2)], $label[rand(2)];
}
