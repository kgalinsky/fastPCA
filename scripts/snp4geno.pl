#!/usr/bin/env perl

use strict;
use warnings;

open my $fh, '<', $ARGV[0] or die $!;
my @stat = stat($fh);
my $line = <$fh>;
close $fh;

my $snps = $stat[7] / length($line);

local $, = "\t";
local $\ = "\n";
for ( my $i = 0 ; $i < $snps ; $i++ ) {
    print sprintf( 'snp%08d', $i+1 ), 1, 0.0, $i * 1000, 'A', 'C';
}
