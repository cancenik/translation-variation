#!/usr/bin/perl -w

use strict; 

## CC 
## Takes one input file and fixes the ID format from "GAPDH-001" to GAPDH

my $file = $ARGV[0]; 
open (IN, $file) || die "Cannot open file\n"; 
while (<IN>) { 
    my @F = split(/-/); 
    print substr ($F[0],1), "\n"; 
}
