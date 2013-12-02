#!/usr/bin/perl -w 

## CC
## May 24 2013

## This script converts variant format to BED format
## The intended use is with zcat pipe

# Variant format
 # 10:60523TNA18526=T|TNA18951=T|TNA12890=T|TNA10847=T|TNA18505=T|T
 #  10:60969CNA18526=C|CNA18951=C|CNA12890=C|CNA10847=C|CNA18505=C|C
 #  10:61005ANA18526=A|ANA18951=A|ANA12890=A|ANA10847=A|ANA18505=A|A
 #  10:61020GNA18526=G|GNA18951=G|GNA12890=G|GNA10847=G|GNA18505=G|G

use strict; 

while (<>) { 
    my @F = split; 
    my $first = shift (@F);
    my @G = split (/:/, $first);
    my $start = $G[1]-1;
    my $name = join ("_", @F);
    print "chr$G[0]\t$start\t$G[1]\t$name\n"; 
}
