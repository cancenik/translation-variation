#!/usr/bin/perl -w 

## CC ## 
## Nov 20 2012 ##

################
# This script adds exon numbers to transcript ids
# As a defaults we ignore strand. 
# We add E001 to first instance and increase to E002, etc until next id
################

use strict;
use constant ID => 3;  

open IN, $ARGV[0] or die $!;
my $exon_count = 1;
my $id = ""; 
my $tmp; 
while (<IN>) {
    my @F  = split;
    if ($F[ID] =~ /(NM_\d+)/) {
	$tmp = $1; 
    }
    else {
	print "Unexpected line: $F[ID]\n";
    }
#    $F[ID] =~ s/"//g;
#    if ($F[ID] eq $id ) {
    if ($tmp eq $id ) {
	$exon_count++; 
	if ($exon_count <10) {
#	    $F[ID] =~ s/;/:E00$exon_count/;
	    $tmp .= ":E00$exon_count";
	}
	else {
#	    $F[ID] =~ s/;/:E0$exon_count/;
	    $tmp .= ":E0$exon_count";
	}
    }
    else { 
	$exon_count = 1; 
#	$id = $F[ID];
	$id = $tmp;
#	$F[ID] =~ s/;/:E00$exon_count/;
	$tmp .= ":E00$exon_count";
#	print "$F[ID]\n";
    }
    $F[ID] = $tmp;
    my $updated_line = join("\t", @F); 
    print "$updated_line\n"; 
}
