#!/usr/bin/perl -w 

## CC 
## This takes results table output from limma and outputs input for radial sets
## Output Format: 
## IND_ID GENE1 GENE2 ... GENE340


use strict; 
use Data::Dumper;
my $results_file = $ARGV[0]; 
open (IN, $results_file) || die "Cannot open input file: $!\n"; 
# Process Header line
# Create Hash with individual IDs 
my $header = <IN>; 
my @G = split(/"\s+"/, $header); 
my (@up, @down);
foreach my $ind (@G) { 
    my @F = split(/-/, $ind); 
    $F[0] =~ s/"//; 
    chomp($F[0]); 
    if ($F[0] =~ /GM/) {
	push @up, "$F[0]";
	push @down, "$F[0]";
    }
}

while (<IN>) { 
    my @F = split; 
    my $gene_id = shift(@F); 
    my @gene_id = split(/-/, $gene_id); 
    $gene_id[0] =~ s/"//; 
    my $i =0; 
    foreach my $result (@F) { 
	if ($result == -1) { 
	    $down[$i] .= "\t$gene_id[0]"; 
	}
	elsif ($result == 1) { 
	    $up[$i] .= "\t$gene_id[0]";
	}
	else { 
	    #No differential
	}
	$i++; 
    }
}
#print Dumper(\@up); 

my $up_file = $ARGV[0] . "_UP"; 
open (UP, ">$up_file") || die "Cannot open UP output: $!\n"; 
foreach (@up ) { 
    print UP $_, "\n"; 
}

my $down_file = $ARGV[0] . "_DOWN"; 
open (DOWN, ">$down_file") || die "Cannot open DOWN output: $!\n"; 
foreach (@down ) { 
    print DOWN $_, "\n"; 
}
