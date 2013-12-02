#!/usr/bin/perl -w

## CC 
## Jun 26 2013

##############
# This script performs simple QC on Damek's omni array data
##############

use strict; 
use constant GC_SCORE => 4; 
use constant ID => 0; 
use constant GT_SCORE => 11; 
use constant LOG_R => 24; 


my $omni_data = "/srv/gs1/projects/snyder/damek/Array1S.snptable.txt";

open (IN, $omni_data) || die "Cannot open input\n"; 

my $low_quality_snp_number = 0; 
my $throw_away_low_quality_in_another_sample = 0 ; 
my $inter_sample_clustering_score_low = 0; 
my $potential_CNV = 0; 
my %low_quality_snp_ids; 
while (<IN>) { 
    my @F = split; 
    # Filter low GC scores
    if ($F[GC_SCORE] < 0.25) { 
	$low_quality_snp_number ++; 
	$low_quality_snp_ids{$F[ID]}++;
	next;
    }
    elsif (exists $low_quality_snp_ids{$F[ID]}) {
	$throw_away_low_quality_in_another_sample ++;  
	next; 
    }
    # Filter low GC scores
    elsif ($F[GT_SCORE] < 0.7 ) { 
	$inter_sample_clustering_score_low ++; 
	next; 
    }
    # Filter potential CNV calls
    elsif (abs($F[LOG_R] ) > 0.5 ) { 
	$potential_CNV++; 
	next; 
    } 
    else { 
	my $line = join ("\t", @F); 
	print "$line\n";
    }
    
} 

print STDERR "Low_Quality_Genotypes_Thrown_Away: $low_quality_snp_number\n"; 
print STDERR "Low_Quality_Genotype_in_Other_Sample: $throw_away_low_quality_in_another_sample\n"; 
print STDERR "Inter_Sample_Clustering_Score_Low: $inter_sample_clustering_score_low \n"; 
print STDERR "Potential_CNV: $potential_CNV\n"
