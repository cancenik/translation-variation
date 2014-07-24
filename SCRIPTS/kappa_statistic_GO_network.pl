#!/usr/bin/perl -w 

## CC 
## May 26 2014

###########
# This script takes a Gene Ontology file from FuncAssociate
# In addition a results file is specified. 
# For all the significant terms; 
# We calculate the kappa similarity between the terms; 
# We output two files: 
# 1- Kappa similarity network among the significant terms
# 2- An attribute table that reformats FuncAssociate Results; p-val LOD 
###########

use strict; 
use Getopt::Long; 
use constant GO_ID => 5; 
use constant PVAL => 4;
use constant LOD =>2;

my ( $go_annotations, $go_results) = ""; 
GetOptions ("go_annotations=s" => \$go_annotations, 
	    "go_results=s" => \$go_results) || die "Invalid arguments!"; 

my $cytoscape_significant_term_attributes = $go_results . ".attr"; 
my $kappa_network = $go_results . "_Kappa_Network.sif"; 

open (INPUT, $go_results ) || die "Cannot open go results: $!\n"; 
open (GO_DAG, $go_annotations ) || die "Cannot open GO DAG:$!\n"; 

open (RESULTS, ">$cytoscape_significant_term_attributes" ); 
open (KAPPA, ">$kappa_network" ) ; 

my %significant_terms; 
while (<INPUT> ) { 
    if (/^#/) { 
	next; 
    }
    elsif(/GO/) { 
	my @F = split; 
	$significant_terms{$F[GO_ID]} = (); 
## CHECK FORMAT
	print RESULTS "$F[GO_ID]\t$F[PVAL]\t$F[LOD]\n"; 
    }
    else { 
	next;
    }
}
close RESULTS; 
close INPUT; 

# DECIDE BEST DATA STRUCTURE FOR KAPPA STAT CALCULATION
# One solution is to keep a hash of arrays. 
# Then loop through the hash elements and calculate 
# Kappa with all others. 
my %go_dag; 
my @all; 
while ( <GO_DAG>) { 
    chomp; 
    if ( /^#/ ) { 
	next; 
    }
    else { 
	my @F = split; 
	if ( exists $significant_terms{$F[0]} ) {
	    foreach (@F) { 
		if (/ENSG/) { 
		    push (@{$go_dag{$F[0]}}, $_);
		    push (@all, $_); 
		}
	    }
	}
    }
}

# CALCULATE ALL UNIQUE ELEMENTS in @all
my %hash   = map { $_, 1 } @all;
my @uniq_all = keys %hash; 
print KAPPA "GO_TERM\t", join("\t", @uniq_all), "\n";
calculate_kappa(\%go_dag); 

# Write function to directly print to file
# The easiest implementation is through set operations

sub calculate_kappa  { 
    my $dag_ref = shift; 
    while ( my ($key, $value) = each (%{$dag_ref} ) ) { 
	my %go_genes  = (); 
	print KAPPA "$key"; 
	foreach (@{$value} )  { 
	    $go_genes{$_} = (); 
	}
	foreach (@uniq_all ) { 
	    if ( exists $go_genes{$_} ) { 
		print KAPPA "\t1"
	    }
	    else { 
		print KAPPA "\t0"; 
	    }
	}
	print KAPPA "\n";
    }
}

