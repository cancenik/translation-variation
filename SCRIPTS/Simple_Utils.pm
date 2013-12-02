#!/usr/bin/perl -w

## CC ## 
## JAN 2013 ## 

##################
# Simple helper functions that are utilized by many scripts
##################

package Simple_Utils;
use strict;

# Parameters: Hash_Ref; Lowest Value to print; lexical filehandle
sub print_hash {
    my $hash_ref = shift;
    my $lowest_value = shift; 
    my $filehandle = shift;
    for my $k (sort keys %$hash_ref) {
	my $v = $hash_ref->{$k};
	if (ref $v ne ref {})  {
	    if (ref $v eq ref []) {
		print $filehandle "\t$k\t@$v\n";
	       
	    }
	    else  {
		print $filehandle "\t$k\t$v\n";
	    }
	}
	else {
	    print $filehandle "$k\n"; 
	    print_hash ($v, $lowest_value, $filehandle); 
	}
    }
}


# Given two input files, the script outputs the fasta sequences from the second file
# only if there is a matching in the first file.
# Currentl selects ENST ids. Can be generalized. 
sub select_fasta {
    my $id_file = shift;
    my $fasta_file = shift; 
    my $select_id = 2;
    open IN, $id_file or die $!;
    my %match_ids;
    while (<IN>) {
	my @F = split;
	$match_ids{$F[$select_id]} = 1;
    }
    close IN;

    open INI, $fasta_file or die $!;
    my $i = 0;
    while (<INI>) {
	if (/>(ENST\d+)/) {
	    if (exists $match_ids{$1}) {
		$i = 1;
	    }
	    else {
		$i = 0;
	    }
	}
	if ($i) {
	    print;
	}
	else {
	    next;
	}
    }
}

# Given two input files, the script outputs the the contents of the second file
# only if there is a matching in the first file.
# Selects RefSeq Ids for now can be generalized if there is a need to 
sub select_subset {
    my $filter_file = shift; 
    my $selection_file = shift; 
    my $select_id = 0; 
    my $use_id = 0 ; 
    open IN, $filter_file or die $!;
    my %match_ids;
    while (<IN>) {
	my @F = split;
	$match_ids{$F[$select_id]} = 1;
    }
    close IN;

    open INI, $selection_file or die $!;
    while (<INI>) {
	my @F = split;
	if ($F[$use_id] =~ /([NX]M_\d+)/) {
	    if (exists $match_ids{$1}) {
		print "@F\n";
	    }
	    else {
		next;
	    }
	}
    }   
}
# Given a string, returns the reverse strand DNA
sub rev_dna {
    my $dna = shift; 
    $dna =~ tr/ACGTacgt/TGCAtgca/;
    return $dna; 
}
1; 
