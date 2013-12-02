#!/usr/bin/perl -w 

## CC 
## Jul 8 2013
## Updated on Nov 8 2013
################
# This script generates the CDS count table for LCL_RIBOSEQ DATASETS
## We will also output coverage and length of the region
# This script will also aggregate species data from appris_aln.sam_transcript_counts
################

use strict; 
use File::Basename;
#use IO::Zlib;
use constant ID => 4; 
use constant COUNT => 6; 
use constant COVERED_BASES => 7; 
use constant FEATURE_LEN => 8; 
use constant REGION => 3; 
use Data::Dumper; 
$Data::Dumper::Terse = 1;
$Data::Dumper::Indent = 0;
require Simple_Utils;

my %counts;
my @files = </srv/gs1/projects/snyder/ccenik/LCL_RIBOSEQ/*/Transcript_Counts_*>;
#my @files = </srv/gs1/projects/snyder/ccenik/LCL_RNASEQ/RIBOZERO_RNA/*/Transcript_Counts*>;
#my @files = </srv/gs1/projects/snyder/ccenik/LCL_RIBOSEQ/*/Species_Counts_UTRs_CDS_CoverageBed.tsv>;
my $counts_file = "Transcript_Counts_All_Libraries.tsv"; 
#my $counts_file = "Species_Counts_All_Libraries.tsv"; 
unless (-e $counts_file) {
open (OUT, ">$counts_file") || die "Cannot open output file: $!\n"; 
print OUT "ID\tREGION\t"; 
foreach my $file (@files) { 
    my $dirname = dirname($file);
    my $basename = basename($dirname);
    print OUT "$basename"."_Counts\t$basename" ."_CoveredBases\t$basename". "_FeatureLen\t"; 
#    open (IN, "gzip -dc $file |"); 
    open (IN, "$file") || die "Cannot open file: $!\n"; 
    while (<IN>) { 
	my @F = split; 
	my @ID = split (/\|/, $F[0]); 
	push (@{$counts{$ID[ID]}->{$F[REGION]}}, ($F[COUNT], $F[COVERED_BASES], $F[FEATURE_LEN])); 
    }
}
print OUT "\n"; 
#print OUT Dumper(\%counts); 

my $out_fh = *OUT;
Simple_Utils::print_hash(\%counts, 0 , $out_fh);
close OUT; 

#REFORMAT OUTPUT FOR R_INPUT
open (REFORMAT, $counts_file) || die "Cannot open file for reformatting\n"; 
my $header = <REFORMAT>; 
my $reformatted_file = "Reformatted_Transcript_Counts_All_Libraries.tsv";
#my $reformatted_file = "Reformatted_Species_Counts_All_Libraries.tsv";
open (FINALOUT, ">$reformatted_file") || die "Cannot open final output: $!\n";
print FINALOUT $header;
my $tmp_id; 
while (<REFORMAT>) {
    chomp;
    my @G = split (/\t/); 
    if ($G[0] eq "") {
	$G[2] =~ s/\s+/\t/g;
	print FINALOUT "$tmp_id\t$G[1]\t$G[2]\n";
    }
    else {
	$tmp_id = $G[0]; 
    }
}
close REFORMAT;
close FINALOUT; 
}

my @species_files = </srv/gs1/projects/snyder/ccenik/LCL_RIBOSEQ/*/appris_aln.sam_transcript_counts>;
my $species_merge = "Appris_Species_Counts"; 
unless (-e $species_merge) {
open (MERGE, ">$species_merge") || die "Cannot open $species_merge"; 
my %transcript_species;
print MERGE "ID\t";  
foreach my $file (@species_files) {
    my $dirname = dirname($file);
    my $basename = basename($dirname);
    print MERGE "$basename\t"; 
    open (SPECIES, $file) || die "Cannop open $file: $!\n"; 
    while (<SPECIES>) { 
	my @G = split;
	my @ID = split (/\|/, $G[0]);
	push (@{$transcript_species{$ID[ID]}}, ($basename, $G[1]));	
    }
}
print MERGE "\n";
#print OUT Dumper(\%counts);                                                                                              
my $out_fh = *MERGE;
Simple_Utils::print_hash(\%transcript_species, 0 , $out_fh);
close MERGE;

my $reformatted_merge_species = "Reformatted_Appris_Species_Counts"; 
open (REFORMAT_SPECIES, ">$reformatted_merge_species") || die "Cannot open :$!\n"; 
open (TOMERGE, "Appris_Species_Counts") ; 
my $library_list = <TOMERGE>; 
print REFORMAT_SPECIES $library_list;
my @library_list = split (/\s+/, $library_list); 
shift(@library_list);
#print "$library_list[0]\n$library_list[1]\n$library_list[2]\n";
while (<TOMERGE>) {
    chomp;
    my @G = split (/\t/, $_);
    if ($G[0] eq "") {
	my @W = split(/\s+/, $G[2]);
#    my $i = 0; 
#    foreach (@W) { 
#	print "$i\n$_\n"; 
#	$i++;
#    }
	print REFORMAT_SPECIES "$G[1]\t"; 
	my %hash = @W;	
	foreach my $lib (@library_list) {
	    if (exists $hash{$lib}) { 
		print REFORMAT_SPECIES "$hash{$lib}\t";
	    }
	    else { 
		print REFORMAT_SPECIES "0\t";
	    }
	} 
	print REFORMAT_SPECIES "\n";
    }
    else {
	print "Error found @G\n";
    }
}
close TOMERGE;
close REFORMAT_SPECIES; 
}
