#!/usr/bin/perl -w

use strict; 
use lib '/home/ccenik/SCRIPTS/';
use Simple_Utils;
use Getopt::Long;
use Data::Dumper;
$Data::Dumper::Indent = 1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Pair = ":";
$Data::Dumper::Sortkeys = 1;

## CC Nov 8
## Modified parameter reading to work with qsub passed environment variables
## Added check for MAPQ >=2 
## Currently all reads are counted not just CDS ones

## CC May 28 
## Added Pileup to stop codon as well as the start codon

######################
## CC May 18 2013
## Discovered bug due to Appris annotation. 
## The name has UTR3, UTR5 and CDS information. 
## If UTR3 is indeed before the CDS; the transcript is on the - strand
## However, UTR3 doesn't mean 3'UTR for transcripts on the - strand
## For - strand transcript UTR3 is indeed 5'UTR
## FIX ANY ISSUES RELATED TO THIS
######################

######################
## CC Feb 6 ## 
## This script takes sam alignment format to nonredundant transcriptome  as input and outputs
## 1- Read Counts per transcript (only + strand alignments) Report only >10
## 2- Read Length distribution of aligning reads
## 3- Distance metrics for each read length class
## 4- Fraction of reads mapping to 5'UTR, CDS, or 3'UTR
######################

## Apr 19 2013 
## CC Update
## Redesign commandline to allow more flexible passing of parameters such as use species
# Specifically introduce --species flag to enable the using of the first occurence of identitcal sequences
# This should be used only whent the sam output is pre-processed using a alignment quality cutoff
######################

use constant ID => 2; 
use constant FLAG =>1; 
use constant SEQ => 9;  
use constant POS => 3;
# READ LENGTH BOUNDARIES WITH RESPECT TO CDS
use constant DOWN => 35;
use constant UP => -35;
use constant MAPQ => 4;  

my $use_species =  1; 
my $last_nt = 'N'; 
my $input_file = $ENV{'file'};
#my $annotation_file = ''; 
#GetOptions ('file=s' =>\$input_file , 'species' => \$use_species, 'nucleotide=s' => \$last_nt, 'annotation=s' => \$annotation_file); 

unless ($input_file)  {
#    print scalar @ARGV;
    die "Usage samtranscriptaln_tocounts.pl --file SAM_FILE [--species | --nucleotide [ACTG] | --annotation_file ANNOT_FILE]\n";
}

# CHANGE DEFINITION OF MODULO DATA TO HASH BY THE LAST TWO NUCLEOTIDES
my (%transcript_counts, %seq_lengths, %dist_hist, %modulo_data, %species, %dist_to_stop); 
open (IN, $input_file) || die "Cannot open SAM file: $!";
while (<IN>) { 
    my @F = split; 
    if ($F[0] =~ /^@/) {
	next;
    }
# CHECK LAST NUCLEOTIDE FIRST
    unless ($last_nt eq "N" || substr ($F[SEQ], -1) eq "$last_nt" ) {
	next;
    }
# CHECK WHETHER WE SAW THE SPECIES BEFORE
    if ($use_species) { 
	if (exists $species{"$F[SEQ]"}) { 
	    next; 
	}
	else { 
	    $species{"$F[SEQ]"}++;
	}
    }

# TRANSCRIPT COUNTS
# If the read maps to rev-comp, seq is also rev-comped and flag is 0. 
# Flag is 16 only when the read is rev NOT rev-comp. 
    unless ( $F[FLAG] ==0 && $F[MAPQ] >= 2 ) {
	next;
    }
    $transcript_counts{$F[ID]} ++; 
# SEQ_LEN_DIST
    my $seq_len = length($F[SEQ]); 
    $seq_lengths{$seq_len}++; 

# DISTANCE TO CDS_START grouped by SEQ_LEN. 
### ID INFORMATION HAS THE UTR, CDS LENGHTS
### FOR SOME_IDS THERE IS NO UTR
### NEED TO USE THE STRAND INFO IN THE ALIGNMENT
### IF NO UTR CHECK FOR ATG

## CC May 18 2013-- NOTE ABOVE ABOUT UTR!
# Commented out previously necessary steps...
#    my $utr3_index = index($F[ID], "UTR3");
#    my $utr5_index = index($F[ID], "UTR5"); 
    my $cds_index = index($F[ID], "CDS");
    if ($cds_index == -1) {
# THERE ARE A FEW (3 in EMI) without CDS INFO
	next;
    }
#    if ($utr5_index == - 1) { 
#	if ($utr3_index== -1) { 
	# No UTR-  ~140 Genes -- CURRENTLY IGNORED
	# IMPLEMENT SOLUTION
#	}
#	elsif ($utr3_index > $cds_index) { 
	    # + STRAND
#	    find_pos(\@F, "+", $last_nt);
#	}
#	else { 
	    # - STRAND
#	    find_pos(\@F, "-",  $last_nt);
#	}
#    }
    else {
#	if ($utr5_index > $cds_index) {
	    # - STRAND
#	    find_pos(\@F, "-", $last_nt);
#	}
#	else { 
	    # + STRAND
#	    find_pos(\@F, "+", $last_nt);
#	}
	find_pos(\@F, "+", $last_nt);
	find_distance_to_stop(\@F, $last_nt); 
    }
} 

## CREATE OUTPUT FILES
my $count_file = $input_file ."_transcript_counts";
my $len_file = $input_file . "_lengths_of_aligned"; 
my $pileup_file = $input_file . "_readpileup_tocds";
my $modulo_file = $input_file . "_modulo_Counts";
my $stop_dist_file = $input_file. "_distance_to_stopcodon"; 
open my $CTF, '>', $count_file || die "Cannot open COUNT file: $!";
open (my $LENF, ">$len_file") || die "Cannot open LEN file: $!";
open (my $PLF, ">$pileup_file") || die "Cannot open PILEUP file: $!";
open (my $MDF , ">$modulo_file") || die "Cannot open MODULO file: $!"; 
open (my $DTF , ">$stop_dist_file") || die "Cannot open STOP file: $!"; 
Simple_Utils::print_hash (\%transcript_counts, 10, $CTF); 
Simple_Utils::print_hash (\%seq_lengths, 0, $LENF);
Simple_Utils::print_hash (\%dist_hist, 0, $PLF); 
Simple_Utils::print_hash(\%modulo_data, 0, $MDF); 
Simple_Utils::print_hash(\%dist_to_stop, 0, $DTF); 

#print Dumper \%dist_hist;
#while( my ($k, $v) = each %dist_hist ) {
#    for my $k1 ( sort keys %$v ) {
#        print "$k\t$k1\t$v->{$k1}\n";
#    }
#}


# Only_Use -35  to 35 FOR CDS
# Takes each sam line and strand and sequence last nucleotide as input
# In addition to detailed output for DOWN to UP; report a modular vector grouped by length and last nucleotide
# Ex. 24 A 1090 2500 4500
# This information is important for testing 3nt periodicity
# We will limit the counting to -20 TO LEN(CDS)
# Modulo data is a hash accessible directly from the function
# We might consider using only species instead of reads for counting modulo
sub find_pos {
    my $line = shift;
    my $strand = shift;
    my $last_nuc = shift; 
# IF LAST NUCLEOTIDE IS NOT CONSIDERED PROCESS EVERYTHING
    my $seq_len = length($$line[SEQ]); 
    my ($dist, $cds_len); 
    if ($strand eq "+") {
	if (${$line}[ID] =~ /CDS:(\d+)-(\d+)/) {
	    $dist = ${$line}[POS] - $1;
	    $cds_len = $2 - $1; 
	    if ($dist < DOWN && $dist > UP ) {
		$dist_hist{"$seq_len"}{"$dist"}++; 
	    }
	}
	else {
	    print "Shouldn't happen\n"; 
	}
    }
    elsif ($strand eq "-") {
	if (${$line}[ID] =~ /CDS:(\d+)-(\d+)/) {
	    $dist = $2 - ${$line}[POS] ;
	    $cds_len = $2 - $1; 
 	    if ($dist < DOWN && $dist > UP ) {
		$dist_hist{"$seq_len"}{"$dist"}++; 
	    }
	}
    }
    else {
	print "WRONG SYNTAX\n";
    }
    if ($dist > UP && $dist < $cds_len) {
	my $mod = $dist % 3;
	my $last_two_nt = substr($$line[SEQ], -2); 
	$modulo_data{"$seq_len"}{"$last_two_nt"}{"$mod"}++;
    }
}

# This function finds the distance of the read to last codon
# It will update hash if the distance is shorter than certain limits
# For now, the last nucleotide is not considered
# The read is always + strand
sub find_distance_to_stop {
    my $line = shift;
    my $last_nuc = shift;
# IF LAST NUCLEOTIDE IS NOT CONSIDERED PROCESS EVERYTHING                      
 #   my $seq_len = length($$line[SEQ]);
    my ($dist);
    if (${$line}[ID] =~ /CDS:(\d+)-(\d+)/) {
	$dist = $2 - ${$line}[POS];
	if ($dist < DOWN && $dist > UP ) {
	    $dist_to_stop{"$dist"}++;
#	    print "$dist_to_stop{$dist}\n";

	}
    }
    else {
	print "Shouldn't happen\n";
    }
}

