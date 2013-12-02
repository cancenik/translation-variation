#!/usr/bin/perl -w

## CC ## 
## Feb 20 2013 ## 

#################
# This script takes a sam file and a kmer length as input
# It outputs either -a all such kmers or -top top 50 kmers
# The counting is based on hashing of overlapping kmers so do not input kmer_len << Seq_Len 
#################

use strict; 
require Simple_Utils; 
use constant SEQ => 9;

unless (@ARGV == 3) {
    die "Usage: kmer_counter SAM_FILE KMER_LEN NUMBER_OF_COUNTS_TO_REPORT\n";
}

my $file = $ARGV[0];
my $kmer_len = $ARGV[1];
my $num_counts_to_report = $ARGV[2];
open (IN, $file) || die "Cannot open SAM file: $!";
my $outfile = "$file" . "_kmer_counts"; 
open my $out_fh, '>' , $outfile || die "Cannot open OUTPUT FILE: $!\n"; 

my %sub_strings; 
### MAPQ, FLAGS, etc are not considered
while (<IN>) { 
    my @F = split; 
    if ($F[0] eq '@SQ') {
	next;
    }
    else {
	for (my $key = 0; $key <= (length($F[SEQ]) - $kmer_len); $key++) {
	    my $str = substr($F[SEQ], $key, $kmer_len); 
	    $sub_strings{$str} ++; 
#	    print "$str\t$sub_strings{$str}\n";
	}
    }
}

Simple_Utils::print_hash(\%sub_strings, $num_counts_to_report, $out_fh);

