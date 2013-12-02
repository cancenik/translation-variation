#!/usr/bin/perl -w 

## CC ## 
## MAY 18 2013 ## 

#######################
# This script takes a set of variants (bed format) intersected with APPRIS.gtf by intersectBED
# We need to find the relative position of the variant on the transcript
# We will also classify variants based on position
# We will only consider SNPs

## We should also keep the direct mapping to genome as an alternative strategy if needed
# OLD DEFINITION
# Variants were in gtf format
# It maps the position of the variant onto the transcript fasta sequence
# The main application is to create a variant of Appris transcript file with variants
# This modified sequence will be used to create enhanced reference for mapping

## Include output in 1ksnp format for use with ase 
# Example Line: 1 72515 C A 0.5 99 1 rs4030300
# Needs to have all fields
# chr number comes from appris_selected_transcripts_numbered 
#######################

# EXAMPLE USAGE: 
#map_variant_to_fasta_using_transcript_gtf.pl --variant variant_test --gtf /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris.gtf --fasta /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/appris_selected_transcripts.fa

use strict; 
use Getopt::Long;
use Data::Dumper;
use List::Util qw (sum); 
require Simple_Utils; 


my $gtf_file = ""; 
my $fasta_file = ""; 
my $variant_file = "";
my $num_file = "" ;
GetOptions ('variant=s' =>\$variant_file , 'gtf=s' => \$gtf_file, 'fasta=s' => \$fasta_file, 'chr_num=s' => \$num_file);

unless ($variant_file)  {
    die "Usage map_variant_to_fasta_using_transcript_gtf.pl --variant VARIANT_FILE --gtf GTF_FILE [--fasta FASTA_FILE\n --chr_num APPRIS_NUM]";
}


## Updated sequence can be printed; Sequence will hashed to check reference allele
my (%seqs, %num); 
if (open ( FASTA, $fasta_file ) ) { 
#my $output_file = $fasta_file . "Modified"; 
#open (OUT, ">$output_file") || die " Cannot open OUTfile: $!\n "; 
    while (<FASTA>) { 
	if (/>(ENST\d+\.\d)/) { 
	    $seqs{"$1"} = <FASTA>; 
	}
	else  { 
	    print "No transcript ID in fasta file\n"
	} 
    }
}

if (open (NUM, $num_file) ) { 
    while (<NUM>) { 
	my @T = split; 
	if ($T[0] =~ /(ENST\d+\.\d)/) { 
	    $num{$1} = $T[1];
	}
	else { 
	    print "Incorrect ID in num_file\n"; 
	}
    }
}

## Go through GTF file and keep track of lengths of all exons
## DATASTRUCT: key(ENST_Strand) value [exon_lengths]
## # We can use the exon count to keep track of the position
open (GTF, $gtf_file ) || die "Cannot open GTF file\n"; 
#my (%trans_len, %trans_exon_counts); 
my (%trans_len); 
while (<GTF> ) { 
    my $line = $_; 
    my @F = split; 
    if ($F[2] eq 'exon' ) { 
	my $length = $F[4] - $F[3]; 
	if ($line =~ /exon_number (\d+)/) { 
	    my $exon_num = $1; 
	    if ($line =~ /transcript_id "(ENST\d+\.\d)/ ) { 
		my $key = $1 . "$F[6]"; 
		push (@{ $trans_len{$key}}, $length+1 ) ; 
		my $sum1 = sum (@{$trans_len{$key}});
#		print "$sum1\n";
	    }
	    else {
		print "No transcript ID\n"; 
	    }
	}
	else { 
	    print $line; 
	}
    }
} 

#print Dumper(\%trans_len); 
# print "###############\n";
# print Dumper(\%trans_exon_counts); 


## WE CAN ASSUME THAT TRANSCRIPT IDS ARE UNIQUE. 
## KEEP STRAND AND EXON LENGTH INTO ACCOUNT FROM THE GTF
# We can use the exon count to keep track of the position
## Filter out positions with no variation
# USE ONLY SNPS
# Each varaitn is stored in a hash. 
# Key is Transcript_ID_Strand; Value is an array of hashes
# For each Variant; We use position on transcript as key; [REF, ALT] pair as value; The position will be calculated from data on each transcript
# SNPs are given on the genome for - strand transcript need to reverse the SNP to find the match. 
# BED coordinates are 0-based GTFs are 1-based. 
# For negative strand take exon_stop - Bed_Stop and reverse SNP to get distance
# chr1 HAVANA exon 894595 894670  -  ... chr1 894603 894604 C means that 42nd nucleotide of this exon is a G
open ( VAR, $variant_file ) || die "Cannot open VARIANT file\n" ; 
my %variants;
while (<VAR>) { 
    my ($key, $dist, $sum); 
    my $line = $_;
    my @G = split; 
    if (check_variation_snp ($G[-2])) { 
#	print "$G[-2]\n"; 
# EXONS LENGTHS ARE ALWAYS ORDERED. GIVEN VARIANT IN EXON N
# SUM LENGTHS OF EXONS i=1 to N-1 and add distance of variant to Exon N start depending on strand
	if ($line =~ /exon_number (\d+)/) {
            my $exon_num = $1 - 2;
            if ($line =~ /transcript_id "(ENST\d+\.\d)/ ) {
		$key= $1 . "$G[6]";
		unless ( exists $trans_len{$key}) { print "ID not found\n"}
		if ($G[6] eq '+') {
		    $dist = $G[-3]  - $G[3];  
		}
		elsif ($G[6] eq '-') { 
		    $dist = $G[4] - $G[-3];  
		}
		else { 
		    print "No strand information\n"; 
		}
		if ($exon_num < 0 ) {
		    $sum = $dist; 
		}
		else {
		    $sum = sum (@{$trans_len{$key}}[0..$exon_num]) + $dist ; 
		}
#		check_ref_allele ($key, $sum, $G[-2]);
#		out_1ksnp ($key, $sum, $G[-2]);
		print "$key\t$sum\t$G[-2]\n";
            }
            else {
                print "No transcript ID\n";
            }
	}
	else {
	    print "No Exon ID\n";
	}
    }
}

sub out_1ksnp {
    open (OUT, ">>1Ksnp") || die "Cannot open: $!\n";
    # 1 72515 C A 0.5 99 1 rs4030300
    my $key = shift;
    my $sum = shift;
    my $variant = shift;
    my $variant_id = "rs"; 
    my $alt_allele; 
    my $strand = chop($key);
    my @Q = split (/_/, $variant);
    my $ref_al = shift(@Q); 
# Identify Variant ID: The genomes that have the variant
# Identify Alternative Alleles
# A_NA18526=G|G_NA18951=G|G_NA12890=G|A_NA10847=G|G_NA18505=G|G
    foreach (@Q) { 
	my @W = split(/=/); 
	my @alleles = split (/\|/, $W[1]); 
	foreach (@alleles) { 
#	    print "$ref_al\t$_\n";
	    if ("$_" eq "$ref_al") {
#		print "YES!\n";
		next; 
	    }
	    else { 
		$variant_id .= "_$W[0]"; 
		$alt_allele = $_; 
	    }
	}
    }
    if (exists $num{$key}) { 
	if ($strand eq '+') { 
	    print OUT "$num{$key}\t$sum\t$ref_al\t$alt_allele\t0.5\t99\t1\t$variant_id\n"
	}
	elsif ($strand eq '-') {
#	    print "$ref_al\t$alt_allele\n"; 
	    $alt_allele =~ tr/ACGT/TGCA/; 
	    $ref_al =~ tr/ACGT/TGCA/;
#	    print "$ref_al\t$alt_allele\n"; 
	    print OUT "$num{$key}\t$sum\t$ref_al\t$alt_allele\t0.5\t99\t1\t$variant_id\n"; 
	}
	else { 
	    print "Incorrect strand information\n"; 
	}
    }
    else {
	print "Cannot find the transcript number\n"; 
    }
    close OUT;
}

sub check_ref_allele { 
    my $key = shift; 
    my $sum = shift; 
    my $variant = shift; 
    my $strand = chop($key);
    my @Q = split (/_/, $variant);
    if (exists $seqs{$key}){
	my $pred = substr ($seqs{$key}, $sum, 1);  
	print "$pred\t$Q[0]\t$strand\t$sum\n";
    }
    else{ 
	print "Cannot find transcript id in fasta file\n"; 
    }
}

sub check_variation_snp { 
    my $value = shift; 
    my @M = split (/_/, $value); 
    my $ref_allele = shift (@M);
#    print "$ref_allele\n"; 
    if (length($ref_allele) > 1 ) { 
	return 0; 
    }
    foreach (@M) { 
	if (/=(.{1})\|(.{1})/) {
	    if ($1 ne $ref_allele || $2 ne $ref_allele) { 
		return 1; 
	    }
	}
	else { 
#	    print "$ref_allele\t$_\n"; 
	    return 0; 
	}
    }
}
