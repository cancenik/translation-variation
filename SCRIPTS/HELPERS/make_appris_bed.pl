#!/usr/bin/perl 

## CC 
## May 17 2013

######################                                                                                                                                                                                    
## CC May 18 2013
## Discovered bug due to Appris annotation.
## The name has UTR3, UTR5 and CDS information.
## If UTR3 is indeed before the CDS; the transcript is on the - strand
## However, UTR3 doesn't mean 3'UTR for transcripts on the - strand
## For - strand transcript UTR3 is indeed 5'UTR
######################                                                                                                                                                                                    

####################
# This script generates bed file for APPRIS annotations.
# The coordinates are with respect to Transcript not genome
# We will generate one region for 5UTR, CDS, 3UTR
####################

my $appris = "/srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/appris_selected_transcripts.fa"; 
open (IN, $appris); 

my $out = "Appris_Regions.bed"; 
open (OUT, ">$out"); 

while (<IN>) {
    if (/>(.*)/) {
	chomp;
	my $id = $1; 
	if ($id =~ /UTR\d:(\d+)-(\d+)\|CDS:(\d+)-(\d+)\|UTR\d:(\d+)-(\d+)/) {
	    print OUT "$id\t$1\t$2\tUTR5\t0\t+\n";
	    print OUT "$id\t$3\t$4\tCDS\t0\t+\n";
	    print OUT "$id\t$5\t$6\tUTR3\t0\t+\n";
	}
	elsif ($id =~ /UTR\d:(\d+)-(\d+)\|CDS:(\d+)-(\d+)/) { 
	    print OUT "$id\t$1\t$2\tUTR5\t0\t+\n";
	    print OUT "$id\t$3\t$4\tCDS\t0\t+\n";
	}
	elsif ($id =~ /CDS:(\d+)-(\d+)\|UTR\d:(\d+)-(\d+)/) { 
	    print OUT "$id\t$1\t$2\tCDS\t0\t+\n";
	    print OUT "$id\t$3\t$4\tUTR3\t0\t+\n";
	}
	elsif ($id =~ /CDS:(\d+)-(\d+)/) { 
	    print OUT "$id\t$1\t$2\tCDS\t0\t+\n";
	}
	else {
	    my $line = $id;
	    @F = split(/\|/, $id);
	    print OUT "$line\t1\t$F[-1]\tAPPRIS\t0\t+\n"; 
	}
    }
    else { 
	next; 
    }
}

