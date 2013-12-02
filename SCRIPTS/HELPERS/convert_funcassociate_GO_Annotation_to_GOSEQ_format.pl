#!/usr/bin/perl -w 

## CC ## 
## Aug 15 2013 ## 

# Potential future improvements
# ID is given as regex from commandline, files from commandline

#################
# This script converts funcassociate input into goseq input
# FUNCASSSOCIATE EX.
# GO:0000002mitochondrial genome maintenanceUC007AFD.1 UC007AFF.1 UC007ZTE.1 UC008WXF.1 UC009HIF.1 UC009HIG.1 UC009HIH.1 UC009HII.1 UC009HIJ.1 UC009HIL.1 UC009IZL.1
# GOSEQ
# UC007AFF.1 GO:0000002, GO:0000003
#################

use strict; 
require Simple_Utils;
my $funcassociate = "funcassociate_go_associations.txt"; 
my $output_file = "goseq_input_gene_go_associations.txt"; 

open (FUNC, $funcassociate) || die "Cannot open $funcassociate: $!"; 
open (OUT, ">$output_file") || die "Cannot open $output_file: $!"; 

while (<FUNC>) { 
    my @F = split; 
    foreach (@F) { 
	if (/UC\d+/i) {
	    my $ucsc_id = lc;
	    print OUT "$ucsc_id\t$F[0]\n"
	}
    }
}

