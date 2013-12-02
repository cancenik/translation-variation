#!/usr/bin/perl -w

### TEST CODE ## 

BEGIN { system ("/bin/bash/source .bash_profile") }

use strict; 
use Getopt::Long;
use List::Util qw(sum);
use List::MoreUtils 'pairwise';

my @G = (1, 2, 3, 12, 24, 44);
my $i = 3;
my $F = sum(@G[0 .. $i-1]); 

print "$F\n";  

my @T = (2,2, 2, 2, 2, 3); 
my @sum = pairwise { $a + $b } @G, @T;
print "@sum\n";
# my  $use_species=  '';
# my $nucleotides_to_consider = 'N';
# my $input_file = '';
# GetOptions ('file=s' =>\$input_file , 'species' => \$use_species, 'nucleotide=s' => \$nucleotides_to_consider);

# print "$use_species\t$nucleotides_to_consider\t$input_file\n";

my $file = $ENV{"file"}; 
print "$file\n";
