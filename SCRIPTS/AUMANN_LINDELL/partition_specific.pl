#!/usr/bin/perl -w
#use local::lib;

# Des

use Pod::Usage; 

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage; 

my ($inputFile, $outputFolder) = ("", "");
my $help;

GetOptions( 
	'i=s'   => \$inputFile, #
	'o=s'    => \$outputFolder, 
	'help'      => \$help 
) or die " unrecognized option(s)!! please refer to the help documentation\n\n";

if ($help) {
	pod2usage( {
		'-verbose' => 2,
		'-exitval' => 1,
	} );
}

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($inputFile eq "") or ($outputFolder eq ""));	


open(my $IH, $inputFile) or die("Can not open the input file $inputFile\n");
`mkdir -p $outputFolder`;
opendir(my $TDIRHANDLE, $outputFolder) or die("Can not open the output directory $outputFolder");
close($TDIRHANDLE);
my ($inputFileName, $inputFileDir, $fileExtension) = fileparse($inputFile, qr/\.[^.]*/);
$fileExtension = $fileExtension ? $fileExtension : ".txt";

# print "\$fileExtension = $fileExtension\n","\$inputFile = $inputFile\n","\$outputFolder = $outputFolder\n";

# Holds the file handlers.
# We have a separate file for each chromosome.
my %chromosomes;
my ($currentFile, $openFile) = ("", "");
my $FH;


MAIN_LOOP: while(my $line = <$IH>){
  my @lineArray = split(/\s+/, $line);
  next MAIN_LOOP if(scalar(@lineArray) < 4);
 
   my @firstColumnContents = split(/\|/, $lineArray[0]);
   
   if( scalar(@firstColumnContents) < 5){
    print "Error in the input file!\nThe first column of this line (in the file $inputFile ) is ", $lineArray[0], "\n";
    print "It does not contain five entries separated by | s\n";
    print "The entries of the first column should be separated by | and we take the fifth entry as the identifier.\n";
    print "Aborting...\n";
    exit;
   }
   
  # This is the identifier of each chromosome or transcript.
  # Each entry in the first column should be uniquely identified by its fifth entry (its index is 4 as the array starts from 0). 
    my $identifier = $firstColumnContents[4]; 
    $currentFile = $outputFolder."/".$identifier.$fileExtension;
 
    if( $currentFile ne $openFile ){
      close($FH) if ($FH);
      open($FH, ">".$currentFile) or die("Can not open the file ". $currentFile . " for writing. The error message is \n" . $!);
      $openFile = $currentFile;
    }
    
    print $FH $line;
 
#   if( exists $chromosomes{$identifier} ){
#     my $IFILE =  $chromosomes{$identifier};
#     print $IFILE $line; 
#   }
#   else{
#     my $fileName = $outputFolder."/".$identifier.$fileExtension;
#     local *FILE;
#     open(FILE, ">".$fileName) or die("Can not open the file ". $fileName . " for writing. The error message is \n" . $!);
#     print FILE $line;
#     $chromosomes{$identifier} = *FILE;
#   }
}

foreach my $key (keys %chromosomes){
  close($chromosomes{$key});
}

close($FH) if($FH); 
close($IH);


__END__

=head1 NAME

partition_into_chromosomes_specific.pl

=head1 SYNOPSIS


Partitions the given file into smaller files according to its first column (usually transcript or chromosomeid).
This program is modified for files that have their first column in the form
 
  ENST00000379339.1|ENSG00000131591.13|OTTHUMG00000000745.8|OTTHUMT00000001851.2|C1orf159-002|C1orf159|2429|UTR3:1-211|CDS:212-1354|UTR5:1355-2429|

Each created file can be identified with the first column.
For technical reasons, when the first column is too long or has special characters (such as |),
it is better to identify the file with a shorter string.
In our partucular case, each transcript can be unqiuely identified with the fifth entry
where each entry is separated by | s. So in the above example, the enty  (C1orf159-002, in the above example).
This should be emphasized:

IMPORTANT: We assume that each entry in the first column can be identified by its fifth entry where
           entries are separated by |

 
 
Run this script as
./partition_into_chromosomes_specific.pl -i <input file> -o <output dir>

For detailed information type

./partition_into_chromosomes_specific.pl --help


=head1 DESCRIPTION

This script has been developed to partition bedgraph files according to the first column.
The first column has an identifier of the chromosome or transcript of the read.
If the reads for all chromosomes are in one file, this script partitions this big file
into smaller files according to the chromosomes (indeed first columns).
So in the end, each file contains the reads of a single chromosome or transcript.
The output files are named according  to the first column (chromosomes).
The output files are created inside the output directory.
For example chr1.bg holds the read counts of the first chromosome.

=head1 EXAMPLE

ex: partition_specific.pl -i input.bg /home/user/outdir

=head1 AUTHORS

 Can Cenik
 
 Hakan Ozadam

 
=head1 FILE FORMAT

A bedGraph file (both input and output) has four entries in each row separated by whitespace.
IMPORTNAT: We assume that the identifier has several entries separated by |.
Each column should be identified by its fifth entry!
The results doesn't make sense if this restriction is not satisfied.

Column 1: identifier
Column 2: start position
Column 3: end position
Column 4: read counts

 Ex:

  chr12 5678 5679 12
  chr12 5679 5680 16
 
=head1 LICENSE AND COPYING

 This program is free software; you can redistribute it and / or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.gnu.org/licenses/licenses.html
