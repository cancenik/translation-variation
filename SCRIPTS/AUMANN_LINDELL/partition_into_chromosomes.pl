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

MAIN_LOOP: while(my $line = <$IH>){
  my @lineArray = split(" ", $line);
  next MAIN_LOOP if(scalar(@lineArray) < 4);
  
  if(exists $chromosomes{$lineArray[0]}){
    my $IFILE =  $chromosomes{$lineArray[0]};
    print $IFILE $line; 
  }
  else{
    my $fileName = $outputFolder."/".$lineArray[0].$fileExtension;
    local *FILE;
    open(FILE, ">".$fileName) or die("Can not open the file".$fileName);
    print FILE $line;
    $chromosomes{$lineArray[0]} = *FILE;
  }
}

foreach my $key (keys %chromosomes){
  close($chromosomes{$key});
}

close($IH);


__END__

=head1 NAME

partition_into_chromosomes.pl

=head1 SYNOPSIS

This program requires two arguments.
The first one is the input file. It is in bedGraph format.
The second one is the ouptput directory.
It partitions the given bedGraph file acording to the chromosomes into smaller files in the output directory.
Run this script as

./separate_Wigs.pl -i <input file> -o <output dir>

For detailed information type

./separate_Wigs.pl --help


=head1 DESCRIPTION

If the reads for all chromosomes are in one file, this script partitions this big file
into smaller files according to the chromosomes.
So in the end, each file contains the reads of a single chromosome.
The output files are named according  to the first column (chromosomes).
For example chr1.bg holds the read counts of the first chromosome.

=head1 EXAMPLE

ex: partition_into_chromosomes.pl -i inpput.bg /home/user/outdir

=head1 AUTHORS

 Can Cenik
 
 Hakan Ozadam

 
=head1 FILE FORMAT

A bedGraph file (both input and output) has four entries in each row separated by whitespace.

Column 1: chromosome
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
