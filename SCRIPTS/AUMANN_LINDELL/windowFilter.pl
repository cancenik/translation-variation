#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);
use File::Basename;
use Class::Struct;
use Benchmark qw(:all) ;
###########################################################################

# Store each line entry in this struct
struct  Entry => [ position => '$', count => '$' ] ;

my %reads;
my ($inputFile, $outPutFolder) = ("", "");
my $minReads = 0;
my $chunkSize = 3;

GetOptions(
  'i=s' => \$inputFile,
  'o=s' => \$outPutFolder,
  'c=i' => \$chunkSize,
  'min=f' => \$minReads,
) or die("Unrecognized optioins.\nFor help, run this script with -help option.\n");


pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) 
    if( ($inputFile eq "") or ($outPutFolder eq "") or ($minReads == 0) );


#print $minReads,"\n"; exit; 
    
#############################################################################

`mkdir -p $outPutFolder`;
opendir(my $TOH, $outPutFolder) or die("Could not open the folder $outPutFolder\n");
close($TOH);

#open(my $IH, $inputFile) or die("Couldn't open the file $inputFile\n"); 
# Get the chromosome name from the input file
my ($chromosome, $inputFileDir, $fileExtension) = fileparse($inputFile, qr/\.([^.]*)/);

my($firstPos, $lastPos) = readFile($inputFile, \%reads);

# foreach my $key (sort {$a <=> $b} keys %reads){
#   print "$key  $reads{$key}\n";
# }


window(\%reads, $minReads, $chromosome, $firstPos, $lastPos, $chunkSize, $outPutFolder);

# We use the notation in aumann-lindell's paper.
# a,b are border values of the output interval where output is made up of intervals of the form [a,b]
# current is a pointer to the given array





###########################################################################
## aumann-lindell window algorithm
sub window{

  my ($reads, $threshold, $chr, $first, $last, $chunkSize, $outDir) = @_;
  
  my ($current, $a, $b) = (0,0,0);
  my ($position, $count) = (0,0);
  my ($Asize, $Bsize, $Asum, $Bsum) = (0,0,0,0);

  
  $current = $first;
  
  my $outBedFile = $outDir."/".$chr.".bed";
  open(my $BEDHANDLE, ">$outBedFile") or die("Couldn't open the file $outBedFile for writing\n");
  
  
  while($current < $last){
    
    ($Bsum, $Bsize) = (0,0);
    ($Asum, $Asize) = (0, 0);
    # Begin with an entry above the threshold

    while( ( getChunk($reads, $current, $chunkSize) < ($threshold * $chunkSize) ) and ($current < $last)    ){ 
      $current += $chunkSize; 
    }
    
    
    $Asize += $chunkSize;
    $Asum = getChunk($reads, $current, $chunkSize);
    
    $a = $current;
    $b = $a + $chunkSize;
    
    while( ( (($Asum + $Bsum) / ($Asize + $Bsize)) >= $threshold) and ($current < $last) ){
      $current += $chunkSize;
      $Bsum += getChunk($reads, $current, $chunkSize);
      $Bsize += $chunkSize;
      
      if( ($Bsum / $Bsize) >= $threshold){
	$Asize += $Bsize;
	$Asum += $Bsum;
	$b = $current + $chunkSize;
	($Bsum, $Bsize) = (0, 0);
      }
    }
    
    
    # Without the last if, we may  report the last entry accidentally.
    if( ($Asum / $Asize) >= $threshold){
    
      print $BEDHANDLE "$chr\t$a\t",$b,"\n";
      my $outBgFile = $outDir."/".$chr."_".$a."_".$b."_.bg";
      open(my $BGHANDLE, ">$outBgFile") or die("Couldn't open the file $outBgFile for writing\n");
      for(my $j = $a; $j < $b; $j++ ){
	print($BGHANDLE "$chr $j ", $j+1, " ", ${$reads}{$j}, " \n") if(exists ${$reads}{$j});  
      }
      close($BGHANDLE);
    }
    
    $current = $b;
  }  
  
  close($BEDHANDLE);
  
}



##############################################################################

# In determining the first and last position,
# We assume that the input data is sorted in ascending order!

# Read the bedGraph file into hash and return the ifirst and last positions as an array.

sub readFile{
  my($file, $hashRef) = @_;
  
  open(my $IH, $file) or die("Can not open the file $file for reading.\n");
  
  my ($firstPosition, $lastPosition) = (-1,-1);
  my @contents;
  
  INITLOOP: while( ($firstPosition < 0) and ($_ = <$IH>)){
    @contents = split(/\s/,$_);
    next READLOOP if(scalar(@contents)) < 4;
   
    $firstPosition = $contents[1];
    for(my $i = $contents[1]; $i < $contents[2]; $i++){
      ${$hashRef}{$i} = $contents[3];
    }
  }
  
  READLOOP: while(<$IH>){
    chomp($_);
    @contents = split(/\s/,$_);
    next READLOOP if(scalar(@contents)) < 4;
    
    for(my $i = $contents[1]; $i < $contents[2]; $i++){
      ${$hashRef}{$i} = $contents[3];
    }
    
  }
  
  $lastPosition = $contents[1];
  
  close($IH);
  
  return ($firstPosition, $lastPosition);
}

#############################################################################

#############################################################################

sub getChunk{
  my($hashRef, $pointer, $chunkSize) = @_;
  
  my $sumResult = 0;
  
  my $currentValue = 0;
  #my $currentPosition 
  
  for(my $i = 0; $i < $chunkSize; $i++){
   $sumResult += (exists(${$hashRef}{$pointer + $i})? ${$hashRef}{$pointer + $i} : 0);
  }
  
  return $sumResult;
}

#############################################################################


__END__


=head1 NAME

windowFilter.pl

=head1 SYNOPSIS

windowFilter.pl -i <input file> -o <output folder> -m <min reads> [-c chunkSize]


windowFilter.pl -help

-c chunksize is optional. Default is 3.

=head1 DESCRIPTION


=head1 EXAMPLE


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
