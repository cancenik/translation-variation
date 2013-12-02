#!/bin/bash

## CC 
## Nov 8 2013

## This script uses to count transcripts and gather stats using species
## We will recount species using coverageBED 
## We will treat each different sequence as different species
## One issue might be 
source $HOME/.bash_profile
LCLRIBO=/srv/gs1/projects/snyder/ccenik/LCL_RIBOSEQ
cd $LCLRIBO
for dir in $(find . -maxdepth 1 -mindepth 1 -name 'GM*' -type d -printf '%f\n')
do
#echo $dir
cd $LCLRIBO/$dir
if [ ! -e ./appris_aln_Q2_species.bam ]
 then
 pwd
 samtools view -S -q 2 -H appris_aln.sam >appris_aln_header.sam
 samtools view -S -q 2 appris_aln.sam >appris_aln_Q2.sam 
 sort -u -k10,10  appris_aln_Q2.sam > appris_aln_Q2_species.sam
 cat appris_aln_header.sam appris_aln_Q2_species.sam >appris_aln_Q2_species_wheader.sam
 samtools view -bS appris_aln_Q2_species_wheader.sam >appris_aln_Q2_species.bam
 rm -f appris_aln_Q2_species.sam
 rm -f appris_aln_Q2.sam
 rm -f appris_aln_Q2_species_wheader.sam
 coverageBed -abam appris_aln_Q2_species.bam -b /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris_Regions.bed -s -split >Species_Counts_UTRs_CDS_CoverageBed.tsv
fi
done

## CUSTOM COUNTER was used before
## However, doesn't assign utr, cds counts separetely
#SCRIPT=~/SCRIPTS 
#qsub -l s_rt=24:00:00,h_rt=24:00:00 -cwd -v file="$LCLRIBO/$dir/appris_aln.sam" $SCRIPT/samtranscriptaln_tocounts_qsub.pl