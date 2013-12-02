#!/bin/bash 

source $HOME/.bash_profile
for i in /srv/gsfs0/projects/snyder/ccenik/Kotaro_RiboSeq/E*/R* ; do
 cd $i
# echo $i
 samtools view -S -q 2 -H uc_knownCanonical.sam >uc_knownCanonical_header.sam
 sort -u -k10,10  uc_knownCanonical_Q2.sam > uc_knownCanonical_Q2_species.sam 
 cat uc_knownCanonical_header.sam uc_knownCanonical_Q2_species.sam >uc_knownCanonical_Q2_species_wheader.sam
 rm -f uc_knownCanonical_Q2_species.sam
 samtools view -bS uc_knownCanonical_Q2_species_wheader.sam >uc_knownCanonical_Q2_species.bam
 rm -f uc_knownCanonical_Q2_species_wheader.sam
 coverageBed -abam uc_knownCanonical_Q2_species.bam -b /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/MOUSE_DATA/UC_Canonical_knownGENE.bed -s -split >Transcript_Counts_UTRs_CDS.tsv        
done