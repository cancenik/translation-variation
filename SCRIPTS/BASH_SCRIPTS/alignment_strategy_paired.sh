#!/bin/bash
#$ -cwd
#$ -l s_rt=24:00:00,h_rt=24:00:00
## CC ## 

source $HOME/.bash_profile

## Modified Alignment pipeline to allow paired end mapping
## Additional commandline parameters:
# $first input_fastq_Pair1.gz
# $second input_fastq_Pair2.gz  
## We want to increase the distance between the pairs to maximize concordance of mapping
## For mRNA (APPRIS)-- Fabian RNAseq has a good tradeoff at -X 2000
## Fabian RNAseq no adapter. Commented out cutadapt
## Add filter for unaligned reads to reduce file size
## Remove -L 18 to increase speed
## Added --no-mixed to rRNA 
## Added APPRIS counting using bedtools+ samtools
## Removed -a parameter from trans remaining alignment as it causes huge increase in run time and file size
## MAIN ALIGNMENT STRATEGY IS FOR THE RIBOSEQ. ALL UNUSED SECTIONS REMOVED FOR SIMPLICITY

## May 4: Current Alignment Pipeline is to the human Reference Genome
## We should think about how to incorporate variants in the alignment to avoid reference bias

## Apr 2 2013 ## 
## Update: ECHO COMMANDLINE ARGUMENTS TO LOG ##

## Mar 19 2013 ## 
#################

################# ALIGNMENT
## Order of Alignment rRNA from UCSC_hg19, transcriptome_v15_Gencode_Protein_Coding_Appris_Principal
# I can add lincRNAs after APPRIS, mtrRNA after rRNA, and EBV at the last step
echo "bowtie2 -1 $first -2 $second  -X 500 --no-mixed  -x rRNA  --no-unal --un-conc-gz Unaligned_rRNA 2>>LOG_$first >/dev/null" >>LOG_$first
bowtie2 -1 $first -2 $second  -X 500 --no-mixed  -x rRNA  --no-unal --un-conc-gz Unaligned_rRNA 2>>LOG_$first >/dev/null
mv Unaligned_rRNA.1 Unaligned_rRNA.1.fastq.gz
mv Unaligned_rRNA.2 Unaligned_rRNA.2.fastq.gz
echo "bowtie2 -X 2000  -x APPRIS -1 Unaligned_rRNA.1.fastq.gz -2 Unaligned_rRNA.2.fastq.gz --no-unal -S appris_aln.sam --un-conc-gz Unaligned_Appris 2>>LOG_$first" >>LOG_$first
bowtie2 -X 2000  -x APPRIS -1 Unaligned_rRNA.1.fastq.gz -2 Unaligned_rRNA.2.fastq.gz --no-unal -S appris_aln.sam --un-conc-gz Unaligned_Appris 2>>LOG_$first
mv Unaligned_Appris.1 Unaligned_Appris.1.fastq.gz
mv Unaligned_Appris.2 Unaligned_Appris.2.fastq.gz
echo "bowtie2 -X 1000  -x PC_Trans -1 Unaligned_Appris.1.fastq.gz -2 Unaligned_Appris.2.fastq.gz --no-unal --un-conc-gz Unaligned_Trans 2>>LOG_$first | samtools view -bS - >trans_aln.bam" >>LOG_$first
bowtie2 -X 1000  -x PC_Trans -1 Unaligned_Appris.1.fastq.gz -2 Unaligned_Appris.2.fastq.gz --no-unal --un-conc-gz Unaligned_Trans 2>>LOG_$first | samtools view -bS - >trans_aln.bam
mv Unaligned_Trans.1 Unaligned_Trans.1.fastq.gz
mv Unaligned_Trans.2 Unaligned_Trans.2.fastq.gz
echo "bowtie2 -X 5000  -x hg19 -1 Unaligned_Trans.1.fastq.gz -2 Unaligned_Trans.2.fastq.gz --no-unal --un-conc-gz Unaligned_Genome 2>>LOG_$first | samtools view -bS - >genome_aln.bam" >>LOG_$first
bowtie2 -X 5000  -x hg19 -1 Unaligned_Trans.1.fastq.gz -2 Unaligned_Trans.2.fastq.gz --no-unal --un-conc-gz Unaligned_Genome 2>>LOG_$first | samtools view -bS - >genome_aln.bam
#################

## FASTQC APPRIS ALIGNMENTS
#echo "fastqc -f appris_aln.sam --noextract" >>LOG_$first
#fastqc -f appris_aln.sam --noextract

############# READ COUNTS & ANALYSIS 
## Convert SAM OUTPUT TO BAM FILE CUTOFF WITH -q
samtools view -bS -q 2 appris_aln.sam >appris_aln_Q2.bam

## Count reads per feature using BedTools force strandedness
## Need to specify --split when junction mapping is possible
coverageBed -abam appris_aln_Q2.bam  -b /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris_Regions.bed -s -split >Transcript_Counts_RNASEQ

