#!/bin/bash
#$ -cwd
#$ -l s_rt=24:00:00,h_rt=24:00:00

source $HOME/.bash_profile

## CC ## 
## June 14 2013: --Change transcriptome alignment to forward (Watson) strand alignment only --norc
## Modified to enable compatibility with qsub and gzip input

## May 13 2013: Added simple modifications to create .bam and .gz when possible
## Current Bedtools count strategy using genomeCov doesn't utilize sam flag correctly. Uses FLAG:16 for example
## However, coverageBed -abam appris_aln.bam -b features.bed -s -split >COUNTS works as expected and gives the same result as
## samtranscriptaln_tocounts.pl

## May 4: Current Alignment Pipeline is to the human Reference Genome
## We should think about how to incorporate variants in the alignment to avoid reference bias

## Apr 2 2013 ## 
## Update: ECHO COMMANDLINE ARGUMENTS TO LOG ##

## Mar 19 2013 ## 
#################
# FILENAMES PASSED THROUGH COMMAND LINE ARGUMENTS. 
# CURRENT INTENDED USE IS SCREEN_QLOGIN_RUN_SCRIPT
# $1 Adapter Sequence
# $2 input_fastq
#################
################# ADAPTER REMOVAL
# Remove 3' adapter using cutadapt
# -g 5'adapter -a 3'adapter -e ERROR_RATE (0.1)

echo "cutadapt -a $1 --overlap=2 --minimum-length=21 --quality-cutoff=33 $2 >Clipped_$2 2>>LOG_$2" >>LOG_$2
cutadapt -a $1 --overlap=2 --minimum-length=21 --quality-cutoff=33 $2 >Clipped_$2 2>>LOG_$2 

# Example
# EMI ADAPTER ATCTCGTATGCCGTCTTCTGCTTG
# cutadapt -a CTGTAGGCACCATCAAT  --overlap=5 --trimmed-only --minimum-length=18 --quality-cutoff=33 
#################
################# ALIGNMENT
## Order of Alignment rRNA from UCSC_hg19, transcriptome_v15_Gencode_Protein_Coding_Appris_Principal
# I can add lincRNAs after APPRIS, mtrRNA after rRNA, and EBV at the last step
echo "bowtie2 -L 18 -x rRNA -q Clipped_$2  --un Unaligned_rRNA_$2 2>>LOG_$2 | samtools view -bS - >rRNA_aln.bam" >>LOG_$2
bowtie2 -L 18 -x rRNA -q Clipped_$2  --un Unaligned_rRNA_$2 2>>LOG_$2 | samtools view -bS - >rRNA_aln.bam
echo "bowtie2 -L 18 --norc -x APPRIS -q Unaligned_rRNA_$2 -S appris_aln.sam --un Unaligned_Appris_$2 2>>LOG_$2" >>LOG_$2
bowtie2 -L 18 --norc -x APPRIS -q Unaligned_rRNA_$2 -S appris_aln.sam --un Unaligned_Appris_$2 2>>LOG_$2
echo "bowtie2 --norc -a -L 18 -x PC_Trans -q Unaligned_Appris_$2 --un Unaligned_Trans_$2 2>>LOG_$2 | samtools view -bS - >trans_aln.bam" >>LOG_$2
bowtie2 --norc -a -L 18 -x PC_Trans -q Unaligned_Appris_$2 --un Unaligned_Trans_$2 2>>LOG_$2 | samtools view -bS - >trans_aln.bam
# NCBI rRNA Alignment -- Includes pre-rRNA
# bowtie2 -L 18 -x NCBI_rRNA -q Unaligned_Trans_$2 -S NCBI_rRNA_Aligned.sam 2>>LOG_$2
echo "bowtie2 -L 18 -x hg19 -q Unaligned_Trans_$2 -S genome_aln.sam --un Unaligned_Genome_$2 2>>LOG_$2" >>LOG_$2
bowtie2 -L 18 -x hg19 -q Unaligned_Trans_$2 -S genome_aln.sam  --un Unaligned_Genome_$2 2>>LOG_$2
#################

############# READ COUNTS & ANALYSIS 
## Convert SAM OUTPUT TO BAM FILE CUTOFF WITH -q
## OR Keep the results in sam format
samtools view -bS -q 2 appris_aln.sam >appris_aln_Q2.bam
#samtools view -S -q 2 appris_aln.sam >appris_aln_Q2.sam

## Count reads per feature using BedTools force strandedness
## Need to specify --split when junction mapping is possible
# coverageBed -abam appris_aln.bam -b features.bed -s -split >COUNTS
coverageBed -abam appris_aln_Q2.bam -b /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris_Regions.bed -s -split >Transcript_Counts_$2
## Extract non-zero with HGNC symbols
# perl -ane 'if ($F[6]>10) {print join("\t", @F); print "\n"}' Transcript_Counts_CoveraBed >test.out

## OR MANUALLY COUNT FOR TRANSCRIPTS, Remove Reads with a non-zero samtools flag
#samtranscriptaln_tocounts.pl --file appris_aln_Q2.sam


##### For Aumann-Lindell Test
### We wrote windowFilter.pl assuming per nucleotide depth
### We are forcing + strand alignment counts only. 
#samtools view -bS -q 2 genome_aln.sam >genome_aln_Q2.bam
#samtools sort genome_aln_Q2.bam sorted
# SPLIT CANNOT BE USED WITH 5'END COUNTING
# bedtools genomecov  -d -5  -strand + -ibam appris_aln_Q2.bam -g /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris_Seq_Lens
# bedtools genomecov  -bg -5  -strand + -ibam appris_aln_Q2_sorted.bam -g /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris_Seq_Lens | gzip >GenomeCoverage_Appris_Q2_PlusStrand_5PEnds_Split.bg.gz 


## FASTQC AFTER ADAPTER REMOVAL
echo "fastqc Clipped_$2" >>LOG_$2
fastqc Clipped_$2
## FASTQC AFTER RRNA ALIGNMENT
echo "fastqc Unaligned_rRNA_$2" >>LOG_$2
fastqc Unaligned_rRNA_$2
## FASTQC TRANS ALIGNMENTS
echo "appris_aln_Q2.bam" >>LOG_$2
fastqc appris_aln_Q2.bam


### GZIP FASTQ FILES
gzip *.fastq


