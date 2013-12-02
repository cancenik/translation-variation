#!/bin/bash
#$ -cwd
#$ -l s_rt=24:00:00,h_rt=24:00:00

source $HOME/.bash_profile

echo "bowtie2 -L 18 -x rRNA -q Merged_Reads.fastq  --un Unaligned_rRNA 2>>LOG | samtools view -bS - >rRNA_aln.bam" >>LOG
bowtie2 -L 18 -x rRNA -q Merged_Reads.fastq  --un Unaligned_rRNA 2>>LOG | samtools view -bS - >rRNA_aln.bam
echo "bowtie2 -L 18 --norc -x APPRIS -q Unaligned_rRNA -S appris_aln.sam --un Unaligned_Appris 2>>LOG" >>LOG
bowtie2 -L 18 --norc -x APPRIS -q Unaligned_rRNA -S appris_aln.sam --un Unaligned_Appris 2>>LOG
echo "bowtie2 --norc -a -L 18 -x PC_Trans -q Unaligned_Appris --un Unaligned_Trans 2>>LOG | samtools view -bS - >trans_aln.bam" >>LOG
bowtie2 --norc -a -L 18 -x PC_Trans -q Unaligned_Appris --un Unaligned_Trans 2>>LOG | samtools view -bS - >trans_aln.bam
# NCBI rRNA Alignment -- Includes pre-rRNA
# bowtie2 -L 18 -x NCBI_rRNA -q Unaligned_Trans -S NCBI_rRNA_Aligned.sam 2>>LOG
echo "bowtie2 -L 18 -x hg19 -q Unaligned_Trans -S genome_aln.sam --un Unaligned_Genome 2>>LOG" >>LOG
bowtie2 -L 18 -x hg19 -q Unaligned_Trans --un Unaligned_Genome 2>>LOG | samtools view -bS - >genome_aln.bam
#################

############# READ COUNTS & ANALYSIS 
## Convert SAM OUTPUT TO BAM FILE CUTOFF WITH -q
## OR Keep the results in sam format
samtools view -bS -q 2 appris_aln.sam >appris_aln_Q2.bam
#samtools view -S -q 2 appris_aln.sam >appris_aln_Q2.sam

## Count reads per feature using BedTools force strandedness
## Need to specify --split when junction mapping is possible
# coverageBed -abam appris_aln.bam -b features.bed -s -split >COUNTS
coverageBed -abam appris_aln_Q2.bam -b /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris_Regions.bed -s -split >Transcript_Counts
## Extract non-zero with HGNC symbols
# perl -ane 'if ($F[6]>10) {print join("\t", @F); print "\n"}' Transcript_Counts_CoveraBed >test.out

## OR MANUALLY COUNT FOR TRANSCRIPTS, Remove Reads with a non-zero samtools flag
#samtranscriptaln_tocounts.pl --file appris_aln_Q2.sam

## FASTQC AFTER RRNA ALIGNMENT
echo "fastqc Unaligned_rRNA" >>LOG
fastqc Unaligned_rRNA
## FASTQC TRANS ALIGNMENTS
echo "appris_aln_Q2.bam" >>LOG
fastqc appris_aln_Q2.bam

### GZIP FASTQ FILES
gzip *.fastq


