## CC ## 
## Apr 2 2013 ## 
## Update: ECHO COMMANDLINE ARGUMENTS TO LOG ##
## Updated to generate mm9 Pipeline

## Mar 19 2013 ## 
#################
# FILENAMES PASSED THROUGH COMMAND LINE ARGUMENTS. 
# CURRENT INTENDED USE IS SCREEN_QLOGIN_RUN_SCRIPT
# $1 Adapter Sequence
# $2 input_fastq
#################
### Start with FASTQ or SRA FILE
### fastq-dump *.sra
# if [[ $2 =~ .sra ]]
#  then
#   fastq-dump $2
#   echo "File converted from sra to fastq"
#   file="$2.fastq"
#  else
#   file="$2"
# fi 
############### LOW QUALITY READ FILTERING
# Some files might have low quality reads unfiltered. Remove those with fastq_illumina_filter
#fastq_illumina_filter --keep N -v FILE 
#################
################# ADAPTER REMOVAL
# Remove 3' adapter using cutadapt
# -g 5'adapter -a 3'adapter -e ERROR_RATE (0.1)
echo "cutadapt -a $1 --overlap=5 --trimmed-only --minimum-length=18 --quality-cutoff=33 $2 >Clipped_$2 2>>LOG_$2" >>LOG_$2
cutadapt -a $1 --overlap=5 --trimmed-only --minimum-length=18 --quality-cutoff=33 $2 >Clipped_$2 2>>LOG_$2 
# Example
# EMI ADAPTER ATCTCGTATGCCGTCTTCTGCTTG
# cutadapt -a CTGTAGGCACCATCAAT  --overlap=5 --trimmed-only --minimum-length=18 --quality-cutoff=33 
#################
################# ALIGNMENT
## Order of Alignment rRNA from UCSC_hg19, transcriptome_v15_Gencode_Protein_Coding_Appris_Principal
# I can add lincRNAs after APPRIS, mtrRNA after rRNA, and EBV at the last step
echo "bowtie2 -5 1 -L 18 -x mm9_rRNA_tRNA_snRNA -q Clipped_$2 -S rRNA_aln.sam --un Unaligned_rRNA_$2 2>>LOG_$2" >>LOG_$2
bowtie2 -5 1 -L 18 -x mm9_rRNA_tRNA_snRNA -q Clipped_$2 --un Unaligned_rRNA_$2 2>>LOG_$2
echo "bowtie2 --norc -5 1 -L 18 -x mm9_knownCanonical -q Unaligned_rRNA_$2 -S ucKnownCanonical.sam --al Aligned_ucKnownCanonical_$2 --un Unaligned_ucKnownCanonical_$2 2>>LOG_$2" >>LOG_$2
bowtie2 --norc -5 1 -L 18 -x mm9_knownCanonical -q Unaligned_rRNA_$2 -S uc_knownCanonical.sam --al Aligned_ucKnownCanonical_$2 --un Unaligned_ucKnownCanonical_$2 2>>LOG_$2
echo "bowtie2 --norc -5 1 -a -L 18 -x mm9_ucknownGeneAll -q Unaligned_Appris_$2 -S trans_aln.sam --un Unaligned_Trans_$2 2>>LOG_$2" >>LOG_$2
bowtie2 --norc -5 1 -a -L 18 -x mm9_ucknownGeneAll  -q Unaligned_ucKnownCanonical_$2 --un Unaligned_Trans_$2 2>>LOG_$2 | samtools view -bS - > trans_aln.bam
# NCBI rRNA Alignment -- Includes pre-rRNA
# bowtie2 -L 18 -x NCBI_rRNA -q Unaligned_Trans_$2 -S NCBI_rRNA_Aligned.sam 2>>LOG_$2
echo "bowtie2 -5 1 -L 18 -x mm9 -q Unaligned_Trans_$2 -S genome_aln.sam --un Unaligned_Genome_$2 2>>LOG_$2" >>LOG_$2
bowtie2 -5 1 -L 18 -x mm9 -q Unaligned_Trans_$2 --un Unaligned_Genome_$2 2>>LOG_$2 | samtools view -bS - > genome_aln.bam
#################

## FASTQC AFTER ADAPTER REMOVAL
echo "fastqc Clipped_$2" >>LOG_$2
fastqc Clipped_$2
## FASTQC AFTER RRNA ALIGNMENT
echo "fastqc Unaligned_rRNA_$2" >>LOG_$2
fastqc Unaligned_rRNA_$2
## FASTQC TRANS ALIGNMENTS
echo "fastqc Aligned_ucKnownCanonical_$2" >>LOG_$2
fastqc Aligned_ucKnownCanonical_$2

############# READ COUNTS & ANALYSIS 
## Convert SAM OUTPUT TO BAM FILE CUTOFF WITH -q
## OR Keep the results in sam format
#samtools view -bS -q 2 appris_aln.sam >appris_aln_Q2.bam
samtools view -S -q 2 uc_knownCanonical.sam >uc_knownCanonical_Q2.sam


## COUNT FOR TRANSCRIPTS, Remove Reads with a non-zero samtools flag
## specific nucleotide can be specified as well
samtranscriptaln_tocounts.pl --file uc_knownCanonical_Q2.sam --species
#samtranscriptaln_tocounts.pl --file uc_knownCanonical_Q2.sam 

# # ALTERNATIVELY WE CAN FILTER BY SPECIES USING 
# # Then count with bedtools
# samtools view -S -q 2 -H uc_knownCanonical.sam >uc_knownCanonical_header.sam
# sort -u -k10,10  uc_knownCanonical_Q2.sam > uc_knownCanonical_Q2_species.sam
# cat uc_knownCanonical_header.sam uc_knownCanonical_Q2_species.sam >uc_knownCanonical_Q2_species_wheader.sam
# rm uc_knownCanonical_Q2_species.sam
# samtools view -bS uc_knownCanonical_Q2_species_wheader.sam >uc_knownCanonical_Q2_species.bam
# rm uc_knownCanonical_Q2_species_wheader.sam
# coverageBed -abam uc_knownCanonical_Q2_species.bam -b /srv/gs1/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/MOUSE_DATA/UC_Canonical_knownGENE.bed -s -split >Transcript_Counts_$2


## GZIP
gzip *.fastq