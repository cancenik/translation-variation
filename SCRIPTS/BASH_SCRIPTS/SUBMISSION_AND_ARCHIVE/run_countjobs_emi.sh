/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/COFOOT_R1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./cofoot_r1_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/COFOOT_R2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./cofoot_r2_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1OVER_R1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./s1over_r1_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1OVER_R2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./s1over_r2_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1SIRNA_R1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./s1sirna_r1_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1SIRNA_R2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./s1sirna_r2_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_CONT1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./rnaseq_c1_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_CONT2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./rnaseq_c2_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_OVER1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./rnaseq_o1_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_OVER2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./rnaseq_o2_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_SiRNA1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./rnaseq_si1_counts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_SiRNA2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/cds_hg18_renamed.bed >./rnaseq_si2_counts
for f in *_counts
do
    cut -f7 "$f" > "out.$f"
done
cut -f4 rnaseq_c1_counts >names
paste names out* >all_cds_counts
# This part doesn't work because the order is not the same
# Another problem is that the libraries are stranded
# paste ~/EMI_RIBOSEQ/cds_hg18_renamed.bed all_cds_counts >~/EMI_RIBOSEQ/cds_hg18_renamed_all_counts.bed
# Best strategy is to merge in R

/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/COFOOT_R1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./cofoot_r1_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/COFOOT_R2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./cofoot_r2_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1OVER_R1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./s1over_r1_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1OVER_R2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./s1over_r2_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1SIRNA_R1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./s1sirna_r1_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RIBOSEQ/S1SIRNA_R2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./s1sirna_r2_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_CONT1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./rnaseq_c1_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_CONT2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./rnaseq_c2_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_OVER1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./rnaseq_o1_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_OVER2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./rnaseq_o2_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_SiRNA1/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./rnaseq_si1_exoncounts
/srv/gs1/software/bedtools/2.16.2/bin/coverageBed -s  -abam /srv/gs1/projects/snyder/ccenik/EMI_RIBOSEQ/RNASEQ/RNASeq_SiRNA2/tophat/accepted_hits.bam -b ~/EMI_RIBOSEQ/exons_hg18_renamed.bed >./rnaseq_si2_exoncounts
for f in *_exoncounts
do
    cut -f7 "$f" > "out.$f"
done
cut -f4 rnaseq_c1_exoncounts >names_exoncounts
paste names_exoncounts out*exoncounts >all_exon_counts
#paste ~/EMI_RIBOSEQ/exons_hg18_renamed.bed all_exoncounts >~/EMI_RIBOSEQ/exons_hg18_renamed_all_exoncounts.bed