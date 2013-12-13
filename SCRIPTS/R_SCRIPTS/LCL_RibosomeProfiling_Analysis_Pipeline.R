library("edgeR")
library("hexbin")
library("limma")
library("qtl")
library ("sva")
library ("isva")

## NOTES
# Number of bases covered is not a very good measure
# The number of species correlates better with length but
# Total read number gives tighter clustering 
# Should use total read number for the analysis

## DATA INPUT
# RData all_rnaseq DF
load('~/project/CORE_DATAFILES/All_RNA_Seq_DF.Rdata')
rna_seq_normalized <- read.table('~/project/CORE_DATAFILES/TMM_VarianceMeanDetrended_CPM_GT1_in40_BatchRemoved_RNASeq_Expression', header=T)

## HGNC_to_ENSG
hgnc_to_ensg <- read.table('~/project/CORE_DATAFILES/HGNCtoENSG.txt', ,header=F,as.is=T,sep="|",fill=T)
ensg_hgnc <- cbind(grep("ENSG", unlist(strsplit(hgnc_to_ensg$V2, "[.]")), value=T), hgnc_to_ensg$V5)
colnames(ensg_hgnc) <- c("ENSG", "HGNC")

## Absolute Protein Amounts
protein_absolute_ibaq <- read.csv('~/project/CORE_DATAFILES/TableS8_Khan_etal.csv')
protein_absolute_ibaq <- merge(protein_absolute_ibaq, ensg_hgnc)

## RNA_SEQ COUNTS 
## GEUVADIS, PICKRELL, PolyA, Ribozero
geuvadis <- read.table(
"/srv/gs1/projects/snyder/ccenik/LCL_RNASEQ/GEUVADIS/Reformatted_Transcript_Counts_All_Libraries.tsv"
, header=T)
pickrell <- read.table(
"/srv/gs1/projects/snyder/ccenik/PICKRELL_RNASEQ/Reformatted_Transcript_Counts_All_Libraries.tsv"
, header=T)
polyA <- read.table(
"/srv/gs1/projects/snyder/ccenik/LCL_RNASEQ/PolyA_RNA/Reformatted_Transcript_Counts_All_Libraries.tsv"
, header=T)
ribozero <- read.table(
"/srv/gs1/projects/snyder/ccenik/LCL_RNASEQ/RIBOZERO_RNA/Reformatted_Transcript_Counts_All_Libraries.tsv"
, header=T)

# For compatibility, use CDS counts only
pickrell_CDS <- split(pickrell, pickrell$REGION)[[2]]
polyA_CDS <- split(polyA, polyA$REGION)[[2]]
ribozero_CDS <- split(ribozero, ribozero$REGION)[[2]]
geuvadis_CDS <- split(geuvadis, geuvadis$REGION)[[2]]
all_rnaseq <- cbind( polyA_CDS[,grep("Counts", colnames(polyA_CDS))], 
ribozero_CDS[,grep("Counts", colnames(ribozero_CDS))], pickrell_CDS[,grep("Counts", colnames(pickrell_CDS))], 
geuvadis_CDS[,grep("Counts", colnames(geuvadis_CDS))])
# Spearman correlation is quite high -> Add colname_suffix
colnames(all_rnaseq)[1:18] <- paste (colnames(all_rnaseq)[1:18], "polyA", sep="_")
colnames(all_rnaseq)[19:44] <- paste (colnames(all_rnaseq)[19:44], "RiboZero", sep="_")
colnames(all_rnaseq)[45:68] <- paste (colnames(all_rnaseq)[45:68], "Pickrell", sep="_")
colnames(all_rnaseq)[69:82] <- paste (colnames(all_rnaseq)[69:82], "Geuvadis", sep="_")
all_rnaseq_counts <- DGEList(counts= all_rnaseq)

pickrell_CDS_counts <- DGEList(counts=pickrell_CDS[,grep("Counts", colnames(pickrell_CDS))])
polyA_CDS_counts <- DGEList(counts=polyA_CDS[,grep("Counts", colnames(polyA_CDS))])
ribozero_CDS_counts <- DGEList(counts=ribozero_CDS[,grep("Counts", colnames(ribozero_CDS))])
geuvadis_CDS_counts <- DGEList(counts=geuvadis_CDS[,grep("Counts", colnames(geuvadis_CDS))])

## RIBOSEQ_COUNTS
#READ COUNTS per appris transcript
# Custom Count all species
species_all <- read.table("Reformatted_Appris_Species_Counts", header=T)
# CDS SPECIES
species <- read.table("Reformatted_Species_Counts_All_Libraries.tsv", header=T)
CDS_species <- split (species, species$REGION)[[2]]

# Total Number of Reads
#data <- read.table ("Reformatted_Transcript_Counts_All_Libraries.tsv", header=T)
data <- read.table ("~/project/CORE_DATAFILES/Reformatted_Transcript_Counts_All_Libraries.tsv", header=T)
CDS <- split (data, data$REGION)[[2]]
UTR3 <- split (data, data$REGION)[[3]]
UTR5 <- split (data, data$REGION)[[4]]
CDS_Counts <- CDS[,grep("Counts", colnames(CDS))]
CDS_Coverage <- CDS[, grep("CoveredBases", colnames(CDS))]
CDS_Len <- CDS[,grep("Len", colnames(CDS))]
CDS_IDs <- CDS[,1]

#covariates <-  read.table ("Sequenced_Ribosome_Profiling_Sample_Information_Batch_Effects.tsv", header=T)
covariates <-  read.table ("~/project/CORE_DATAFILES/Sequenced_Ribosome_Profiling_Sample_Information_Batch_Effects.tsv", header=T)

## DATA ANALYSIS 

# TOTAL COUNTS
colSums(CDS_Counts[,-c(1,2)])
dim (CDS_Counts[keep(CDS_Counts[,-c(1,2)], 100),])
cm <- cor(CDS_Counts[keep(CDS_Counts[,-c(1,2)], 100),-c(1,2)])
dd <- dist (t(log10(CDS_Counts[keep(CDS_Counts[,-c(1,2)], 100),-c(1,2)]+1)) )
hc <- hclust (dd, "ward") 
hc <- hclust (dd)

# SPECIES
colSums (species[, -c(1)])
dim(species[keep(species[, -c(1)], 100), ])
cm <- cor(species[keep(species[, -c(1)], 100), -c(1)])
dd <- dist(t (species[keep(species[, -c(1)], 100), -c(1)] ) ) 
hc <- hclust(dd, "ward")
hc <- hclust(dd)

## SPECIES TO COUNTS COMPARISON
m1_counts <- CDS_species[,c(1,grep("Counts", colnames(CDS_species)))]
m1 <- merge(m1_counts, CDS, by="ID")
m1_species_sum <- rowSums(m1[,2:35])
cds_count_sum <- rowSums(m1[,grep("Counts", colnames(m1))])
ratios <- cds_count_sum/m1_species_sum
ratios_ids <- data.frame(ID=as.vector(m1[,1]), ratio=ratios)
ratios_dataframe <- data.frame(ID=as.vector(m1[,1]), ReadCount=log10(cds_count_sum+1) , SpeciesCount=log10(m1_species_sum+1) )
# Perform loess regression between read_count to species_count
# Call outliers as 2*SE away from the fit-- Need to think about length here. 
# A very short gene with very high expression might be expected to have a skewed ratio
## This needs a lot of memory
count_to_species <- predict(loess(ReadCount~SpeciesCount, data=ratios_dataframe, statistics="approximate", trace.hat="approximate"), se=T)
c1 <- count_to_species$fit+10^2*count_to_species$s
length(which(ratios_dataframe$ReadCount- c1 > 0))
# c2 <- count_to_species$fit-10^2*count_to_species$s
# plot(ratios_dataframe$ReadCount, ratios_dataframe$SpeciesCount, pch=19, cex=0.2)
#lines(ratios_dataframe$ReadCount,count_to_species$fit, col="red")
# lines(ratios_dataframe$ReadCount,count_to_species$fit+3*count_to_species$s, lty=2)
# lines(ratios_dataframe$ReadCount,count_to_species$fit-3*count_to_species$s, lty=2)


# RNASEQ NORMALIZATION AND VOOM
rnaexpr <- rowSums(cpm(all_rnaseq_counts) > 1) >= 40
all_rnaseq_counts <- all_rnaseq_counts[rnaexpr,]
all_rnaseq_counts <- calcNormFactors (all_rnaseq_counts, method= "TMM")
all_rnaseq_counts$samples
cor (all_rnaseq[rnaexpr,], method="spearman")
pdf("Mean_Variance_Modelling_RNASEQ.pdf")
v2 <- voom (all_rnaseq_counts, plot=T)
dev.off()
# BATCH CORRECTION IS ESSENTIAL HERE
rnaseq_batch <- c (rep(1,18), rep(2, 26), rep(3,24), rep(4, 14) ) 
sample_id <- unlist(strsplit(colnames(v2$E), split= "_"))
sample_id <- sample_id[grep("GM", sample_id)]
mod <- model.matrix(~as.factor(sample_id))
batch_removed <- ComBat (v2$E, batch=rnaseq_batch, mod=mod)

# Update v2$E with batch_removed; this keeps weights
v2$E <- batch_removed

write.table(batch_removed,
file ="TMM_VarianceMeanDetrended_CPM_GT1_in40_BatchRemoved_RNASeq_Expression",
sep="\t", row.names=geuvadis_CDS[rnaexpr,1])

norm_dd_rnaseq <- dist( t( v2$E) )
hc_rnaseq <- hclust( norm_dd_rnaseq)

batch_dd_rnaseq <- dist( t( batch_removed) )
hc_batch_rnaseq <- hclust( batch_dd_rnaseq)

# VOOM- TMM Normalization- Ribosome Profiling
cds_counts <- DGEList(counts=CDS_Counts)
isexpr <- rowSums(cpm(cds_counts) > 1) >= 36
cds_counts <- cds_counts[isexpr,]
cds_counts$samples$lib.size <- colSums(cds_counts$counts)
cds_counts <- calcNormFactors (cds_counts, method= "TMM")
cds_counts$samples
#apply(cds_counts$counts, 2, quantile)
s1 <- unlist(strsplit(colnames(CDS_Counts), "_"))
sample_labels <- s1[grep('GM', s1)]
table(sample_labels)
design <- model.matrix(~sample_labels)
pdf("Mean_Variance_Modeling_voom.pdf")
v <- voom(cds_counts,design, plot=T)
dev.off()
norm_dd <- dist(t (v$E ) ) 
norm_hc <- hclust (norm_dd)
replicate_present <- duplicated(sample_labels) |duplicated(sample_labels, fromLast = TRUE)
norm_dd_rep <- dist(t (v$E[,replicate_present] ) )
norm_hc_rep <- hclust (norm_dd_rep)

detrended_matrix <- cbind(as.vector(CDS_IDs[isexpr]), as.data.frame (v$E) )
colnames(detrended_matrix)[1] <- "ID"
detrended_matrix_with_ratio <- merge (detrended_matrix,ratios_ids , by="ID")

write.table(detrended_matrix_with_ratio, 
file ="TMM_VarianceMeanDetrended_CPM_GT1_in20_RiboSeq_Expression_Read_Species_Ratio", 
sep="\t", row.names=F) 
write.table(v$E, file ="TMM_VarianceMeanDetrended_CPM_GT1_in20_RiboSeq_Expression", 
row.names=CDS_IDs[isexpr], sep="\t")
# Added comment
# SVA/Batch Correction
p1 <- paste (covariates$Date_Cells_Frozen, covariates$Date_Ribosome_Footprint, covariates$Data_Gel_Purified, sep="_")
p1 <- covariates$Data_Gel_Purified
p2 <- paste (covariates$Sequencing_ID, "Counts", sep="_")
p1 <- p1[p2 %in% colnames(v$E)]
p2 <- p2[p2 %in% colnames(v$E)]
p1 <- p1[sort (p2, index.return=T)$ix]
mod <- model.matrix(~as.factor(sample_labels), data=as.data.frame(v$E))

svobj <- sva (v$E, mod=mod, B=20)
fit = lmFit(v$E, svobj$sv)
norm_expr <- residuals (fit, v$E)
v$E <- norm_expr
# Replace v$E with norm_expr

sva_dd <- dist( t ( norm_expr[,replicate_present]))
sva_hc <- hclust (sva_dd) 

#
# Compare absolute levels of protein with rna and ribo
grand_mean_rna <- apply (rna_seq_normalized, 1, median)
grand_mean_rna  <- data.frame(HGNC=rownames(rna_seq_normalized), grand_mean_rna)
grand_mean_ribo <- apply(norm_expr, 1, median)
grand_mean_ribo <- data.frame (HGNC=CDS[isexpr,1], grand_mean_ribo)
CDS_Lens <- data.frame(HGNC=CDS_IDs, CDS_Len[,1])
merge_ribo_prot <- merge(grand_mean_ribo,protein_absolute_ibaq, by="HGNC" )
merge_ribo_rna_prot <- merge (merge_ribo_prot, grand_mean_rna, by="HGNC")
merge_ribo_rna_prot_len <- merge(merge_ribo_rna_prot, CDS_Lens, by="HGNC")
dim(merge_ribo_rna_prot)
rna_cor <- cor.test(merge_ribo_rna_prot$grand_mean_rna, log10(merge_ribo_rna_prot$ibaq.human))
ribo_cor <- cor.test(merge_ribo_rna_prot$grand_mean_ribo, log10(merge_ribo_rna_prot$ibaq.human))
}

#

# Differential Expression Analysis and Translation Efficiency
# Need to merge normalized RNA-Seq with Ribo
# Identify design matrix
# Two predictors: Sample Label + Ribo vs RNA
### NEED TO IMPLEMENT THIS!!! EASIEST WAY IS TO MODIFY ONE OF THE ELIST
# It is possible to cbind ELIST GIVEN THAT THEY ARE SAME ORDER
# EXTRACT VECTOR BINDING BETWEEN THE TWO AND CBIND THE VOOM ELISTs
#### We should bring the voom-derived weights to this calculation
norm_expr_with_ids <- data.frame(HGNC=CDS[isexpr, 1], norm_expr)
rna_seq_normalized_with_ids <- data.frame(HGNC=rownames(rna_seq_normalized), rna_seq_normalized)
all_rna_ribo <- merge (norm_expr_with_ids, rna_seq_normalized_with_ids, by="HGNC" )

# intersect ids
common_ids <- intersect(rownames(rna_seq_normalized), CDS_IDs[isexpr])

sample_id_all <- unlist(strsplit(colnames(all_rna_ribo), split= "_"))
sample_id_all <- as.factor(sample_id_all[grep("GM", sample_id_all)])
treatment <- as.factor(c(rep("Ribo", dim(norm_expr)[2]), rep("RNA", dim(rna_seq_normalized)[2] ) ) )
# Include an interaction term
#factor_all <- as.factor(paste (treatment, sample_id_all, sep="."))
#design <- model.matrix(~0+factor_all)
#colnames(design) <- levels(factor_all)
design <- model.matrix(~treatment+treatment:sample_id_all)
fit <- lmFit(all_rna_ribo[,-1], design)

# One approach is to have two different linear models
# One for different RNA, one for different Ribo and one for different Ribo-RNA
# This way the moderated F-statistic can be used as the measure of any difference
# This might be achieved by the subset argument to lm
# design2 <- model.matrix(~sample_id_all)
# fit2 <- lmFit(all_rna_ribo[,-1], design2, subset= Treatment=="Ribo")
# To achieve ribo rna contrast only, we can use another design with just treatment

# EBayes also returns a moderated F-statistic, $F.p.value
# Test for presence of any nonzero contrast 
## FIGURES
# Absolute Protein to Grand_Mean_RNA or RIBO
pdf ("Absolute_Protein_IBAQ_RNA_RIBO.pdf", width=4, height=8)
par (mfrow = c(2,1))
plot(merge_ribo_rna_prot$grand_mean_ribo, log10(merge_ribo_rna_prot$ibaq.human), 
     pch=20, cex=0.3, xlab="Mean_Ribosome_Profiling_Expression", 
     ylab="Log10_iBaq_ProteinExpression")
legend(0, 8.5 , paste("Pearson Cor", round(ribo_cor$estimate, 3), sep=": "), bty="n" )
plot(merge_ribo_rna_prot$grand_mean_rna, log10(merge_ribo_rna_prot$ibaq.human), 
     pch=20, cex=0.3, xlab="Mean_RNASeq_Expression", 
     ylab="Log10_iBaq_ProteinExpression")
legend(0, 8.5 , paste("Pearson Cor", round(rna_cor$estimate, 3), sep=": "), bty="n" )
dev.off()

# Species to Counts correlation
pdf("Species_Counts_log10.pdf")
plot(log10(cds_count_sum+1), log10(m1_species_sum+1), pch=19, cex=0.4)
dev.off()

# All-by_all correlation plot becomes infeasible
pdf("LCL_Ribo_Pairs_Normalized.pdf")
pairs(v$E[,replicate_present], diag.panel=panel.hist, upper.panel=NULL, lower.panel=NULL)
dev.off()

pdf ("LCL_Ribo_MDS_RepPresent.pdf")
plotMDS(v$E[,replicate_present], labels=sample_labels[replicate_present])
dev.off()

pdf ("LCL_Ribo_Hierarchical_Clustering_Normalized_Counts.pdf", height=11, width=9)
plot(norm_hc)
dev.off()

pdf("Distribution_of_Normalized_RiboExpression.pdf")
hist(norm_expr,100)
dev.ofF()


pdf ("LCL_RiboWithReps_Hierarchical_Clustering_Normalized_Counts.pdf", height=11, width=9)
plot(norm_hc_rep)
dev.off()

pdf ("LCL_Ribo_Hierarchical_Clustering_Normalized_Counts_Residuals_SVA.pdf")
plot(sva_hc)
dev.off()

plot(norm_expr[,2], v$E[,2], pch=19, cex=0.3, xlab="SVA_GM12878_Rep1", ylab="GM12878_Rep1")
plot(norm_expr[,21], v$E[,21], pch=19, cex=0.3, xlab="SVA_GM18526_MiSeq", ylab="GM18526_MiSeq")
##

## FUNCTIONS
panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="red", ...)
}

panel.cor <- function(x,y, ...) {
    par(new=TRUE)
    cor_val <- cor.test(x, y, method="spearman", na.rm=T)$estimate
    cor_val2 <- cor.test(x, y, na.rm=T)$estimate
    cor_val <- round(cor_val, digits=2)
    cor_val2 <- round(cor_val2, digits=2)
    legend("center", cex=0.75, bty="n", paste(cor_val, "\n", cor_val2))
}

panel.smoothScatter <- function (x, y, ...) {
    par(new=TRUE)
    smoothScatter(x, y, nrpoints=0)
}

keep <- function(x, n) {rowSums(x) > n}
##


## MISC
# TOTAL COVERED BASES
colSums(CDS_Coverage[,-c(1,2)])
dim (CDS_Coverage[keep(CDS_Coverage[,-c(1,2)], 100),])
cm <- cor(CDS_Coverage[keep(CDS_Coverage[,-c(1,2)], 100),-c(1,2)])
dd <- dist (t(log10(CDS_Coverage[keep(CDS_Coverage[,-c(1,2)], 100),-c(1,2)]+1)) )
hc <- hclust (dd, "ward") 
hc <- hclust (dd)

# Length Correlations
cor.test(CDS_Len[keep(CDS_Counts,1000),1], rowSums(CDS_Counts[keep(CDS_Counts,1000),]))
cor.test(m1[keep(m1[,2:35], 100),39], m1_species_sum[keep(m1[,2:35], 100)])
cor.test(m1[keep(m1[,2:35], 100),39], m1_species_sum[keep(m1[,2:35], 100)], method="spearman")

# Length to Read Count Correlation
pdf ("Length_vs_TotalCDSCoverage.pdf")
plot(log10(CDS_Len[keep(CDS_Coverage,1000),1]), log10(rowSums(CDS_Coverage[keep(CDS_Coverage,1000),])), pch=19, cex=0.4)
dev.off()


