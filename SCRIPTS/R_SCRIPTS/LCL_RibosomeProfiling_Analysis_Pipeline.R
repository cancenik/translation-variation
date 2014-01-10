library("edgeR")
library("hexbin")
library("limma")
library("qtl")
library ("sva")
library("MASS")
#library ("isva")
# In addition to differential expression analysis at various levels, RNA, Ribo, TE
# We need to add an analysis of variance for the samples with replicates 
# These will reveal gene specific differences in inter-individual variance vs. within individual
# One way to do this is F-test of equality of variances

## NOTES
# Number of bases covered is not a very good measure
# The number of species correlates better with length but
# Total read number gives tighter clustering 
# Should use total read number for the analysis

# Data Directory
data_dir <- '~/project/CORE_DATAFILES/'

## HGNC_to_ENSG
hgnc_to_ensg <- read.table(paste (data_dir, 'HGNCtoENSG.txt', sep=""), ,header=F,as.is=T,sep="|",fill=T)
ensg_hgnc <- cbind(grep("ENSG", unlist(strsplit(hgnc_to_ensg$V2, "[.]")), value=T), hgnc_to_ensg$V5)
colnames(ensg_hgnc) <- c("ENSG", "HGNC")

## Absolute Protein Amounts
protein_absolute_ibaq <- read.csv(paste (data_dir,'TableS8_Khan_etal.csv', sep=""))
protein_absolute_ibaq <- merge(protein_absolute_ibaq, ensg_hgnc)

## RNA_SEQ COUNTS 
## GEUVADIS, PICKRELL, PolyA, Ribozero
geuvadis <- read.table(
  paste (data_dir,"Reformatted_Transcript_Counts_All_Libraries_GEUVADIS.tsv", sep="")
, header=T)
pickrell <- read.table(
  paste (data_dir,"Reformatted_Transcript_Counts_All_Libraries_PICKRELL.tsv", sep="")
  , header=T)
polyA <- read.table(
  paste (data_dir,"Reformatted_Transcript_Counts_All_Libraries_PolyA.tsv", sep="")
  , header=T)
ribozero <- read.table(
  paste (data_dir,"Reformatted_Transcript_Counts_All_Libraries_RZ.tsv", sep="")
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
colnames(all_rnaseq)[45:70] <- paste (colnames(all_rnaseq)[45:70], "Pickrell", sep="_")
colnames(all_rnaseq)[71:86] <- paste (colnames(all_rnaseq)[71:86], "Geuvadis", sep="_")
all_rnaseq_counts <- DGEList(counts= all_rnaseq)


## RIBOSEQ_COUNTS
# Species Counts -- Can this be substituted for Species_Read comparison?
species_all <- read.table(paste (data_dir,"Reformatted_Appris_Species_Counts", sep=""), header=T)
# CDS SPECIES
#species <- read.table(paste(data_dir, "Reformatted_Species_Counts_All_Libraries.tsv", sep=""),header=T)
#CDS_species <- split (species, species$REGION)[[2]]

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

######## DATA ANALYSIS ##################################

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
# Identify differences in PolyA to RiboZero and remove
polyA_mean <- apply(all_rnaseq_counts[,grep("polyA", colnames(all_rnaseq_counts))],1, mean)
RZ_mean <- apply(all_rnaseq_counts[,grep("RiboZero", colnames(all_rnaseq_counts))],1, mean)
plot(log10(polyA_mean), log10(RZ_mean), cex=0.2, pch=19)
fit.rna <- lm(log10(RZ_mean)~log10(polyA_mean))
# Identify outliers
# rna_seq_normalized$ID[abs(stdres(fit.rna)) > 3]
plot(log10(polyA_mean),stdres(fit.rna) , pch=19, cex=0.2)
polyA_RZ_inconsistent <- (abs(stdres(fit.rna)) > 3)
all_rnaseq_counts <- all_rnaseq_counts[!polyA_RZ_inconsistent,]

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
sep="\t", row.names=as.character(rna_seq_normalized_with_ids[,1])[!abs(stdres(fit.rna)) > 3])

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
row.names(v) <- CDS_IDs[isexpr]
save (v, file='~/project/CORE_DATAFILES/RiboProfiling_TMM_Voom_normalizedEList.RData')
# Replace v$E with norm_expr
sva_norm_weights <- v
sva_norm_weights$E <- norm_expr
save (v, file='~/project/CORE_DATAFILES/RiboProfiling_TMM_Voom_normalized_SVA_EList.RData')
sva_dd <- dist( t ( norm_expr[,replicate_present]))
sva_hc <- hclust (sva_dd) 
# We can look sva correlations with sequencing depth, batch, etc

#ISVA seems to remove almost all variance from the data. However, the clustering is very good
isvs <- isvaFn(v$E, sample_labels)
fit2 = lmFit(v$E, isvs$isv)
isva_expr <- residuals(fit2, v$E)
hist(isva_expr)
isva_dd <- dist(t(isva_expr[,replicate_present]))
isva_hc <- hclust(isva_dd)
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


#

# Differential Expression Analysis and Translation Efficiency
# Need to merge normalized RNA-Seq with Ribo
# Identify design matrix
# Two predictors: Sample Label + Ribo vs RNA
#### We should bring the voom-derived weights to this calculation
rna_id_in_ribo <- rownames(rna_seq_normalized) %in% CDS[isexpr, 1]
ribo_id_in_rna <- CDS[isexpr, 1] %in% rownames(rna_seq_normalized)
ribo_rep_present <- unlist(strsplit(colnames(v[,replicate_present]), split= "_"))
ribo_rep_present <- ribo_rep_present[grep("GM", ribo_rep_present)]
rna_with_ribo_replicate <- c()
for (i in 1:length(ribo_rep_present)) { 
  rna_with_ribo_replicate <- c(rna_with_ribo_replicate, grep(ribo_rep_present[i], colnames(v2)))
}


# We can switch v with sva_norm_weights 
## SUBSETTING MESSES UP THE DESIGN MATRIX
#all_expr_elist <- cbind(sva_norm_weights[ribo_id_in_rna,replicate_present], v2[rna_id_in_ribo,unique(sort(rna_with_ribo_replicate))])
all_expr_elist <- cbind(v[ribo_id_in_rna,replicate_present], v2[rna_id_in_ribo,unique(sort(rna_with_ribo_replicate))])
sample_id_all <- unlist(strsplit(colnames(all_expr_elist), split= "_"))
sample_id_all <- as.factor(sample_id_all[grep("GM", sample_id_all)])
all_expr_elist$design <- model.matrix(~sample_id_all)
treatment <- as.factor(c(rep("Ribo", dim(v[,replicate_present])[2]), rep("RNA", dim(v2[,unique(sort(rna_with_ribo_replicate))])[2] ) ) )
treatment <- relevel(treatment,ref="RNA")
# Questions of interest
# What are the genes with differential RNA expression across individuals
# What are the genes with differential Ribo expression across individuals
# What are the genes with differential Ribo-RNA expression across individuals
# What are the genes with differential Ribo-RNA expression in a given individual
# Related to this is the coefficient associate with Ribo-RNA for each gene in each individual
# The moderated F-statistic can be used as the measure of any difference
# ## EBayes also returns a moderated F-statistic, $F.p.value


##### SUBSETTING IS NOT WORKING -- NEED TO FIX
ribo_rna_design <- model.matrix(~sample_id_all)
ribo_fit <- lmFit (all_expr_elist, ribo_rna_design, subset= treatment=="Ribo")
rna_fit <- lmFit (all_expr_elist, ribo_rna_design, subset= treatment=="RNA")
ribo_fit2 <- eBayes(ribo_fit)
rna_fit2 <- eBayes(rna_fit)
topTable(ribo_fit2, coef=2,number=300)
topTable(rna_fit2, coef=2,number=300)
results.ribo <- decideTests(ribo_fit2, p.value=0.01, lfc=1)
results.rna <- decideTests(rna_fit2, p.value=0.01, lfc=1)

# We need pretty visualizations to show relationship between RNA, Ribo, TE across individuals
# We need to do some GO Analysis

# FOR TE CALCULATION NEED TO NORMALIZE EVERYTHING TOGETHER OR 
# NEED UNIMODAL DISTRIBUTION OF COEFFICIENTS WHEN FITTING TREATMENT ONLY MODEL
te_design <- model.matrix(~treatment+treatment:sample_id_all)
te_fit <- lmFit(all_expr_elist, te_design)
te_fit2 <- eBayes(te_fit)
# Check that IDs are correct
#te_fit2$genes$ID <- as.character(CDS[isexpr ,1][ribo_id_in_rna])
# Here the significance testing for difference in TE within a given individual, ie
# Any of the coef of te_fit2 is problematic as testing against 0 is strange. 
# One idea might be to voom the entire table of RNA_Seq + Ribo_Seq for this
# This nested interaction formula should allow us to calculate within individual diffs. 


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


