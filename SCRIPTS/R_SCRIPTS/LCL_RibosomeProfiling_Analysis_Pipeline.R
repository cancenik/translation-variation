library("edgeR")
library("hexbin")
library("limma")
library("qtl")
library ("sva")
library("MASS")
library("kohonen")
library("pgirmess")
#source('~/project/kohonen_hexagonal/kohonnen_hex.R')

# Data Directory
data_dir <- '~/project/CORE_DATAFILES/'

## HGNC_to_ENSG
hgnc_to_ensg <- read.table(paste (data_dir, 'HGNCtoENSG.txt', sep=""), ,header=F,as.is=T,sep="|",fill=T)
ensg_hgnc <- cbind(grep("ENSG", unlist(strsplit(hgnc_to_ensg$V2, "[.]")), value=T), hgnc_to_ensg$V5)
enst_hgnc <- cbind(hgnc_to_ensg$V1, hgnc_to_ensg$V5)
colnames(ensg_hgnc) <- c("ENSG", "HGNC")
colnames(enst_hgnc) <- c("ENST", "HGNC")
substr(enst_hgnc[,1],1,1) <- ""

## Absolute Protein Amounts
protein_absolute_ibaq <- read.csv(paste (data_dir,'TableS8_Khan_etal.csv', sep=""))
protein_absolute_ibaq <- merge(protein_absolute_ibaq, ensg_hgnc)

## Linfenfg Proteins
linfeng_protein <- read.table(paste(data_dir, 'ldataMerged95x5953ENSGwithCovar.txt', sep=""))
linfeng_protein <- t(linfeng_protein[,-c(5955:5961)])
colnames(linfeng_protein) <- linfeng_protein[5954,]
linfeng_protein <- linfeng_protein[-5954,]
row.names(linfeng_protein) <- linfeng_protein[,1]

## RNA_SEQ COUNTS -- Drop 18508
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
all_rnaseq <- all_rnaseq[,-c(55, 77)]
all_rnaseq_counts <- DGEList(counts= all_rnaseq)

## Cufflinks Ratios
#cuff_ratios <- read.table(paste(data_dir, "Transcript_CufflinksRatios_wIDs", sep=""), header=T)
#cuff_ratios <- merge(cuff_ratios, ensg_hgnc, by.x="gene.identifier", by.y="ENSG")

## RIBOSEQ_COUNTS
# CDS SPECIES 
species <- read.table(paste(data_dir, "Reformatted_Species_Counts_All_Libraries.tsv", sep=""),header=T)
CDS_species <- split (species, species$REGION)[[2]]
CDS_species <- CDS_species[,grep("Counts", colnames(CDS_species))]

# Total Number of Reads
data <- read.table ("~/project/CORE_DATAFILES/Reformatted_Transcript_Counts_All_Libraries.tsv", header=T)
CDS <- split (data, data$REGION)[[2]]
UTR3 <- split (data, data$REGION)[[3]]
UTR5 <- split (data, data$REGION)[[4]]
CDS_Counts <- CDS[,grep("Counts", colnames(CDS))]
CDS_Coverage <- CDS[, grep("CoveredBases", colnames(CDS))]
CDS_Len <- CDS[,grep("Len", colnames(CDS))]
CDS_IDs <- CDS[,1]

#covariates <-  read.table ("~/project/CORE_DATAFILES/Sequenced_Ribosome_Profiling_Sample_Information_Batch_Effects.tsv", header=T)
######## DATA ANALYSIS ##################################
################### NORMALIZATION OF COUNTS AND INITIAL QC
## SPECIES TO COUNTS COMPARISON
species_sum <- rowSums(CDS_species)
cds_count_sum <- rowSums(CDS_Counts)
ratios <- cds_count_sum/species_sum
plot(log10(cds_count_sum+1), log10(species_sum+1), cex=0.2)
ratios_dataframe <- data.frame(ID=CDS_IDs, ReadCount=log10(cds_count_sum+1) , SpeciesCount=log10(species_sum+1) )

#### Perform loess regression between read_count to species_count
# Call outliers as 2*SE away from the fit
## This needs a lot of memory
#count_to_species <- predict(loess(ReadCount~SpeciesCount, data=ratios_dataframe, statistics="approximate", trace.hat="approximate"), se=T)
#c1 <- count_to_species$fit+10^2*count_to_species$s
#length(which(ratios_dataframe$ReadCount- c1 > 0))
# c2 <- count_to_species$fit-10^2*count_to_species$s
# plot(ratios_dataframe$ReadCount, ratios_dataframe$SpeciesCount, pch=19, cex=0.2)
#lines(ratios_dataframe$ReadCount,count_to_species$fit, col="red")
# lines(ratios_dataframe$ReadCount,count_to_species$fit+3*count_to_species$s, lty=2)
# lines(ratios_dataframe$ReadCount,count_to_species$fit-3*count_to_species$s, lty=2)

# RNASEQ NORMALIZATION
# Identify differences in PolyA to RiboZero and remove
rnaexpr <- rowSums(cpm(all_rnaseq_counts) > 1) >= 40
all_rnaseq_counts <- all_rnaseq_counts[rnaexpr,]
polyA_mean <- apply(all_rnaseq_counts[,grep("polyA", colnames(all_rnaseq_counts))],1, mean)
RZ_mean <- apply(all_rnaseq_counts[,grep("RiboZero", colnames(all_rnaseq_counts))],1, mean)
plot(log10(polyA_mean), log10(RZ_mean), cex=0.2, pch=19, ylab="log10(RiboZero_ReadCount)", xlab="log10(PolyA_ReadCount)", main="Comparing RNASeq Methods")
fit.rna <- lm(log10(RZ_mean)~log10(polyA_mean))
# rna_seq_normalized$ID[abs(stdres(fit.rna)) > 3]
outlier_colors <- rep("Black", length(log10(polyA_mean)))
outlier_colors[abs(stdres(fit.rna)) > 3] <- "Red"
plot(log10(polyA_mean),stdres(fit.rna) , pch=19, cex=0.2, xlab="Log10(PolyA_Reads)", ylab="Standardized Residuals", col=outlier_colors)
polyA_RZ_inconsistent <- (abs(stdres(fit.rna)) > 3)
all_rnaseq_counts <- all_rnaseq_counts[!polyA_RZ_inconsistent,]
row.names(all_rnaseq_counts) <- (geuvadis_CDS$ID[rnaexpr])[!polyA_RZ_inconsistent]
all_rnaseq_counts <- calcNormFactors (all_rnaseq_counts, method= "TMM")
sample_id_rna <- unlist(strsplit(colnames(all_rnaseq_counts), split= "_"))
sample_id_rna <- sample_id_rna[grep("GM", sample_id_rna)]
design_rnaseq <- model.matrix(~sample_id_rna)

# VOOM/TMM Normalization of Ribosome Profiling
cds_counts <- DGEList(counts=CDS_Counts)
isexpr <- rowSums(cpm(cds_counts) > 1) >= 36
cds_counts <- cds_counts[isexpr,]
row.names(cds_counts) <- CDS_IDs[isexpr]
cds_counts <- calcNormFactors (cds_counts, method= "TMM")
s1 <- unlist(strsplit(colnames(CDS_Counts), "_"))
sample_labels <- s1[grep('GM', s1)]
table(sample_labels)

############################################
### ALTERNATIVE AND PREFFERED APPROACH IS TO PROCESS ALL RNA_SEQ AND RIBO_SEQ
# RNA - all_rnaseq_counts
# Ribo - cds_counts
joint_counts <- merge(all_rnaseq_counts$counts, cds_counts$counts,  by="row.names")
s1 <- unlist(strsplit(colnames(joint_counts), "_"))
sample_labels_joint <- s1[grep('GM', s1)]
type <- c(rep("RNA", 84), rep("Ribo", 50) )
table(sample_labels_joint)
full_design <- model.matrix(~sample_labels_joint + type)
joint_count_ids <- joint_counts[,1]
joint_counts <- DGEList(counts= joint_counts[,-1])
joint_counts <- calcNormFactors (joint_counts, method= "TMM")
# Quantile normalization vs none, lowest correlation is 0.995 
#v3 <- voom(joint_counts, full_design, plot=T)
v3 <- voom(joint_counts, full_design, plot=T, normalize.method="quantile" )
# Joint voom vs sepearate lowest correlation .991 / .995 if no QN
# Ribosome profiling is similar
rnaseq_batch <- c (rep(1,18), rep(2, 26), rep(3,25), rep(4, 15) ) 
batch_removed_joint <- ComBat (v3$E[,1:84], batch=rnaseq_batch, mod=design_rnaseq)
#Update
v3$E[,1:84] <- batch_removed_joint

reps_ribo_joint <- duplicated(sample_labels_joint[85:134]) | duplicated(sample_labels_joint[85:134], fromLast = TRUE)
joint_dd_rep <- dist(t (v3$E[,85:134][,reps_ribo_joint] ) )
joint_ribo_rep <- hclust (joint_dd_rep)
plot(joint_ribo_rep)
reps_rna_joint <- duplicated(sample_labels_joint[1:84]) | duplicated(sample_labels_joint[1:84], fromLast = TRUE)
joint_dd_rna <- dist(t (v3$E[,1:84][,reps_rna_joint]))
joint_rna_hc <- hclust(joint_dd_rna)
plot(joint_rna_hc)
plotMDS(v3$E, labels=type)

# Extract Joint Replicated Subset
cor_coefs_rna_median <- apply(cor(v3$E[,1:84]),1, median)
cor_coefs_ribo_median <- apply(cor(v3$E[,85:134]),1, median)
hist(cor_coefs_ribo_median)
hist(cor_coefs_rna_median)
which(cor_coefs_rna_median < .9)
# GM18504_Rep2_Counts_Pickrell is weird. It's unlike all the other datasets
## UPDATE v3 by dropping GM18504_Rep2_Counts_Pickrell
v3 <- v3[,-51]
sample_labels_joint <- sample_labels_joint[-51]
type <- type[-51]
v3$design <- model.matrix(~sample_labels_joint + type)
row.names(v3) <- joint_count_ids

#save (v3, file="~/project/CORE_DATAFILES/Joint_Expression_RNA_Ribo_Voom")
#save (type, file="~/project/CORE_DATAFILES/Joint_Expression_RNA_Ribo_Classification_Factor")
# EXTRACT JOINT
joint_replication_ribo <- duplicated(sample_labels_joint[type=="Ribo"]) | duplicated(sample_labels_joint[type=="Ribo"], fromLast=TRUE)
joint_replication_rna <- duplicated(sample_labels_joint[type=="RNA"]) | duplicated(sample_labels_joint[type=="RNA"], fromLast=TRUE)
sample_labels_joint_wreps <- sample_labels_joint[c(joint_replication_rna, joint_replication_ribo)]
type_wreps <- type[c(joint_replication_rna, joint_replication_ribo)]

common_rna_ribo <- c((sample_labels_joint_wreps[type_wreps=="RNA"] %in% sample_labels_joint_wreps[type_wreps=="Ribo"]),
                     (sample_labels_joint_wreps[type_wreps=="Ribo"] %in% sample_labels_joint_wreps[type_wreps=="RNA"])) 

sample_labels_joint_common <- sample_labels_joint_wreps[common_rna_ribo]
type_common <- type_wreps[common_rna_ribo]
joint_expression_common <- v3[,c(joint_replication_rna, joint_replication_ribo)][,common_rna_ribo]
joint_expression_common$design <- model.matrix(~sample_labels_joint_common+type_common)
row.names(joint_expression_common) <- joint_count_ids

# # UPDATE v3 to remove GM18504_Rep2
# write.table(v3$E[,type=="RNA"], file =paste (data_dir, "TMM_VarianceMeanDetrended_QN_FullModel_RNA_Expression", sep=""), 
# row.names=joint_count_ids, sep="\t")
# write.table(v3$E[,type=="Ribo"], file =paste (data_dir, "TMM_VarianceMeanDetrended_QN_FullModel_RiboProfiling_Expression", sep=""), 
#             row.names=joint_count_ids, sep="\t")            


##### ANALYSIS REQUIRING JOINT RNA-RIBOSEQ WITH REPLICATES
### COMPARING VARIATION IN THE EXPRESSION VALUES
# Alternative much simpler approach is to use ratio of CVs
# It is unclear if CV makes sense as the relationship between Var Mean is different
# It might make more sense to multiply by the associated weights
weighted_sd <- function (y) { 
  return ( sd(y[,2]*y[,4]) )
}
# We can add a permutation scheme here to get p-values on these differences
# There doesn't seem to be any significant association between expression and the cv ratios as expected
# The easiest permutation scheme is to permute type_common to generate new RNA, Ribo Classification
# Run through existing code and compare the actual difference in CV witht the permutation p-value
# We might rely on nominal p-value threshold for selecting significant ones
# type_common will be modified making sure that each individual is preserved as in sample_labels_joint_common
# Loop over unique(sample_labels_joint_common)

rna_replicate_mean <- apply (joint_expression_common$E[,type_common=="RNA"], 1, function(x) {
aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="RNA"])), mean)  
} )
ribo_replicate_mean <- apply (joint_expression_common$E[,type_common=="Ribo"], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="Ribo"])), mean)  
} )
rna_replicate_weight_mean <- apply (joint_expression_common$weights[,type_common=="RNA"], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="RNA"])), mean)  
} )
ribo_replicate_weight_mean <- apply (joint_expression_common$weights[,type_common=="Ribo"], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="Ribo"])), mean)  
} )
rna_replicate_mean_weights <- mapply(cbind, rna_replicate_mean, rna_replicate_weight_mean, SIMPLIFY=F)
ribo_replicate_mean_weights <- mapply(cbind, ribo_replicate_mean, rna_replicate_weight_mean, SIMPLIFY=F)

rna_cv_between_individuals <- as.numeric(lapply(rna_replicate_mean_weights, weighted_sd))
ribo_cv_between_individuals <- as.numeric(lapply(ribo_replicate_mean_weights, weighted_sd))

# Calculate median CV of replicate CVs
rna_replicatecvs <- apply (joint_expression_common$E[,type_common=="RNA"]*joint_expression_common$weights[,type_common=="RNA"], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="RNA"])), sd)  
} )
ribo_replicatecvs <- apply (joint_expression_common$E[,type_common=="Ribo"]*joint_expression_common$weights[,type_common=="Ribo"], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="Ribo"])), sd)  
} )

rna_repcv_median <- as.numeric(lapply(rna_replicatecvs, function(z){median(z$x)}))
ribo_repcv_median <- as.numeric(lapply(ribo_replicatecvs, function(z){median(z$x)}))

rnacv <- (rna_cv_between_individuals/rna_repcv_median) 
ribocv <- (ribo_cv_between_individuals/ribo_repcv_median)

pdf(file = "~/Google_Drive/Manuscript Figures/RNA_Between_Individual_Variance.pdf", width=4, height=4)
#99% of the dataset is less than 8; so limit xlim to 0_to_8
# We changed number of breaks so that the number of breaks in x-axis is similar
hist(rnacv, 100, main= "RNA Expression", xlab = "Between/ Within Individual Coefficient of Variation", xlim=c(0,8))
dev.off()
pdf(file = "~/Google_Drive/Manuscript Figures/Ribo_Between_Individual_Variance.pdf", width=4, height=4)
hist(ribocv, 200 , main= "Ribosome Occupancy", xlab = "Between/ Within Individual Coefficient of Variation", xlim=c(0,8))
dev.off()

#low_sig_to_noise <- which( rnacv < 1 & ribocv < 1)
length(which(rnacv/ribocv > 2))
length(which(rnacv/ribocv < .5))
sorted_cvs <- sort(rnacv/ribocv , index.return=T)
#write.table(row.names(joint_expression_common)[sorted_cvs$ix], file = "~/project/CORE_DATAFILES/Sorted_InterIndividualCV_RNA_to_Ribo.txt", row.names=F)

p1 <- hist(rnacv/ribocv, 100 )
plot(p1, col=rgb(0,0,1,1/4),  tck=.02, xlab="RNA Expression to Ribosome Occupancy Ratio", main="Between/ Within Individual CV", ylim=c(0,500))
abline(v=1, lwd=3)


############################################################
# Differential Expression Analysis and Translation Efficiency
# Two predictors: Sample Label + Ribo vs RNA
# Questions of interest
# What are the genes with differential RNA expression across individuals
# What are the genes with differential Ribo expression across individuals
# What are the genes with differential Ribo-RNA expression across individuals
# What are the genes with differential Ribo-RNA expression in a given individual
# Related to this is the coefficient associate with Ribo-RNA for each gene in each individual
# The moderated F-statistic can be used as the measure of any difference
# ## EBayes also returns a moderated F-statistic, $F.p.value

# We can switch joint_expression_common, v3, and SVA 
all_expr_elist <- joint_expression_common
treatment <- relevel(as.factor(type_common),ref="RNA")
all_expr_elist$design <- model.matrix(~treatment + as.factor(sample_labels_joint_common))

ribo_expr_elist <- all_expr_elist[,treatment=="Ribo"]
ribo_expr_elist$design <- model.matrix(~0+as.factor(sample_labels_joint_common[treatment=="Ribo"]))
colnames(ribo_expr_elist$design) <- sort(unique(sample_labels_joint_common))

rna_expr_elist <- all_expr_elist[,treatment=="RNA"]
rna_expr_elist$design <- model.matrix(~0+as.factor(sample_labels_joint_common[treatment=="RNA"]))
colnames(rna_expr_elist$design) <- sort(unique(sample_labels_joint_common))

# GENERATE CONTRAST MATRIX
contrast_strings <- c()
for (i in 1:length(colnames(rna_expr_elist$design))) { 
contrast_strings[i] <- 
 paste ( 
paste (colnames(rna_expr_elist$design)[i], 
       paste (colnames(rna_expr_elist$design)[-i], collapse="+"), sep="- (")
, length(colnames(rna_expr_elist$design))-1 , sep=")/" )
}
# I hate to do this but will enumerate all strings here
contrast.matrix<- (makeContrasts(contrast_strings[1], 
                                 contrast_strings[2], 
                                 contrast_strings[3], 
                                 contrast_strings[4], 
                                 contrast_strings[5], 
                                 contrast_strings[6], 
                                 contrast_strings[7], 
                                 contrast_strings[8], 
                                 contrast_strings[9], 
                                 contrast_strings[10], 
                                 contrast_strings[11], 
                                 contrast_strings[12], 
                                 contrast_strings[13], 
                                 contrast_strings[14], 
                                levels=rna_expr_elist$design))

# MODEL FITTING
ribo_fit <- lmFit (ribo_expr_elist, design=ribo_expr_elist$design, weights=ribo_expr_elist$weights)
rna_fit <- lmFit (rna_expr_elist, design=rna_expr_elist$design, weights=rna_expr_elist$weights)

ribo_fit2 <- contrasts.fit(ribo_fit, contrast.matrix)
ribo_fit2 <- eBayes(ribo_fit2)
rna_fit2 <- contrasts.fit(rna_fit, contrast.matrix)
rna_fit2 <- eBayes(rna_fit2)

topTable(ribo_fit2)
topTable(rna_fit2)
results.ribo <- decideTests(ribo_fit2, p.value=0.01, lfc=log2(1.5))
results.rna <- decideTests(rna_fit2, p.value=0.01, lfc=log2(1.5))
as.numeric(apply(abs(results.ribo), 2, sum))
as.numeric(apply(abs(results.rna), 2, sum))
# We can also look at B-statistic for the propensity of a gene to be differentially expressed
# We need pretty visualizations to show relationship between RNA, Ribo, TE across individuals
# We need to do some GO Analysis

riborna_fit <- lmFit(all_expr_elist)
riborna_fit2 <- eBayes(riborna_fit)
# Almost everything is significant for treatRibo if no LFC threshold
# Using ordered FuncAssociate Pos_Ribo associates with ER, Mito, Extracellular space/organelle
# Neg-Ribo enriched in translation, ribosome, viral transcription, etc
apply(abs(decideTests(riborna_fit2, p.value=0.01, lfc=1)), 2, sum)
#write.table(topTable(riborna_fit2, coef=2, lfc=1, number=2326), file=paste(data_dir, 'Differential_Ribosome_Occupancy', sep=""), row.names=F ) 
# Here the significance testing for difference in TE within a given individual, ie
# Any of the coef of te_fit2 is problematic as testing against 0 is strange. 
# One idea might be to voom the entire table of RNA_Seq + Ribo_Seq for this
# This nested interaction formula should allow us to calculate within individual diffs. 
te_factor <- paste(treatment, sample_labels_joint_common, sep=".")
te_factor <- factor(te_factor, levels=unique(te_factor))
tedesign <- model.matrix(~0+ te_factor)
colnames(tedesign) <- levels(te_factor)
te_fit <- lmFit(all_expr_elist, design=tedesign)
# Contrast Matrix - generate matrix manually
# Interested in individual specific TE
# We are also interested in differential TE in given individual vs others
# Second one is similar to the contrast matrix in ribo/rna only models.
cont.matrix.te <- matrix(0,nrow=28, ncol=14, 
dimnames=list(Levels=levels(te_factor), Contrasts=unique(sample_labels_joint_common)))
sample_order <- c(unique(sample_labels_joint_common), unique(sample_labels_joint_common, fromLast=T))
for (j in 1:14) { 
  cont.matrix.te[which(sample_order == sample_order[j])[1],j] <- -1
  cont.matrix.te[which(sample_order == sample_order[j])[2],j] <- 1
}
te_fit4 <- contrasts.fit(te_fit, cont.matrix.te)
te_fit4 <- eBayes(te_fit4)
te.results <- decideTests (te_fit4, p.value=0.01, lfc=1)
as.numeric(apply(abs(te.results), 2, sum))


cont.matrix.diff.te <- matrix(0,nrow=28, ncol=14, dimnames=list(Levels=levels(te_factor), Contrasts=unique(sample_labels_joint_common)))
sample_order <- c(unique(sample_labels_joint_common), unique(sample_labels_joint_common, fromLast=T))
for (j in 1:14) {
  c1 <- which(sample_order == sample_order[j])
  c2 <- seq(1,28,1)[-c1]
  cont.matrix.diff.te[c1[1],j] <- -1
  cont.matrix.diff.te[c1[2],j] <- 1
  cont.matrix.diff.te[c2[which(c2<=14)],j] <- 1/13
  cont.matrix.diff.te[c2[which(c2>14)],j] <- -1/13
}
te_fit3 <- contrasts.fit(te_fit, cont.matrix.diff.te)
te_fit3 <- eBayes(te_fit3)
te.diff.results <- decideTests (te_fit3, p.value=0.01, lfc=log2(1.5))
# Here the weird behavior of GM18504 RNA-Seq causes a large number of differences in TE
as.numeric(apply(abs(te.diff.results), 2, sum))

#### We will visualize the results using Radial Sets
#### We will need the following information
#### An Attributes table
#### GENE_ID RNA_EXP RIBO_EXP TE PROT LEVEL THIS COULD BE IND SEPARATE OR MEANS
#### THREE TABLES OF DIFFERENTIAL EXPRESSION
#### RIBO TABLE - EXAMPLE ONE FOR UP ONE FOR DOWN
#### IND_ID GENE_LIST
# Results of differentials are in te.diff.results, results.ribo, results.rna
# This best achieved in perl so simply output these tables
#write.table(results.ribo, file =paste (data_dir, "Differential_Ribo_Occupancy", sep=""))
#write.table(results.rna, file =paste (data_dir, "Differential_RNA_Occupancy", sep=""))
#write.table(te.diff.results, file =paste (data_dir, "Differential_TE", sep=""))

# In general genes exhibit great variation in their translation efficiency
# However, across individuals the translation efficiency of genes is much less variable than RNA
# For absolute levels, the translation efficiency component captured by ribosome profiling is important to predict protein levels
# However, across individual differences in TE have lower correlation between across individual protein measurements perhaps due to less variance


## ANALYSIS ON ALL DATA INCLUDING RNA-RIBO IRRESPECTIVE OF REPLICATION
# Comparing across Individual Correlation to Linfeng's proteomics
# ensg_hgnc
linfeng_common <- colnames(linfeng_protein) %in% unique(sample_labels_joint)
linfeng_protein_common <- linfeng_protein[,linfeng_common]
linfeng_protein_common <- merge(linfeng_protein_common, ensg_hgnc, by.x="row.names", by.y="ENSG")
# Subset linfeng protein to only no NAs. If we want we can also include samples with missing values sqrt(#individuals)?
# For correlation we can use use="pairwise.complete.obs"
number_NAs <- 1
linfeng_protein_na <- linfeng_protein_common[apply(is.na(linfeng_protein_common), 1, sum) < number_NAs, ]
linfeng_protein_ribo_rna <- merge (v3$E, linfeng_protein_na, by.x="row.names", by.y="HGNC")
row.names(linfeng_protein_ribo_rna) <- linfeng_protein_ribo_rna[,1]
linfeng_protein_ribo_rna <- linfeng_protein_ribo_rna[,-c(1,135)]
type_prot <- c(type, rep("Prot", dim(linfeng_protein_ribo_rna)[2] - length(type)))
sample_labels_joint_prot <- c(sample_labels_joint, colnames(linfeng_protein_ribo_rna)[134:161])

# Calculate across individual correlation between protein-rna-ribo
rna_replicate_mean_prot <- apply (linfeng_protein_ribo_rna[,type_prot=="RNA"], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_joint_prot[type_prot=="RNA"])), mean)  
} )
ribo_replicate_mean_prot <- apply (linfeng_protein_ribo_rna[,type_prot=="Ribo"], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_joint_prot[type_prot=="Ribo"])), mean)  
} )

rna_samples <- match(as.character(rna_replicate_mean_prot[[1]]$Group.1), sample_labels_joint_prot[type_prot=="Prot"],)
rna_samples <- rna_samples[!is.na(rna_samples)]
rna_in_prot <- as.character(rna_replicate_mean_prot[[1]]$Group.1) %in% sample_labels_joint_prot[type_prot=="Prot"]
across_ind_rna_correlation <- c()
across_ind_rna_correlation_pval <- c()

for (i in 1:length(rna_replicate_mean_prot)) { 
#  cor1 <-  cor.test(rna_replicate_mean_prot[[i]]$x[rna_in_prot], as.numeric(as.matrix(linfeng_protein_ribo_rna[i,type_prot=="Prot"]))[rna_samples], use="pairwise.complete.obs")
  cor1 <-  cor.test(rna_replicate_mean_prot[[i]]$x[rna_in_prot], as.numeric(as.matrix(linfeng_protein_ribo_rna[i,type_prot=="Prot"]))[rna_samples], use="pairwise.complete.obs", method="spearman")
  across_ind_rna_correlation <- c(across_ind_rna_correlation, cor1$estimate)
  across_ind_rna_correlation_pval <- c(across_ind_rna_correlation_pval, cor1$p.value)
}  

ribo_samples <- match(as.character(ribo_replicate_mean_prot[[1]]$Group.1), sample_labels_joint_prot[type_prot=="Prot"],)
ribo_samples <- ribo_samples[!is.na(ribo_samples)]
ribo_in_prot <- as.character(ribo_replicate_mean_prot[[1]]$Group.1) %in% sample_labels_joint_prot[type_prot=="Prot"]
across_ind_ribo_correlation <- c()
across_ind_ribo_correlation_pval <- c()

for (i in 1:length(ribo_replicate_mean_prot)) { 
  # Get two numeric matching vectors
#  cor1 <-cor.test(ribo_replicate_mean_prot[[i]]$x[ribo_in_prot], as.numeric(as.matrix(linfeng_protein_ribo_rna[i,type_prot=="Prot"]))[ribo_samples], use="pairwise.complete.obs")
  cor1 <-cor.test(ribo_replicate_mean_prot[[i]]$x[ribo_in_prot], as.numeric(as.matrix(linfeng_protein_ribo_rna[i,type_prot=="Prot"]))[ribo_samples], use="pairwise.complete.obs", method="spearman")  
  across_ind_ribo_correlation <- c(across_ind_ribo_correlation, cor1$estimate)
  across_ind_ribo_correlation_pval <- c(across_ind_ribo_correlation_pval, cor1$p.value)
} 

length(across_ind_ribo_correlation)
median(across_ind_ribo_correlation)
median(across_ind_rna_correlation)
color_by_pval <- rep(0, length(ribo_replicate_mean_prot))
# FDR ~ 25%
pval_cutoff <- 0.0001
color_by_pval[across_ind_ribo_correlation_pval < pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff] <- 1
color_by_pval[across_ind_ribo_correlation_pval< pval_cutoff & across_ind_rna_correlation_pval >= pval_cutoff] <- 2
color_by_pval[across_ind_ribo_correlation_pval>=pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff] <- 3
plot(across_ind_ribo_correlation, across_ind_rna_correlation, pch=19, cex=.65, tck=.02, col=c("Black", "Red", "Blue", "Gold")[as.factor(color_by_pval)], 
     xlab="Between Individual Ribosome Occupancy-Protein Expression Correlation", ylab = "Between Individual RNA-Protein Expression Correlation")
cor.test(across_ind_ribo_correlation, across_ind_rna_correlation)
#### ADD FISHER'S TEST FOR ENRICHMENT -- Huge Enrichment for Both correlating significantly
fmat <- matrix(nrow=2, ncol=2)
fmat[1,1] = length(which(across_ind_ribo_correlation_pval < pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff))
fmat[1,2] =  length(which(across_ind_ribo_correlation_pval < pval_cutoff & across_ind_rna_correlation_pval > pval_cutoff))
fmat[2,1] = length(which(across_ind_ribo_correlation_pval > pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff))
fmat[2,2] = length(which(across_ind_ribo_correlation_pval > pval_cutoff & across_ind_rna_correlation_pval > pval_cutoff))
fisher.test(fmat)

# Remove one column which is not shared
ribo_replicate_mean_rna <- lapply(ribo_replicate_mean_prot, function(x){x <- x[-24,]})
c2 <- mapply (cbind, rna_replicate_mean_prot, ribo_replicate_mean_rna, SIMPLIFY=F)
across_ind_rna_ribo <- as.numeric(lapply(c2, function(x){ cor(x[,2], x[,4],method="spearman") }))
# Histograms of Across Ind Ribo-Prot, RNA-Prot and Ribo-RNA correlations
p1 <- hist(across_ind_ribo_correlation,40)
p2 <- hist(across_ind_rna_correlation,40)
p3 <- hist(across_ind_rna_ribo, 40)
pdf(file = "~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Across_Individual_Correlations.pdf", width=9, height=6.5)
plot(p1, col=rgb(0,0,1,1/4), xlim=c(-1,1), xlab="Spearman Correlation Coefficient", main="Correlation Coefficient Distribution")
plot(p2, col=rgb(1,0,0,1/4), xlim=c(-1,1), add=T)
plot(p3, col=rgb(0,1,0,1/4), xlim=c(-1,1), add=T)
legend(.4,170,c("RNA Expression-\nProtein Expression", "Ribosome Occupancy-\nProtein Expression", "RNA Expression-\nRibosome Occupancy"), 
       yjust =0.5, x.intersp=0.2, y.intersp=1.5,bty="n", border="white", fill=c(rgb(1,0,0,1/4), rgb(0,0,1,1/4), rgb(0,1,0,1/4)), cex=.9)
dev.off()
quantile(across_ind_rna_ribo)
quantile(across_ind_ribo_correlation)
quantile(across_ind_rna_correlation)

# We can do these sorted so we can run ordered analysis
# Another Idea is to merge all RNA significant and all Ribo Significant
# RNA significant chemical response; extracellular region
# Ribo Response Regulation of phosphorylation + DAVID:Post-translational modifications
# Extracted p < 1e-4 or p<1e-3
#write.table( row.names(linfeng_protein_ribo_rna)[across_ind_ribo_correlation_pval < pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff], file=paste (data_dir,'Ribo_RNA_Prot_Cor_IDs' ,sep=""), row.names=F ) 
#write.table(row.names(linfeng_protein_ribo_rna)[across_ind_ribo_correlation_pval< pval_cutoff & across_ind_rna_correlation_pval >= pval_cutoff], file=paste (data_dir,'Ribo_Prot_Cor_IDs',sep=""), row.names=F ) 
#write.table(row.names(linfeng_protein_ribo_rna), file= paste (data_dir,'Prot_RNA_Ribo_Common_IDs',sep=""), row.names=F)
#write.table(row.names(linfeng_protein_ribo_rna)[across_ind_ribo_correlation_pval>=pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff], file=paste (data_dir,'RNA_Prot_Cor_IDs',sep=""), row.names=F)
## COMPARISION TO CHRISTINE'S PROTEOMICS
# Christine Quantification correlation to SILAC proteomics is pretty low. May be better not use this
#### Compare absolute levels of protein with rna and ribo -- Overall correlation is better with ribosome profiling
gm12878_prot <- read.csv('~/project/CORE_DATAFILES/GM12878_B0_FDR5_140120_shortforCan.csv', stringsAsFactors=F)
gm12878_prot <- merge(gm12878_prot, ensg_hgnc, by.x="ID.1", by.y="ENSG")
gm12878_prot = gm12878_prot[,c(5,14)]
all_prot_data = merge(gm12878_prot, protein_absolute_ibaq, by="HGNC")
all_prot_data = all_prot_data[-which(all_prot_data$USE == 0), ]

plot(log10(as.numeric(all_prot_data$USE)), log10(all_prot_data$ibaq.human), pch = 19, cex=.65, xlab="GM12878_Label_Free", ylab="SILAC", main="log10 iBAQ")
abline(lm(log10(all_prot_data$ibaq.human) ~ log10(as.numeric(all_prot_data$USE))), lwd = 2)
lines(lowess(log10(all_prot_data$ibaq.human) ~ log10(as.numeric(all_prot_data$USE))), lwd=2)
cor.test(all_prot_data$ibaq.human, as.numeric(all_prot_data$USE), method="spearman")

fit.prot <- lm(log10(all_prot_data$ibaq.human) ~ log10(as.numeric(all_prot_data$USE)))

outlier_colors <- rep("Black",times = length(log10(all_prot_data$ibaq.human) ))
outlier_colors[abs(stdres(fit.prot)) > 1] <- "Red"
plot( log10(as.numeric(all_prot_data$USE)),stdres(fit.prot) , pch=19, cex=0.2, ylab="Standardized Residuals", col=outlier_colors)
prot_inconsistent <- (abs(stdres(fit.prot)) > 1)
all_prot_data = all_prot_data[!prot_inconsistent, ]

gm12878_ribo <- v3$E[,c(85,86)]
gm12878_ribo_mean <- apply(gm12878_ribo, 1, median)
gm12878_ribo_mean <- data.frame(HGNC=row.names(v3), gm12878_ribo_mean)
gm12878_rna <- v3$E[,c(1,2,3,21,22)]
gm12878_rna_mean  <- apply(gm12878_rna, 1, median)
gm12878_rna_mean <- data.frame(HGNC=row.names(v3), gm12878_rna_mean)

ribo_rna_12878 <- merge(gm12878_ribo_mean, gm12878_rna_mean, by="HGNC")
ribo_rna_prot_12878 <- merge(ribo_rna_12878, gm12878_prot, by="HGNC")
ribo_rna_prot_consistent_12878 <- merge(ribo_rna_12878, all_prot_data, by="HGNC")
ribo_rna_prot_12878 <- ribo_rna_prot_12878[as.numeric(ribo_rna_prot_12878$USE) > 1000 & ribo_rna_prot_12878$gm12878_ribo_mean > 5, ] 

plot(ribo_rna_prot_12878$gm12878_ribo_mean, ribo_rna_prot_12878$gm12878_rna_mean, cex=0.2, pch=19, xlim=c(5,15), ylim=c(3,15))
plot(ribo_rna_12878$gm12878_ribo_mean, ribo_rna_12878$gm12878_rna_mean, cex=0.2, pch=19, xlim=c(5,15), ylim=c(3,15))
plot(ribo_rna_prot_12878$gm12878_ribo_mean,  log10(as.numeric(ribo_rna_prot_12878$USE)), pch=19, cex=.2, xlim=c(4, 15))
plot(ribo_rna_prot_12878$gm12878_rna_mean,  log10(as.numeric(ribo_rna_prot_12878$USE)), pch=19, cex=.2, xlim=c(4, 15))

#Spearman - Grand Mean is BEST
cor.test(ribo_rna_prot_consistent_12878$gm12878_rna_mean, log10(as.numeric(ribo_rna_prot_consistent_12878$USE)))
cor.test(ribo_rna_prot_consistent_12878$gm12878_ribo_mean, log10(as.numeric(ribo_rna_prot_consistent_12878$USE)))

cor.test(ribo_rna_prot_12878$gm12878_ribo_mean, log10(as.numeric(ribo_rna_prot_12878$USE)), method="spearman")
cor.test(ribo_rna_prot_12878$gm12878_rna_mean, log10(as.numeric(ribo_rna_prot_12878$USE)), method="spearman")

cor.test(ribo_rna_prot_12878$gm12878_ribo_mean, log10(as.numeric(ribo_rna_prot_12878$USE)))
cor.test(ribo_rna_prot_12878$gm12878_rna_mean, log10(as.numeric(ribo_rna_prot_12878$USE)))

#### COMPARISION WITH ACROSS GENE QUANTIFICATION FROM SILAC
#CDS_Len <- data.frame(LENGTH=CDS_Len[,1])
#row.names(CDS_Len) <- CDS_IDs
#CDS_Len <- CDS_Len[row.names(CDS_Len) %in% joint_count_ids,]

grand_mean_rna <- apply (v3$E[,type=="RNA"], 1, median)
grand_mean_rna  <- data.frame(HGNC=joint_count_ids, grand_mean_rna)
grand_mean_ribo <- apply(v3$E[,type=="Ribo"], 1, median)
grand_mean_ribo <- data.frame (HGNC=joint_count_ids, grand_mean_ribo)

# Calculated TE correlates weakly with Protein amount but better than ratio
# It might be again due to polysome profile shape
# The amount of protein TE correlation is likely cell type specific
grand_mean_te <- apply(te_fit4$coefficients, 1, median)
mean_te_df <- data.frame(HGNC=row.names(te_fit4), grand_mean_te)
m1 <- merge(mean_te_df, gm12878_prot, by="HGNC")
c1 <- m1$grand_mean_te > -1
plot(m1$grand_mean_te[c1], log10(as.numeric(m1$USE))[c1], pch=19, cex=.2)
cor.test(m1$grand_mean_te[c1], log10(as.numeric(m1$USE)+1)[c1], method="spearman")

ribo_rna_te <- merge(merge(grand_mean_ribo, grand_mean_rna, by="HGNC"), mean_te_df, by="HGNC")
ribo_rna_te_prot <- merge (ribo_rna_te, protein_absolute_ibaq, by="HGNC", all.x=T)
ribo_rna_te_prot <- ribo_rna_te_prot[,-c(5,6,8,9)]
#write.table(ribo_rna_te_prot, file=paste(data_dir, "Radial_Sets_Attribute_Table_4Levels_Gene_Expression.txt", sep=""), row.names=F)
# We might want to standardize the measurements
ribo_rna_te_prot$ibaq.human <- log10(ribo_rna_te_prot$ibaq.human)
all_cors = cor(ribo_rna_te_prot[,-1], use="complete.obs", method="spearman")
plot(ribo_rna_te_prot$grand_mean_te, ribo_rna_te_prot$ibaq.human, pch= 19, cex =.4, tck = .02, xlim = c(-3,3), xlab="Median Translation Efficiency", ylab="log10 iBAQ protein expression")
text(par("usr")[2] - 0.75, par("usr")[4] - 0.75, labels= paste("Spearman rho", signif(all_cors[3,4], 2) , sep=" = "))
ribo_rna_te_prot[,-1] <- scale(ribo_rna_te_prot[,-1], scale=F)
row.names(ribo_rna_te_prot) <- ribo_rna_te_prot[,1]
ribo_rna_te_prot <- as.matrix(ribo_rna_te_prot[,2:5])

CDS_Lens <- data.frame(HGNC=CDS_IDs, CDS_Len[,1])
merge_ribo_prot <- merge(grand_mean_ribo,protein_absolute_ibaq, by="HGNC" )
merge_ribo_rna_prot <- merge (merge_ribo_prot, grand_mean_rna, by="HGNC")
merge_ribo_rna_prot_len <- merge(merge_ribo_rna_prot, CDS_Lens, by="HGNC")
dim(merge_ribo_rna_prot)
# Length normalization increases spearman correlation but lowers Pearson correlation
# In all possible comparisons riboseq is better correlated
rna_cor <- cor.test(merge_ribo_rna_prot$grand_mean_rna, log10(merge_ribo_rna_prot$ibaq.human))
ribo_cor <- cor.test(merge_ribo_rna_prot$grand_mean_ribo, log10(merge_ribo_rna_prot$ibaq.human))
ribo_rna <- cor.test(merge_ribo_rna_prot$grand_mean_ribo, merge_ribo_rna_prot$grand_mean_rna)

rna_cor <- cor.test(merge_ribo_rna_prot$grand_mean_rna, log10(merge_ribo_rna_prot$ibaq.human), method="spearman")
ribo_cor <- cor.test(merge_ribo_rna_prot$grand_mean_ribo, log10(merge_ribo_rna_prot$ibaq.human), method="spearman")

plot(merge_ribo_rna_prot$grand_mean_rna, log10(merge_ribo_rna_prot$ibaq.human), xlab="RNA Expression", ylab="log10 iBAQ protein expression", pch=19, cex=.4)
text(par("usr")[2]-0.5, par("usr")[4]-0.5, labels=paste ("R^2", round(rna_cor$estimate^2,2) , sep="=") , adj=c(1,1), cex=2)
plot(merge_ribo_rna_prot$grand_mean_ribo, log10(merge_ribo_rna_prot$ibaq.human), xlab="Ribosome Profiling Expression", ylab="log10 iBAQ protein expression", pch=19, cex=.4)
text(par("usr")[2]-0.5, par("usr")[4]-0.5, labels=paste ("R^2 ", round(ribo_cor$estimate^2,2) , sep="=") , adj=c(1,1), cex=2)
plot(merge_ribo_rna_prot$grand_mean_ribo, merge_ribo_rna_prot$grand_mean_rna, xlab="Ribosome Profiling Expression", ylab="RNA Expression", pch=19, cex=.4)
text(par("usr")[2]-0.5, par("usr")[4]-0.5, labels=paste ("R^2 ", round(ribo_rna$estimate^2,2) , sep="=") , adj=c(1,1), cex=2)
#

## APPLY SOMs to different versions of data. 
# APPLY BDK to ribo_rna_te_prot
# We can do a bidirectional kohonen predicting protein levels and plot the cor per cell with either or all
# We can cluster cells by which parameter correlates best with protein levels in each cell. 
total_cells_bdk <- floor(sqrt(dim(ribo_rna_te_prot)[2]/2) * sqrt (dim(ribo_rna_te_prot)[1]))
if (floor(sqrt(total_cells_bdk/1.3333)) %% 2 == 0) { 
  ydim = floor(sqrt(total_cells_bdk/1.3333))
} else { 
  ydim = floor(sqrt(total_cells_bdk/1.3333)) + 1
}
xdim = floor(total_cells_bdk/ydim + 0.5)
abs.som.data.noNA <- ribo_rna_te_prot[!apply(is.na(ribo_rna_te_prot), 1, any),]

# Run the SOM, 1000 times and keep track of the distances pick the one with the min 75% distance
my_75th_distance <- 1
for (i in 1:100) { 
  set.seed(i*3433)
  absolute.som <- bdk(abs.som.data.noNA[,1:3], Y= abs.som.data.noNA[,4] , toroidal=T, xweight=.8, contin=T, grid=somgrid(xdim, ydim, "hexagonal"))
  if (quantile(absolute.som$distances)[4] < my_75th_distance) { 
    my_seed <- i*3433
    my_75th_distance <- quantile(absolute.som$distances)[4]
  }
}
# best_seed is 130454
set.seed(130454)
absolute.som <- bdk(abs.som.data.noNA[,1:3], Y= abs.som.data.noNA[,4] , toroidal=T, xweight=.8, contin=T, grid=somgrid(xdim, ydim, "hexagonal"))
plot(absolute.som)
plot(absolute.som, type="quality")
plot(absolute.som, type="count")
plot(absolute.som, type="changes")

# One version is te, rna, ribo logFCs with Linfeng's proteomics. This is a 4x14x9000 matrix

# superSOM Data Structure is a list of matrices including Linfeng Proteomics with NAs -> This will use individuals in the plot.
# We can also do a general one with just absolute ribo,rna, te, and absolute SILAC amounts (Here SILAC can be coded as a classification parameter for BDF)
# Following Xie, Boyle, et al. The grid is hexagonal, toroid
# Total number of cells, sqrt(DATA_TYPE/2) x sqrt (genes x individuals )
linfeng_prot_common_with_te <-  colnames(linfeng_protein) %in% colnames(te_fit3$coefficients)
linfeng_te_columns <- linfeng_protein[,linfeng_prot_common_with_te]
class(linfeng_te_columns) <- "numeric"
linfeng_te_match <- cbind ( linfeng_te_columns[,1], rep(NA, dim(linfeng_te_columns)[1]), linfeng_te_columns[,2:5], rep(NA, dim(linfeng_te_columns)[1]), linfeng_te_columns[,6:12])
colnames(linfeng_te_match) <- sort(colnames(te_fit3$coefficients))
linfeng_te_match <- merge(linfeng_te_match, ensg_hgnc, by.x="row.names", by.y="ENSG")
linfeng_te_match <- merge (data.frame(HGNC=joint_count_ids), linfeng_te_match, by="HGNC", all.x=T)
row.names(linfeng_te_match) <- linfeng_te_match[,1]
linfeng_te_match <- linfeng_te_match[,-c(1,2)]
### We might want to scale the matrix before feeding it to SOM 
quantile(ribo_fit2$coefficients )
quantile(rna_fit2$coefficients )
quantile(te_fit3$coefficients )
quantile(linfeng_te_match, na.rm=T)

colnames(ribo_fit2$coefficients) <- sort(colnames(te_fit3$coefficients))
colnames(rna_fit2$coefficients) <- sort(colnames(te_fit3$coefficients))
som.data = list ( ribo= ribo_fit2$coefficients , rna = rna_fit2$coefficients, te = te_fit3$coefficients[,sort(colnames(te_fit3$coefficients), index.return=T)$ix])
som.data.prot = list ( ribo= ribo_fit2$coefficients , rna = rna_fit2$coefficients, te = te_fit3$coefficients[,sort(colnames(te_fit3$coefficients), index.return=T)$ix], prot=as.matrix(linfeng_te_match))

total_cells <- floor(sqrt(length(som.data)/2) * sqrt (dim(ribo_fit2$coefficients)[1] * dim(ribo_fit2$coefficients)[2]))
if (floor(sqrt(total_cells/1.3333)) %% 2 == 0) { 
  ydim = floor(sqrt(total_cells/1.3333))
} else { 
  ydim = floor(sqrt(total_cells/1.3333)) + 1
}
xdim = floor(total_cells/ydim + 0.5)
# We can play with weights Increased weight to RNA, Ribo compared to TE - Change weights to increase quality of SOM
som.exp = supersom(data =som.data, grid=somgrid(xdim, ydim, "hexagonal"), toroidal=T, contin=T)
som.exp.prot = supersom(data =som.data.prot, grid=somgrid(xdim, ydim, "hexagonal"), toroidal=T, contin=T, maxNA.fraction=9/14, weights=c(.4,.4,.15,.05))

plot(som.exp, type="codes")
plot(som.exp, type="quality")
plot(som.exp, type="mapping", pch=19, cex=.3)
plot(som.exp, type="changes")
plot(som.exp, type="counts")
## Hexagonal plotting 
# som.exp$unit.classif has the info about where each gene went
# Create Matrix of the quantity of interest and pass this directly to plotCplane and remove componentPlaneMatrix function
# We also need some visually appealing colors
## Write function to apply function to each cell of SOM.
## We will plot these once we have the good tools for plotting
# som.exp$unit.classif keeps track of the cell each object belongs to
# Left bottom is 1, the row above that is 1+xdim 
# Need vector of length xdim * ydim

# Calculate cor.estimate for each level grouped by unit.classif
# X is passed as a dataframe to function
#problem x is dataframe need to apply cor to each row
cor_between_ind <- function (x) { 
 cors <- apply(x, 1, function(y){cor(y[1:14], y[15:28], use="pairwise.complete.obs", method="spearman")})
 return (median(cors, na.rm=T))
}

# It might make more sense to make these plots with the som that includes the protein levels
ribo_prot_cor_in_som <- by(data.frame(cbind(som.exp$data$ribo,linfeng_te_match)), som.exp$unit.classif, FUN = cor_between_ind )
rna_prot_cor_in_som <- by(data.frame(cbind(som.exp$data$rna,linfeng_te_match)), som.exp$unit.classif, FUN = cor_between_ind )
te_prot_cor_in_som <- by(data.frame(cbind(som.exp$data$te,linfeng_te_match)), som.exp$unit.classif, FUN = cor_between_ind )

median (ribo_prot_cor_in_som, na.rm=T)
median(across_ind_ribo_correlation)
median (rna_prot_cor_in_som, na.rm=T)
median(across_ind_rna_correlation)
median(te_prot_cor_in_som)

plot.kohonen(som.exp, type="property", property=ribo_prot_cor_in_som, main= "Between Individual Ribosome Occupancy Protein Level Correlation")
plot.kohonen(som.exp, type="property", property=rna_prot_cor_in_som, main = "Between Individual RNA Occupancy Protein Level Correlation")
plot.kohonen(som.exp, type="property", property=te_prot_cor_in_som, main = "Between Individual Translation Efficiency Protein Level Correlation")

plot.kohonen(som.exp, type="counts" )

# absolute.som$data -> Numeric Matrix
# abs.som.data.noNA -> Numeric Matrix, 4th column is ibaq.human
ribo_prot_cor_across_genes_som <- by (data.frame(abs.som.data.noNA[,c(1,4)]), absolute.som$unit.classif, FUN = function(x){cor(x)[1,2]} )
rna_prot_cor_across_genes_som <- by (data.frame(abs.som.data.noNA[,c(2,4)]), absolute.som$unit.classif, FUN = function(x){cor(x)[1,2]} )
te_prot_cor_across_genes_som <- by (data.frame(abs.som.data.noNA[,c(3,4)]), absolute.som$unit.classif, FUN = function(x){cor(x)[1,2]} )
prot_mean <-  by (data.frame(abs.som.data.noNA[,4]), absolute.som$unit.classif, FUN = colMeans)

plotCplane(absolute.som, variable=ribo_prot_cor_across_genes_som)
plotCplane(absolute.som, variable=rna_prot_cor_across_genes_som)
plotCplane(absolute.som, variable=te_prot_cor_across_genes_som)

plot.kohonen(absolute.som, property=ribo_prot_cor_across_genes_som, type="property", main ="Ribosome Occupancy Protein Correlation")
plot.kohonen(absolute.som, property=rna_prot_cor_across_genes_som, type="property", main = "RNA Expression Protein Correlation")
plot.kohonen(absolute.som, property=te_prot_cor_across_genes_som, type="property", main="Translation Efficiency Protein Correlation")
plot.kohonen(absolute.som, property=prot_mean, type = "property")
#abs.som.data.noNA
# absolute.som

#### ANALYSES BASED ON JUST RIBOSOME PROFILING
### KOZAK SEQUENCE ANALYSIS
# Each transcript will contribute a delta Kozak and delta expression and type of variant
# We dropped 18508 from the analysis but it can still be part of the individuals. 
# There are a bunch of ATG variants; our PWM scores for these are identical for variant and ref 
# These are also filtered out

# Do not test ones with MAF < XX 
kozak_scores <- read.table('~/project/CORE_DATAFILES/Kozak_Reference_Sequence_Scores.txt')
kozak_score_variants <- read.table('~/project/CORE_DATAFILES/Kozak_Variant_Sequence_Scores.txt',
stringsAsFactors=FALSE, fill=T, col.names=paste ("V", seq(1:60), sep=""))
kozak_score_variants_hapmap <- read.table('~/project/CORE_DATAFILES/Kozak_Variant_Sequence_Scores_HapMap.txt', 
stringsAsFactors=F, fill=T, col.names=paste ("V", seq(1:10), sep=""))
all_kozak_score_variants <- merge(kozak_score_variants, kozak_score_variants_hapmap, all=T, by=c("V1", "V2", "V3", "V4") )
all_kozak_score_variants[is.na(all_kozak_score_variants)] <- ""
for ( i in 1:dim(all_kozak_score_variants)[2]) { 
  for (j in 1:dim(all_kozak_score_variants)[1]) { 
    if (all_kozak_score_variants[j,i] == "NA18508" ) { 
      all_kozak_score_variants[j,i] = ""
    }
  }
}
gm18508_unique_variants = 
  which (apply(all_kozak_score_variants, 1, function(x) { length(which(x==""))} ) ==62)
all_kozak_score_variants <- all_kozak_score_variants[-gm18508_unique_variants, ]
# Filter out variants that change the ATG
atg_variants <- grep ("6|7|8", all_kozak_score_variants$V2)
all_kozak_score_variants <- all_kozak_score_variants[-atg_variants,]
kozak_var_ind <- all_kozak_score_variants[,-c(2,3,4)]
kozak_score_variants<- all_kozak_score_variants[,c(1,4)]
kozak_merge <- merge(kozak_score_variants, kozak_scores, by="V1")

multi <- duplicated(kozak_merge[,1]) | duplicated(kozak_merge[,1], fromLast=T)
number_alleles <- apply (kozak_var_ind, 1, function(x){length(grep('NA', x))})
number_alleles <- number_alleles[!multi]
MAF <- floor(0.05 * 60)
hist(number_alleles, 60, xlab="Number of Alleles", main="")

# Comaparing the strength of the wild-type vs the variant Kozak strength
# Let's not worry about allele_frequency. Very few cases with alternative allele very high frequency
p1 <- hist(kozak_merge$V2,20)
#p2 <- hist(kozak_merge$V4[number_alleles <= 30], 50)
p2 <- hist(kozak_merge$V4,20)
plot(p1, col=rgb(0,0,1,1/4), xlim=c(-16,-6), xlab="Kozak Score", main="Distribution of Kozak Scores")
plot(p2, col=rgb(1,0,0,1/4), xlim=c(-16,-6), add=T)
boxplot(kozak_merge$V4, kozak_merge$V2, col=c(rgb(1,0,0,1/4), rgb(0,0,1,1/4)), notch=T, range=1, xlab=NULL)
wilcox.test(kozak_merge$V4, kozak_merge$V2)
pdf(file="~/Google_Drive/Manuscript Figures/Kozak_Analysis/Difference_in_Kozak_Scores.pdf", width=3, height=4, bg= "transparent")
boxplot(kozak_merge$V2 -kozak_merge$V4, notch=T)
abline(h = 0)
dev.off()
rna_only <- v3[,type=="RNA"]
sample_labels_rna <- unlist(strsplit(colnames(rna_only), split= "_"))
sample_labels_rna <- sample_labels_rna[grep("GM", sample_labels_rna)]
ribo_only <- v3[,type=="Ribo"]
sample_labels_ribo <- unlist(strsplit(colnames(ribo_only), split= "_"))
sample_labels_ribo <- sample_labels_ribo[grep("GM", sample_labels_ribo)]

#Kozak multi -- There is one position where all 60 is alternate
multi_diff = kozak_merge$V2[multi] - kozak_merge$V4[multi]
multi_ind = kozak_var_ind[multi,]
multi_alleles <- apply(multi_ind, 1, function(x){length(grep('NA', x))})
total_alleles <- by(multi_alleles, as.factor(multi_ind[,1]), sum)  

# One idea is to treat everything as quantitative. 
# Hence the Beta in lm is a modifier of the change in Kozak strength
# This can be applied to all but 0 change ones; 
# For the single variant ones this is identitical to 0,1,2 coding

# "ENST00000263461.5" is significant ribo not significant rna WDR11-001 which is also implicated in cancer => FDR ~ .5
# I should order WDR11-001 and from single ones:  65, 145 FDR .5;   316, 319 FDR < .1 
multi_pval <- c()
# total_number_of_alleles can be greater than max(number_alleles) as it is the sum in two positions
for (i in names(total_alleles[total_alleles > 60*.1 ]) ) { 
 my_index = which ( row.names(ribo_only) == enst_hgnc[grep (i, enst_hgnc), 2])
 # Add other requirements
 if ( length(my_index) ) { 
 index_factor <- rep(0, times=length(sample_labels_ribo))
 # Modify index factor using the corresponding multi_diff
 rna_index_factor <- rep (0, times=length(sample_labels_rna))
 
  all_variants <-  (which(multi_ind[,1] == i))
  for ( j in all_variants) { 
   # Find corresponding ribo 
    all_ind <- grep("NA",multi_ind[j,-1], value=T)
    allele_num <- rle(all_ind)$lengths
    ind_unique <- unique(all_ind)
    ind_unique <- sub("NA", "GM", ind_unique)
    ribo_index <- grep(paste(ind_unique , collapse="|"), sample_labels_ribo)
    ribo_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_ribo, value=T)
    index_factor[ribo_index] = index_factor[ribo_index] + (multi_diff[j] * allele_num[match(ribo_index_values, ind_unique)])
    
    rna_index <-  grep(paste(ind_unique , collapse="|"), sample_labels_rna)
    rna_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_rna, value=T)
    rna_index_factor[rna_index] <- rna_index_factor[rna_index] +  (multi_diff[j] * allele_num[match(rna_index_values, ind_unique)] )
  }
# print (index_factor)
 multi_pval <- c(multi_pval, summary(lm(ribo_only$E[my_index,] ~ index_factor, weights=ribo_only$weights[my_index,]))$coefficients[2,4])
  boxplot(ribo_only$E[my_index,]~ index_factor, ylab= "Ribosome Occupancy", xlab="Kozak Strength" , names=unique(sort(round(index_factor, digits=2))) )
legend("topright", paste("p-val = ",  round(summary.lm(lm(ribo_only$E[my_index,] ~ index_factor, weights=ribo_only$weights[my_index,]))$coefficients[2,4], digits=4 ), sep="" ), inset=0.05, bty= "n" )
  boxplot(rna_only$E[my_index,]~rna_index_factor, ylab="RNA Expression", xlab= "Kozak Strength", names=unique(sort(round(index_factor, digits=2))))
legend("topright", paste("p-val = ",  round(summary.lm(lm(rna_only$E[my_index,] ~ rna_index_factor, weights=rna_only$weights[my_index,]))$coefficients[2,4], digits=2 ), sep="" ), inset=0.05, bty= "n" )
summary.lm(lm(ribo_only$E[my_index,] ~ index_factor, weights=ribo_only$weights[my_index,]))
summary.lm(lm(rna_only$E[my_index,] ~ rna_index_factor, weights=rna_only$weights[my_index,]))
 }
}

# Kozak Diff is WT - VARIANT
# Positive score means WT kozak strength is better
# Negative score means VARIANT kozak Strength is better
# grep ("10847|19240", sample_labels) gives the index of the expr value
# Go over the individuals, grep the samples and calculate difference in mean

kozak_diff <- kozak_merge$V2[!multi] - kozak_merge$V4[!multi]
single_var_ind <- kozak_var_ind[!multi,]
ribo_diff <- c()
list_of_pval <- c()
for (i in 1:length(single_var_ind[,1])) { 
  all_ind <- grep("NA",single_var_ind[i,-1], value=T)
  allele_num <- rle(all_ind)$lengths
  ind_unique <- unique(all_ind)
  ind_unique <- sub("NA", "GM", ind_unique)
  # CHECK ENST EQUIVALENT IS PRESENT IN V3$E, IF NOT ADD NA
  my_index <- which ( row.names(ribo_only) == enst_hgnc[grep (single_var_ind[i,1], enst_hgnc),2] )
  
  ribo_index <- grep(paste(ind_unique , collapse="|"), sample_labels_ribo)
  ribo_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_ribo, value=T)

  rna_index <-  grep(paste(ind_unique , collapse="|"), sample_labels_rna)
  rna_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_rna, value=T)
  if (length(my_index) & length(ribo_index) %% 50 != 0 & length(ind_unique) > 1 & length(ind_unique) < 29 & sum(allele_num) > 60*0.1 ) {
# Ribo_Diff is a weighted mean difference
    ribo_diff <- c(ribo_diff, weighted.mean(ribo_only$E[my_index,ribo_index],ribo_only$weights[my_index,ribo_index] ) - weighted.mean(ribo_only$E[my_index, -ribo_index], ribo_only$weights[my_index, -ribo_index] ))
    index_factor <- rep(0,times=length(sample_labels_ribo))
    index_factor[ribo_index] <- allele_num[match(ribo_index_values, ind_unique)]
    my_pval <- summary(lm(ribo_only$E[my_index,] ~ index_factor, weights=ribo_only$weights[my_index,]))$coefficients[2,4]
    list_of_pval <- c(list_of_pval, my_pval)
    if ( my_pval < 0.01) {
    pdf(file=paste("~/Google_Drive/Manuscript Figures/Kozak_Analysis/Ribo_", row.names(ribo_only)[my_index],".pdf", sep=""), width=5, height=5   )
    boxplot(ribo_only$E[my_index,]~ index_factor, ylab= "Ribosome Occupancy", xlab="Allele Number" , names=unique(sort(index_factor)), main=row.names(ribo_only)[my_index] )
    legend("topright", paste("p-val = ",  signif(my_pval, digits=3 ), sep="" ), inset=0.05, bty= "n" )
    dev.off()
    rna_index_factor <- rep (0, times=length(sample_labels_rna))
    rna_index <-  grep(paste(ind_unique , collapse="|"), sample_labels_rna)
    rna_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_rna, value=T)
    rna_index_factor[rna_index] <-  allele_num[match(rna_index_values, ind_unique)] 
    rna_pval <- summary(lm(rna_only$E[my_index,]~ rna_index_factor, weights=rna_only$weights[my_index,]))$coefficients[2,4]
    pdf(file=paste("~/Google_Drive/Manuscript Figures/Kozak_Analysis/RNA_", row.names(rna_only)[my_index], ".pdf",  sep=""), width=5, height=5   )
    boxplot(rna_only$E[my_index,]~ rna_index_factor, ylab= "RNA Expression", xlab="Allele Number" , names=unique(sort(index_factor)), main= row.names(rna_only)[my_index] )
    legend("topright", paste("p-val = ",  signif(rna_pval, digits=3 ), sep="" ), inset=0.05, bty= "n" )  
    dev.off()
    }
  }
  else { 
    ribo_diff <- c(ribo_diff,NA)
    list_of_pval <- c(list_of_pval,NA)
  }
}

# Type of mutation and associated ribosome diff
ribo_diff_mutation_type = data.frame(DIFF=ribo_diff, TYPE=all_kozak_score_variants[!multi,2])
ribo_diff_mutation_type = ribo_diff_mutation_type[!is.na(ribo_diff_mutation_type[,1]),]
i1 <- ribo_diff_mutation_type[,2] %in% names(table(ribo_diff_mutation_type[,2]))[table(ribo_diff_mutation_type[,2]) > 2]
boxplot(ribo_diff_mutation_type[i1,1] ~ factor(ribo_diff_mutation_type[i1,2]))
abline(h=0)

# When we look at the pvalues, All the significant changes have Kozak strength changes in the upper half
# The direction of the effect is inconsistent in some
# MAF 10% 101 , c1 < 1 => 14; c1 < .5 => 10; c1 < .01 => 4; c1 < .25 =>6
# FDR < .1
sum(!is.na(list_of_pval))
length(which(p.adjust(list_of_pval) < .05))
significant_ribo_diff <- which(p.adjust(list_of_pval) < .05)
color_by_pval <- rep(0, length(list_of_pval))
pval_cutoff <- max(list_of_pval[significant_ribo_diff])
color_by_pval[list_of_pval <= pval_cutoff] <- 1
plot(ribo_diff, kozak_diff, cex=0.65, pch=19, xlim=c(-1.2, 1.2), tck = .02, col=c("Black", "Red")[as.factor(color_by_pval)])
abline(v=c(0), h=c(0,log(2), -log(2)))

# Some figure to show that the significantly different Ribo Diff Ones have significant effect on Kozak
# We can just state this; when we take transcripts with ribo difference significant at 5% FDR
# Wilcox.test p-value is 0.02
boxplot(abs(kozak_diff[significant_ribo_diff]), abs(kozak_diff[!is.na(p.adjust(list_of_pval))]))
wilcox.test(abs(kozak_diff[significant_ribo_diff]), abs(kozak_diff[!is.na(p.adjust(list_of_pval))]))


# Take the translation efficiency table and for each position look for significant diff with Kruskal
kozak_seq_score_table <- read.table('~/project/CORE_DATAFILES/Kozak_IDs_PWM_Strand_Seq_HGNC_TE_Ribo.bed')
# Translation Efficiency and Ribosome Occupancy are different somewhat
f <- function(s, letter) strsplit(s, "")[[1]][letter]
g <- function(s) strsplit(s, "")[[1]][c(4,10)]
for ( j in c(1:6, 10:11)) { 
seq_factor <- sapply(as.character(kozak_seq_score_table$V7), f, letter=j)
k1 <- kruskal.test(kozak_seq_score_table$V9 ~ as.factor(seq_factor))
print (k1$p.value * 8)
print(kruskalmc(kozak_seq_score_table$V9 ~ as.factor(seq_factor), probs=.01))
boxplot(kozak_seq_score_table$V10 ~ as.factor(seq_factor), varwidth=T, ylim=c(-1,1), notch=T, range=.001, cex=.2)
#boxplot(kozak_seq_score_table$V10 ~ as.factor(seq_factor), varwidth=T, ylim=c(4,6), notch=T, range=.001, cex=.2)
}
# Formalize the across gene analysis and maybe do the testing with TE


###


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

# This function would crash if there is only one data point in a cluster
output_clusters <- function (k, data_original, hc, names)
{
  desired_clusters = k
  dir <- paste("./Cluster" ,k, sep= "")
  dir.create(dir)
  cluster_membership <- cutree(hc, k=desired_clusters)
  for (i in 1:k ) {
    if (length(data_original[which(cluster_membership == i),]) > 3) {
      filename_means <-
        paste(round(colMeans(data_original[which(cluster_membership == i),])),
              collapse= "_")
    }
    else {
      filename_means <- paste(round (data_original[which(cluster_membership ==
                                                           i),]), collapse= "_")
    }
    writeLines(as.character(names[which(cluster_membership == i)]),
               paste(dir, paste("cluster", i, filename_means, sep="_"), sep ="/") )
  }
}

### HEATMAP FIGURE 
# ranked_data_na_omit <- apply(as.matrix(na.omit(cbind(hela_rna_mean, mass_spec_data$HeLa_prot, mass_spec_data$Guo_FP))), 2,
#                              function(x){rank(x, na.last=F, ties.method="min")
#                              })
# omitted_lines <- as.vector(na.action(na.omit(cbind(
#   hela_rna_mean, mass_spec_data$HeLa_prot, mass_spec_data$Guo_FP))))
# names_post_omit <- mass_spec_data$REFSEQ[-omitted_lines]
# ranked_data_na_omit_normalized <- ranked_data_na_omit - apply(ranked_data_na_omit, 1, median)

#### SAMPLE HEATMAP
#hum_heat_norm <- heatmap.2(ranked_data_na_omit_normalized[1:200,], Colv=NA,
#labRow=NA, scale="none",
#col = redgreen(128),labCol=
#c("RNA", "Prot", "FP"), trace="none", key=TRUE, density.info="none",
#keysize=1, lmat=rbind( c(3, 4), c(2,1)), lwid=c(1.25, 0.5), lhei= c(0.2,1),
#hclustfun= hclust_ward )
#####

# # Uncomment region to output various clusters
#hcr <- hclust(dist(ranked_data_na_omit_normalized))
#output_clusters(4, ranked_data_na_omit_normalized, hcr, names_post_omit)
##

### DIFFERENCE IN VARIATION
### THIS APPROACH IS COMMENTED OUT
# 
# # Test for difference in variance per gene across individuals
# # Partition sum of squares to estimate between individual variance take out variance component due to rep to rep variance
# # Test statistic is the difference between F values, significance is by permutation testing
# # Need to think about issues with respect to degrees of freedom associated with the calculated F-value
# 
# #joint_expression_matrix <- merge(v$E[,replicate_present][,sample_labels[replicate_present] %in% sample_id[replicate_present_rnaseq]], v2$E[,replicate_present_rnaseq][,sample_id[replicate_present_rnaseq] %in% sample_labels[replicate_present]], by="row.names")
# #joint_expression_matrix <- joint_expression_matrix[,-1]
# joint_expression_matrix <- joint_expression_common$E
# #gene_names_joint_expression_matrix <- merge(v$E[,replicate_present][,sample_labels[replicate_present] %in% sample_id[replicate_present_rnaseq]], v2$E[,replicate_present_rnaseq][,sample_id[replicate_present_rnaseq] %in% sample_labels[replicate_present]], by="row.names")[,1]
# #sample_id_all <- unlist(strsplit(colnames(joint_expression_matrix), split= "_"))
# #sample_id_all <- sample_id_all[grep("GM", sample_id_all)]
# sample_id_all <- sample_labels_joint_common
# sample_id_all[1:51] <- paste(sample_id_all[1:51], "RNA", sep="_")
# sample_id_all[52:84] <- paste(sample_id_all[52:84], "Ribosome_Profiling", sep="_")
# 
# # GO over each gene. Calculate F-value for ribo and rna separately
# # One issue is that mean diff is highly correlated with p-value.
# # The higher the mean diff, the higher the p-value
# Mean_diff <- apply(joint_expression_matrix, 1, function(x) {mean(x[1:52] - mean(x[53:84]))} ) 
# joint_expression_matrix_mean_subtracted <- t(apply(joint_expression_matrix, 1, function(x) {c(x[1:52]-mean(x[1:52]), x[53:84]-mean(x[53:84])  )}))
# F_diff <- apply(joint_expression_matrix_mean_subtracted, 1, function (x) { 
#   (summary (aov(x[52:84] ~ as.factor(sample_id_all[52:84])))[[1]]$F[1] ) - (summary(aov(x[1:51] ~ as.factor(sample_id_all[1:51])))[[1]]$F[1])
#   })
# F_diff_pval <- c()
# individuals <- unique(sample_id_all)
# # Even with 100 permutations, this is extremely slow
# # If I have the time, I will rewrite this with apply
# # Something like: 
# #perm_values <- apply(joint_expression_matrix_mean_subtracted, 1, function (x) { 
# #  (summary (aov(x[ribo] ~ as.factor(sample_id_all[ribo])))[[1]]$F[1] ) - (summary(aov(x[!ribo] ~ as.factor(sample_id_all[!ribo])))[[1]]$F[1])
# #})
# # FINISH THE F_PVAL CALCULATION
# for (i in 1:length(F_diff)) {
#     # Need some permutation scheme-
#     # Go over individuals and assign the ribo/rna label
#     perm_values <- c()
#     for (k in 1:100) { 
#     ribo <- c(rep(FALSE, 51), rep(TRUE, 84-51))
#     for (j in 1: length(individuals)) {
#       if (runif(1) > 0.5) { 
#         ribo[sample_id_all== individuals[j]] <- !ribo[sample_id_all== individuals[j]]  
#       }  
#     }
#     ribo_F_perm <- summary (aov(joint_expression_matrix_mean_subtracted[i,ribo] ~ as.factor(sample_id_all[ribo])))[[1]]$F[1]
#     rna_F_perm <- summary(aov(joint_expression_matrix_mean_subtracted[i,!ribo] ~ as.factor(sample_id_all[!ribo])))[[1]]$F[1]
#     perm_values <- c(perm_values, ribo_F_perm - rna_F_perm)
#     }
#     p1 <- min (length(which( perm_values > F_diff[i] ) ) /100 , length(which( perm_values < F_diff[i] ) ) /100 )
#     F_diff_pval <- c(F_diff_pval, 2*p1)
# }
# hist(F_diff, 300)
# hist(F_diff, 300, xlim=c(-20,20))
# hist(F_diff[F_diff_pval<0.05], 200, xlab="F_Diff", main="Significant Differences in Variance")
# hist(Mean_diff, 100)
# plot(F_diff_pval, Mean_diff, pch=19, cex=.2, xlab="Pvalue", ylab="Mean Expression Difference")
# # RNA expression is more variable for most things consistent with previous reports that suggests buffering
# # Extract_ids and run FuncAssociate. 
# 
# save (F_diff, file= "~/project/CORE_DATAFILES/FValue_Differences")
# save (F_diff_pval, file="~/project/CORE_DATAFILES/FValue_Differences_Pvals")
# save (joint_expression_matrix_mean_subtracted, file="~/project/CORE_DATAFILES/Joint_Expression_Matrix")
# save(low_pval_indices, file="~/project/CORE_DATAFILES/low_pval_indices")
# # For the set of transcripts where F_diff_pval < 0.01, do more extensive permutation -- run this overnight
# low_pval_indices <- which(F_diff_pval < 0.05)
# # rna_variable <- F_diff[low_pval_indices] < 0
# # ribo_variable <- F_diff[low_pval_indices] > 0
# # write.table(gene_names_joint_expression_matrix[low_pval_indices][rna_variable], file="~/Desktop/RNA_variable.txt", row.names=F)
# # write.table(gene_names_joint_expression_matrix[low_pval_indices][ribo_variable], file="~/Desktop/Ribo_variable.txt", row.names=F)
# # write.table(gene_names_joint_expression_matrix, file="~/Desktop/All_Tested_IDs", row.names=F)
# for (i in low_pval_indices) { 
#   perm_values <- c()
#   for (k in 1:10000) { 
#     ribo <- c(rep(FALSE, 51), rep(TRUE, 84-51))
#     for (j in 1: length(individuals)) {
#       if (runif(1) > 0.5) { 
#         ribo[sample_id_all== individuals[j]] <- !ribo[sample_id_all== individuals[j]]  
#       }  
#     }
#     ribo_F_perm <- summary (aov(joint_expression_matrix_mean_subtracted[i,ribo] ~ as.factor(sample_id_all[ribo])))[[1]]$F[1]
#     rna_F_perm <- summary(aov(joint_expression_matrix_mean_subtracted[i,!ribo] ~ as.factor(sample_id_all[!ribo])))[[1]]$F[1]
#     perm_values <- c(perm_values, ribo_F_perm - rna_F_perm)
#   }
#   p1 <- min (length(which( perm_values > F_diff[i] ) ) /10000 , length(which( perm_values < F_diff[i] ) ) /10000 )
#   F_diff_pval[i] <-  2*p1  
# }
#####################################
# END OF OLD VARIATION ANALYSIS


# Calculate correlation between APPRIS USAGE AND TRANSLATION
cuff_ratios <- cuff_ratios[cuff_ratios$HGNC %in% row.names(v3) ,]
# Batch correct APPRIS AND THE REMAINING EXPRESSIONS
appris_expression <- as.matrix(cuff_ratios[,seq(2,ncol(cuff_ratios)-1, by=3)])
nonappris_expression <- as.matrix(cuff_ratios[,seq(3,ncol(cuff_ratios), by=3)]  )
cuff_ids <- cuff_ratios[, 179]

cuff_batch <- c (rep(1,16), rep(2, 18), rep(3,25)) 
s1 <- unlist(strsplit(colnames(appris_expression), "_"))
sample_labels_cuff <- s1[grep('GM', s1)]
cuff_design <- model.matrix(~sample_labels_cuff)
# Need to remove genes with all zeros
s1 <- apply(nonappris_expression, 1, function(x){sum(x==0)}) == 59
s2 <- apply(appris_expression, 1, function(x){sum(x==0)})== 59
nonappris_expression <- nonappris_expression[!(s1 | s2), ]
appris_expression <- appris_expression[!(s1 | s2), ]
cuff_ids <- cuff_ids[!(s1 | s2)]
# Filter samples with no variation or very low/high ratio in all samples

batch_removed_appris <- ComBat (appris_expression, batch=cuff_batch, mod=cuff_design)
batch_removed_nonappris <- ComBat (nonappris_expression, batch=cuff_batch, mod=cuff_design)
cuff_ratios_batch_corrected <- batch_removed_appris / (batch_removed_appris + batch_removed_nonappris)
# Even after batch correction GEUVADIS CORRELATES POORLY WITH SNYDER
row.names(cuff_ratios_batch_corrected) <- cuff_ids
one_frames <- read.table('~/project/CORE_DATAFILES/gencode.v15.protein_coding.one_frame_only.ids', header=F)
one_frames <- merge (one_frames, ensg_hgnc, by.x="V1", by.y="ENSG")
# Most of the one-frame transcripts have poor expression and hence are filtered out
cuff_ratio_one_frame_prots <- merge(one_frames, cuff_ratios_batch_corrected, by.x="HGNC", by.y="row.names")

cuff_replicate_median <- apply (cuff_ratio_one_frame_prots[,-c(1,2)], 1, function(x) {
  aggregate(x, by= list(as.factor(sample_labels_cuff)), median)  
} )
# We can try to use everything or we can subset to Snyder
# Calculate per sample median ratio

# redo using Hagen's ratio appris sum cor ribo nonappris to ribo
ribo_replicate_mean <- apply(v3$E[row.names(v3)%in% cuff_ratio_one_frame_prots[,1],type=="Ribo"],1 ,function(x){
  aggregate(x, by=list(as.factor(sample_labels_joint[type=="Ribo"])), mean)  
})

ribo_replicate_mean <- lapply(ribo_replicate_mean, function(x){x <- x[-c(8,10,12,15,16,23),]})
cuff_replicate_median <- lapply(cuff_replicate_median, function(x){x <- x[-10,]})

c2 <- mapply (cbind, cuff_replicate_median, ribo_replicate_mean, SIMPLIFY=F)
across_ind_cuff_cor <- as.numeric(lapply(c2, function(x){ cor.test(x[,2], x[,4], method="spearman")$p.value }))
hist(across_ind_cuff_cor, 50)
# There are very few cases of signal at p < 0.005. The two cases with FDR < .25 are not UTR variants
# There is a bunch RP stuff that messes this up. Better leave this as a future project
# "ENSG00000006453"    "BAIAP2L1-001" 
#############################################################

#########OLD KOZAK ALLELIC EFFECT PLOTS:
# Look at the effect of number of alleles in the significant ones
# Look at the multiple ones; the analysis of multiple ones could be similar to number of alleles
# Multiple variants always occur on different positions
# It seems like that all the significant ones are high allele frequency. We might exclude singletons from testing

# Allelic effect boxplots --- KOZAK_VAR_IND <=> SINGLE_VAR_IND
significant_transcripts <- kozak_var_ind[significant_ribo_diff, 1]
# 316 -> ENST00000362031.4  SNX6 
# 319 -> ENST00000366628.4  NTPCR
# 65 -> ENST00000254759.3 COQ3 -> This another mitochondrial gene with alternative translation initiation site that is Kozak dependent switch
# 145 -> ENST00000295641.10 STK11IP -- Let's not order this one
# ODC1 might also be interesting in the future but the FDR is too high to follow up now. 
significant_gene_ids <- enst_hgnc[grep(paste(significant_transcripts, collapse="|"), enst_hgnc),2]
ribo_indicies <- grep(paste(significant_gene_ids, collapse="|"), row.names(ribo_only))
# 316, 319 have clear boxplots with expected direction
# 500 is also expected direction but the boxplot is less clear
# 52, 234 are unexpected direction -> 234 boxplot is not clear with low expression
# 52, 234, 500 have the same effect at RNA
# 316, 319 are definitely at TE; 
# When we get to 50% FDR, 65, 145 (no Kozak effect) is also TE

for ( i in significant_ribo_diff) { 
  all_ind <- grep("NA",kozak_var_ind[i,-1], value=T)
  allele_num <- rle(all_ind)$lengths
  ind_unique <- unique(all_ind)
  ind_unique <- sub("NA", "GM", ind_unique)
  my_index <- which ( row.names(ribo_only) == enst_hgnc[grep (kozak_var_ind[i,1], enst_hgnc),2] )
  
  ribo_index <- grep(paste(ind_unique , collapse="|"), sample_labels_ribo)
  ribo_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_ribo, value=T)
  
  rna_index <-  grep(paste(ind_unique , collapse="|"), sample_labels_rna)
  rna_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_rna, value=T)
  rna_index_factor <- rep (0, times=length(sample_labels_rna))
  rna_index_factor[rna_index] <- allele_num[match(rna_index_values, ind_unique)]
  
  index_factor <- rep(0,times=length(sample_labels_ribo))
  index_factor[ribo_index] <- allele_num[match(ribo_index_values, ind_unique)]
  boxplot(ribo_only$E[i,]~ index_factor, ylab= "Ribosome Occupancy", xlab="Number of Alleles", main = enst_hgnc[grep(kozak_var_ind[i,1], enst_hgnc), 2] )
  legend("topright", paste("p-val = ",  signif(summary.lm(lm(ribo_only$E[my_index,] ~ index_factor, weights=ribo_only$weights[my_index,]))$coefficients[2,4], digits=2 ), sep="" ), inset=0.05, bty= "n" )
  boxplot(rna_only$E[i,]~rna_index_factor, ylab="RNA Expression", xlab="Number of Alleles", main = enst_hgnc[grep(kozak_var_ind[i,1], enst_hgnc), 2])
  legend("topright", paste("p-val = ",  signif(summary.lm(lm(rna_only$E[my_index,] ~ rna_index_factor, weights=rna_only$weights[my_index,]))$coefficients[2,4], digits=2 ), sep="" ), inset=0.05, bty= "n" )
  summary.lm(lm(ribo_only$E[my_index,] ~ index_factor, weights=ribo_only$weights[my_index,]))
  summary.lm(lm(rna_only$E[my_index,] ~ rna_index_factor, weights=rna_only$weights[my_index,]))
}

##### OLD RNA-SEQ ONLY VOOM -- NOT USED ANYMORE
#pdf("Mean_Variance_Modelling_RNASEQ.pdf")
v2 <- voom (all_rnaseq_counts, design_rnaseq, plot=T)
#dev.off()
replicate_present_rnaseq <- duplicated(sample_id) | duplicated(sample_id, fromLast = TRUE)
norm_dd_rnaseq <- dist( t( v2$E) )
hc_rnaseq <- hclust( norm_dd_rnaseq)
# BATCH CORRECTION IS ESSENTIAL HERE
rnaseq_batch <- c (rep(1,18), rep(2, 26), rep(3,25), rep(4, 15) ) 
batch_removed <- ComBat (v2$E, batch=rnaseq_batch, mod=design_rnaseq)
# Update v2$E with batch_removed; this keeps weights
v2$E <- batch_removed
#write.table(batch_removed,
#file ="TMM_VarianceMeanDetrended_CPM_GT1_in40_BatchRemoved_RNASeq_Expression",
#sep="\t", row.names=as.character(rna_seq_normalized_with_ids[,1])[!abs(stdres(fit.rna)) > 3])
batch_dd_rnaseq <- dist( t( batch_removed) )
hc_batch_rnaseq <- hclust( batch_dd_rnaseq)

###############
# OLD_RIBO_ONLY VOOM NOT USED ANYMORE
#pdf("Mean_Variance_Modeling_voom.pdf")
v <- voom(cds_counts,design, plot=T)
#dev.off()
row.names(v) <- CDS_IDs[isexpr]
norm_dd <- dist(t (v$E ) ) 
norm_hc <- hclust (norm_dd)
replicate_present <- duplicated(sample_labels) | duplicated(sample_labels, fromLast = TRUE)
norm_dd_rep <- dist(t (v$E[,replicate_present] ) )
norm_hc_rep <- hclust (norm_dd_rep)
#write.table(v$E, file ="TMM_VarianceMeanDetrended_CPM_GT1_in36_RiboProfiling_Expression", 
#row.names=CDS_IDs[isexpr], sep="\t")

### SVA- Correction; we decided to not use this correction 
# SVA/Batch Correction for Joint data
#### COMMENTED OUT SVA NORMALIZATION
# # We can look sva correlations with sequencing depth, batch, etc
# # GM19139 - HAS no RNA index 124
# # GM19139 is most similar to GM19137 - For SVA purposes use GM19137
# full_design_nosingular <- model.matrix(~sample_labels_joint[-124] + type[-124])
# #svobj_joint <- sva (v3$E[,-124], mod=full_design_nosingular, B=50)
# svobj_joint <- sva (v3$E[,-124], mod=full_design_nosingular, B=50, n.sv=3)
# fit_joint <- lmFit(v3$E[,-124], svobj_joint$sv) 
# fit2 <- eBayes(fit_joint)
# sst <- rowSums(v3$E^2)
# ssr <- sst-fit2$df.residual*fit2$sigma^2
# hist((ssr/sst), 50)
# quantile(ssr/sst)

# With 20 components
# 0%         25%         50%         75%        100% 
# 0.006377750 0.009714234 0.013078694 0.023458592 0.408221116 
# With 3 components
# 0.0002159566 0.0072133405 0.0084482581 0.0103326042 0.1358229828 
# With 10 components
# 0.002256528 0.008297705 0.009897312 0.014138372 0.265767843 

# norm_expr_joint <- residuals(fit_joint, v3$E[,-124])
# cor_coefs <- c()
# for (i in 1:133) {
#   cor_coefs[i] <- (cor.test(norm_expr_joint[,i], v3$E[,-124][,i]))$estimate
# }
# hist(as.numeric(cor_coefs), 40, xlab="Pearson Cor Coef", main="SVA vs No-SVA Expression")
# hist(as.numeric(cor_coefs)[type[-124]=="Ribo"], 40, xlab="Pearson Cor Coef", main="SVA vs No-SVA Expression-Ribo")
# hist(as.numeric(cor_coefs)[type[-124]=="RNA"], 40, xlab="Pearson Cor Coef", main="SVA vs No-SVA Expression-RNA")
# max(as.numeric(cor_coefs))
# plot(norm_expr_joint[,94], v3$E[,-124][,94], pch=19, xlab="SVA", ylab="QN_TMM_Voom", main=colnames(v3$E)[94], cex=.2)
# plot(norm_expr_joint[,126], v3$E[,-124][,126], pch=19, xlab="SVA", ylab="QN_TMM_Voom", main=colnames(v3$E)[127], cex=.2)
# write.table(norm_expr_joint[,1:84], file =paste (data_dir, "Top3_SVA_Removed_QN_FullModel_RNA_Expression", sep=""),             
#             row.names=joint_count_ids, sep="\t")
# write.table(norm_expr_joint[,85:133], file =paste (data_dir, "Top3_SVA_Removed_QN_FullModel_RiboProfiling_Expression", sep=""),             
#             row.names=joint_count_ids, sep="\t")

