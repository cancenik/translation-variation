library("edgeR")
library("hexbin")
library("limma")
library("qtl")
library ("sva")
library("MASS")
#library ("isva")
# I decided not to use isva after initial tests

## NOTES
## Test for difference in variance
## Calculate F-value using anova on lm for each gene. Compare Ribo RNA
## Permute the Ribo & RNA labels within each individual. If we have 4 reps of RNA 2 reps of Ribo
## We can get 4 reps of Ribo and 2 reps of RNA

# Number of bases covered is not a good measure
# The number of species correlates better with length but
# Total read number gives tighter clustering 
# Should use total read number for the analysis

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


## RIBOSEQ_COUNTS
# CDS SPECIES
species <- read.table(paste(data_dir, "Reformatted_Species_Counts_All_Libraries.tsv", sep=""),header=T)
CDS_species <- split (species, species$REGION)[[2]]
CDS_species <- CDS_species[,grep("Counts", colnames(CDS_species))]

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
colSums(CDS_Counts)
dim (CDS_Counts[keep(CDS_Counts, 100),])
cm <- cor(CDS_Counts[keep(CDS_Counts, 100),])
dd <- dist (t(log10(CDS_Counts[keep(CDS_Counts, 100),]+1)) )
hc <- hclust (dd, "ward") 
hc <- hclust (dd)

# SPECIES
colSums (CDS_species)
dim(CDS_species[keep(CDS_species, 100), ])
cm <- cor(CDS_species[keep(CDS_species, 100),])
dd <- dist(t (log10(CDS_species[keep(CDS_species, 100), ]+1) ) ) 
hc <- hclust(dd, "ward")
hc <- hclust(dd)


################### THIS SECTION NEEDS TO BE FIXED
## SPECIES TO COUNTS COMPARISON
species_sum <- rowSums(CDS_species)
cds_count_sum <- rowSums(CDS_Counts)
ratios <- cds_count_sum/species_sum
plot(log10(cds_count_sum+1), log10(species_sum+1), cex=0.2)
#ratios_ids <- data.frame(ID=as.vector(m1[,1]), ratio=ratios)
#ratios_dataframe <- data.frame(ID=as.vector(m1[,1]), ReadCount=log10(cds_count_sum+1) , SpeciesCount=log10(m1_species_sum+1) )

#### Perform loess regression between read_count to species_count
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
plot(log10(polyA_mean), log10(RZ_mean), cex=0.2, pch=19, ylab="log10(RiboZero_ReadCount)", xlab="log10(PolyA_ReadCount)", main="Comparing RNASeq Methods")
fit.rna <- lm(log10(RZ_mean)~log10(polyA_mean))
# Identify outliers
# rna_seq_normalized$ID[abs(stdres(fit.rna)) > 3]
outlier_colors <- rep("Black", length(log10(polyA_mean)))
outlier_colors[abs(stdres(fit.rna)) > 3] <- "Red"
plot(log10(polyA_mean),stdres(fit.rna) , pch=19, cex=0.2, xlab="Log10(PolyA_Reads)", ylab="Standardized Residuals", col=outlier_colors)
polyA_RZ_inconsistent <- (abs(stdres(fit.rna)) > 3)
all_rnaseq_counts <- all_rnaseq_counts[!polyA_RZ_inconsistent,]
row.names(all_rnaseq_counts) <- (geuvadis_CDS$ID[rnaexpr])[!polyA_RZ_inconsistent]

all_rnaseq_counts <- calcNormFactors (all_rnaseq_counts, method= "TMM")
all_rnaseq_counts$samples
cor (all_rnaseq[rnaexpr,], method="spearman")
sample_id <- unlist(strsplit(colnames(all_rnaseq_counts), split= "_"))
sample_id <- sample_id[grep("GM", sample_id)]
design_rnaseq <- model.matrix(~sample_id)

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

# VOOM- TMM Normalization- Ribosome Profiling
cds_counts <- DGEList(counts=CDS_Counts)
isexpr <- rowSums(cpm(cds_counts) > 1) >= 36
cds_counts <- cds_counts[isexpr,]
row.names(cds_counts) <- CDS_IDs[isexpr]
dim(cds_counts)
cds_counts <- calcNormFactors (cds_counts, method= "TMM")
cds_counts$samples
s1 <- unlist(strsplit(colnames(CDS_Counts), "_"))
sample_labels <- s1[grep('GM', s1)]
table(sample_labels)
design <- model.matrix(~sample_labels)
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

#

### An alternative way to process RNA-Seq Ribo-Seq jointly
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
v3 <- voom(joint_counts, full_design, plot=T)
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

# SVA/Batch Correction for Joint data
# p1 <- paste (covariates$Date_Cells_Frozen, covariates$Date_Ribosome_Footprint, covariates$Data_Gel_Purified, sep="_")
# p1 <- covariates$Data_Gel_Purified
# p2 <- paste (covariates$Sequencing_ID, "Counts", sep="_")
# p1 <- p1[p2 %in% colnames(v$E)]
# p2 <- p2[p2 %in% colnames(v$E)]
# p1 <- p1[sort (p2, index.return=T)$ix]
mod <- model.matrix(~as.factor(sample_labels), data=as.data.frame(v$E))
svobj <- sva (v$E, mod=mod, B=20)
fit = lmFit(v$E, svobj$sv)
norm_expr <- residuals (fit, v$E)
row.names(v) <- CDS_IDs[isexpr]
save (v, file='~/project/CORE_DATAFILES/RiboProfiling_TMM_Voom_normalizedEList.RData')
# Replace v$E with norm_expr
sva_norm_weights <- v
sva_norm_weights$E <- norm_expr
save (sva_norm_weights, file='~/project/CORE_DATAFILES/RiboProfiling_TMM_Voom_normalized_SVA_EList.RData')
sva_dd <- dist( t ( norm_expr[,replicate_present]))
sva_hc <- hclust (sva_dd) 

# We can look sva correlations with sequencing depth, batch, etc
#GM19139 - HAS no RNA index 124
# GM19139 is most similar to GM19137 - For SVA purposes use GM19137
full_design_nosingular <- model.matrix(~sample_labels_joint[-124] + type[-124])
svobj_joint <- sva (v3$E[,-124], mod=full_design_nosingular, B=50)
fit_joint <- lmFit(v3$E[,-124], svobj_joint$sv) 
norm_expr_joint <- residuals(fit_joint, v3$E[,-124])
cor_coefs <- c()
for (i in 1:133) {
  cor_coefs[i] <- (cor.test(norm_expr_joint[,i], v3$E[,-124][,i]))$estimate
}
hist(as.numeric(cor_coefs), 40, xlab="Pearson Cor Coef", main="SVA vs No-SVA Expression")
hist(as.numeric(cor_coefs)[type[-124]=="Ribo"], 40, xlab="Pearson Cor Coef", main="SVA vs No-SVA Expression-Ribo")
hist(as.numeric(cor_coefs)[type[-124]=="RNA"], 40, xlab="Pearson Cor Coef", main="SVA vs No-SVA Expression-RNA")
max(as.numeric(cor_coefs))
plot(norm_expr_joint[,94], v3$E[,-124][,94], pch=19, xlab="SVA", ylab="QN_TMM_Voom", main=colnames(v3$E)[94], cex=.2)
plot(norm_expr_joint[,126], v3$E[,-124][,126], pch=19, xlab="SVA", ylab="QN_TMM_Voom", main=colnames(v3$E)[127], cex=.2)


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

# OUTPUT SVA RNA and normal RNA
write.table(v3$E[,1:84], file =paste (data_dir, "TMM_VarianceMeanDetrended_QN_FullModel_RNA_Expression", sep=""), 
row.names=joint_count_ids, sep="\t")
write.table(v3$E[,85:134], file =paste (data_dir, "TMM_VarianceMeanDetrended_QN_FullModel_RiboProfiling_Expression", sep=""), 
            row.names=joint_count_ids, sep="\t")            
write.table(norm_expr_joint[,1:84], file =paste (data_dir, "Top20_SVA_Removed_QN_FullModel_RNA_Expression", sep=""),             
            row.names=joint_count_ids, sep="\t")
write.table(norm_expr_joint[,85:133], file =paste (data_dir, "Top20_SVA_Removed_QN_FullModel_RiboProfiling_Expression", sep=""),             
            row.names=joint_count_ids, sep="\t")


##### ANALYSIS REQUIRING JOINT RNA-RIBOSEQ WITH REPLICATES
### DIFFERENCE IN VARIATION

# Test for difference in variance per gene across individuals
# Partition sum of squares to estimate between individual variance take out variance component due to rep to rep variance
# Test statistic is the difference between F values, significance is by permutation testing
# Need to think about issues with respect to degrees of freedom associated with the calculated F-value

#joint_expression_matrix <- merge(v$E[,replicate_present][,sample_labels[replicate_present] %in% sample_id[replicate_present_rnaseq]], v2$E[,replicate_present_rnaseq][,sample_id[replicate_present_rnaseq] %in% sample_labels[replicate_present]], by="row.names")
#joint_expression_matrix <- joint_expression_matrix[,-1]
joint_expression_matrix <- joint_expression_common$E
#gene_names_joint_expression_matrix <- merge(v$E[,replicate_present][,sample_labels[replicate_present] %in% sample_id[replicate_present_rnaseq]], v2$E[,replicate_present_rnaseq][,sample_id[replicate_present_rnaseq] %in% sample_labels[replicate_present]], by="row.names")[,1]
#sample_id_all <- unlist(strsplit(colnames(joint_expression_matrix), split= "_"))
#sample_id_all <- sample_id_all[grep("GM", sample_id_all)]
sample_id_all <- sample_labels_joint_common
sample_id_all[1:52] <- paste(sample_id_all[1:52], "RNA", sep="_")
sample_id_all[53:84] <- paste(sample_id_all[53:84], "Ribosome_Profiling", sep="_")

# GO over each gene. Calculate F-value for ribo and rna separately
# One issue is that mean diff is highly correlated with p-value.
# The higher the mean diff, the higher the p-value
F_diff <- c()
F_diff_pval <- c()
Mean_diff <- c()
individuals <- unique(sample_id_all)
# Even with 100 permutations, this is extremely slow
# If I have the time, I will rewrite this with apply
for (i in 1:dim(joint_expression_matrix)[1]) {
# Subtract mean 
    joint_expression_matrix[i,1:52] <- joint_expression_matrix[i,1:52] - mean(joint_expression_matrix[i,1:52])
    joint_expression_matrix[i,53:84] <- joint_expression_matrix[i,53:84] - mean(joint_expression_matrix[i,53:84])
    ribo_F <- anova (lm(joint_expression_matrix[i,1:52] ~ as.factor(sample_id_all[1:52])))$F[1]
    rna_F <- anova (lm(joint_expression_matrix[i,53:84] ~ as.factor(sample_id_all[53:84])))$F[1]
    F_diff <- c(F_diff,ribo_F-rna_F)
    Mean_diff <- c(Mean_diff, mean(joint_expression_matrix[i,1:52]) - mean(joint_expression_matrix[i,53:84]))
    # Need some permutation scheme-
    # Go over individuals and assign the ribo/rna label
    perm_values <- c()
    for (k in 1:100) { 
    ribo <- c(rep(FALSE, 52), rep(TRUE, 84-52))
    for (j in 1: length(individuals)) {
      if (runif(1) > 0.5) { 
        ribo[sample_id_all== individuals[j]] <- !ribo[sample_id_all== individuals[j]]  
      }  
    }
    ribo_F_perm <- anova (lm(joint_expression_matrix[i,ribo] ~ as.factor(sample_id_all[ribo])))$F[1]
    rna_F_perm <- anova (lm(joint_expression_matrix[i,!ribo] ~ as.factor(sample_id_all[!ribo])))$F[1]
    perm_values <- c(perm_values, ribo_F_perm - rna_F_perm)
    }
    p1 <- min (length(which( perm_values > (ribo_F-rna_F) ) ) /100 , length(which( perm_values < (ribo_F-rna_F) ) ) /100 )
    F_diff_pval <- c(F_diff_pval, 2*p1)
}
hist(F_diff, 300)
hist(F_diff[F_diff_pval<0.05], 300)
hist(Mean_diff, 100)
plot(F_diff_pval, Mean_diff)
# RNA expression is more variable for most things consistent with previous reports that suggests buffering
# Extract_ids and run FuncAssociate. 

#save (F_diff, file= "~/project/CORE_DATAFILES/FValue_Differences")
#save (F_diff_pval, file="~/project/CORE_DATAFILES/FValue_Differences_Pvals")
#save (joint_expression_matrix, file="~/project/CORE_DATAFILES/Joint_Expression_Matrix")

# For the set of transcripts where F_diff_pval < 0.01, do more extensive permutation -- run this overnight
low_pval_indices <- which(F_diff_pval < 0.05)
# rna_variable <- F_diff[low_pval_indices] < 0
# ribo_variable <- F_diff[low_pval_indices] > 0
# write.table(gene_names_joint_expression_matrix[low_pval_indices][rna_variable], file="~/Desktop/RNA_variable.txt", row.names=F)
# write.table(gene_names_joint_expression_matrix[low_pval_indices][ribo_variable], file="~/Desktop/Ribo_variable.txt", row.names=F)
# write.table(gene_names_joint_expression_matrix, file="~/Desktop/All_Tested_IDs", row.names=F)
for (i in low_pval_indices) { 
  perm_values <- c()
  for (k in 1:10000) { 
    ribo <- c(rep(TRUE, 33), rep(FALSE, 84-33))
    for (j in 1: length(individuals)) {
      if (runif(1) > 0.5) { 
        ribo[sample_id_all== individuals[j]] <- !ribo[sample_id_all== individuals[j]]  
      }  
    }
    ribo_F_perm <- anova (lm(as.numeric(joint_expression_matrix[i,ribo]) ~ as.factor(sample_id_all[ribo])))$F[1]
    rna_F_perm <- anova (lm(as.numeric(joint_expression_matrix[i,!ribo]) ~ as.factor(sample_id_all[!ribo])))$F[1]
    perm_values[k] <- ribo_F_perm - rna_F_perm
  }
  p1 <- min (length(which( perm_values > F_diff[i] ) ) /10000 , length(which( perm_values < F_diff[i] ) ) /10000 )
  F_diff_pval[i] <-  2*p1  
}


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


# We can switch v3 with sva_norm 
## SUBSETTING MESSES UP THE DESIGN MATRIX
#all_expr_elist <- cbind(sva_norm_weights[ribo_id_in_rna,replicate_present], v2[rna_id_in_ribo,unique(sort(rna_with_ribo_replicate))])
all_expr_elist <- joint_expression_common
treatment <- relevel(as.factor(type_common),ref="RNA")
all_expr_elist$design <- model.matrix(~as.factor(sample_labels_joint_common)+treatment)

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


## ANALYSIS ON ALL DATA INCLUDING RNA-RIBO IRRESPECTIVE OF REPLICATION
## COMPARISION TO CHRISTINE'S PROTEOMICS
#### Compare absolute levels of protein with rna and ribo -- Overall correlation is better with ribosome profiling
gm12878_prot <- read.csv('~/project/CORE_DATAFILES/GM12878_B0_FDR5_140120_shortforCan.csv', stringsAsFactors=F)
gm12878_prot <- merge(gm12878_prot, ensg_hgnc, by.x="ID.1", by.y="ENSG")

gm12878_ribo <- grep("GM12878", colnames(v))
gm12878_ribo_mean <- apply(v$E[,gm12878_ribo], 1, mean)
gm12878_ribo_mean <- data.frame(HGNC=CDS[isexpr,1], gm12878_ribo_mean)

gm12878_rna <- grep("GM12878", colnames(v2))
gm12878_rna_mean  <- apply(v2$E[,gm12878_rna], 1, mean)
gm12878_rna_mean <- data.frame(HGNC=rownames(v2), gm12878_rna_mean)

ribo_rna_12878 <- merge(gm12878_ribo_mean, gm12878_rna_mean, by="HGNC")
ribo_rna_prot_12878 <- merge(ribo_rna_12878, gm12878_prot, by="HGNC")
ribo_rna_prot_12878 <- ribo_rna_prot_12878[as.numeric(ribo_rna_prot_12878$USE) > 1000 & ribo_rna_prot_12878$gm12878_ribo_mean > 5, ] 

plot(ribo_rna_prot_12878$gm12878_ribo_mean, ribo_rna_prot_12878$gm12878_rna_mean, cex=0.2, pch=19, xlim=c(5,15), ylim=c(3,15))
plot(ribo_rna_12878$gm12878_ribo_mean, ribo_rna_12878$gm12878_rna_mean, cex=0.2, pch=19, xlim=c(5,15), ylim=c(3,15))
plot(ribo_rna_prot_12878$gm12878_ribo_mean,  log10(as.numeric(ribo_rna_prot_12878$USE)), pch=19, cex=.2, xlim=c(4, 15))
plot(ribo_rna_prot_12878$gm12878_rna_mean,  log10(as.numeric(ribo_rna_prot_12878$USE)), pch=19, cex=.2, xlim=c(4, 15))

#Spearman - Grand Mean is BEST
cor.test(ribo_rna_prot_12878$gm12878_ribo_mean, log10(as.numeric(ribo_rna_prot_12878$USE)+1), method="spearman")
cor.test(ribo_rna_prot_12878$gm12878_rna_mean, log10(as.numeric(ribo_rna_prot_12878$USE)+1), method="spearman")

cor.test(ribo_rna_prot_12878$gm12878_ribo_mean, log10(as.numeric(ribo_rna_prot_12878$USE)+1))
cor.test(ribo_rna_prot_12878$gm12878_rna_mean, log10(as.numeric(ribo_rna_prot_12878$USE)+1))

# Ribo -- MEAN
# USE -> 0.485/0.35, USE.1 0.471/0.37, MAYBE.USE 0.461/0.382, MAYBE.USE2 0.485/0.35
# RNA -> Ever so slightly lower in spearman, equal in pearson

#### CMPARISION WITH ACROSS GENE QUANTIFICATION FROM SILAC
grand_mean_rna <- apply (v3$E[,type=="RNA"], 1, median)
grand_mean_rna  <- data.frame(HGNC=joint_count_ids, grand_mean_rna)
grand_mean_ribo <- apply(v3$E[,type=="Ribo"], 1, median)
grand_mean_ribo <- data.frame (HGNC=joint_count_ids, grand_mean_ribo)
CDS_Lens <- data.frame(HGNC=CDS_IDs, CDS_Len[,1])
merge_ribo_prot <- merge(grand_mean_ribo,protein_absolute_ibaq, by="HGNC" )
merge_ribo_rna_prot <- merge (merge_ribo_prot, grand_mean_rna, by="HGNC")
merge_ribo_rna_prot_len <- merge(merge_ribo_rna_prot, CDS_Lens, by="HGNC")
dim(merge_ribo_rna_prot)
rna_cor <- cor.test(merge_ribo_rna_prot$grand_mean_rna, log10(merge_ribo_rna_prot$ibaq.human))
ribo_cor <- cor.test(merge_ribo_rna_prot$grand_mean_ribo, log10(merge_ribo_rna_prot$ibaq.human))

#



#### ANALYSES BASED ON JUST RIBOSOME PROFILING
### KOZAK SEQUENCE ANALYSIS
kozak_scores <- read.table('~/project/CORE_DATAFILES/Kozak_Reference_Sequence_Scores.txt')
kozak_score_variants <- read.table('~/project/CORE_DATAFILES/Kozak_Variant_Sequence_Scores.txt',
stringsAsFactors=FALSE, fill=T, col.names=paste ("V", seq(1:59), sep=""))
kozak_score_variants_hapmap <- read.table('~/project/CORE_DATAFILES/Kozak_Variant_Sequence_Scores_HapMap.txt', stringsAsFactors=F, fill=T, col.names=paste ("V", seq(1:9), sep=""))
all_kozak_score_variants <- merge(kozak_score_variants, kozak_score_variants_hapmap, all=T, by=c("V1", "V2", "V3") )
all_kozak_score_variants[is.na(all_kozak_score_variants)] <- ""
kozak_var_ind <- all_kozak_score_variants[,-c(2,3)]
number_alleles <- apply (kozak_var_ind, 1, function(x){length(grep('NA', x))})
# Note that there is a strong enrichment for single individual variants

# There is a correlation between number of alleles and difference. 
# If there are a lot of alleles than the difference is less likely to be negative
# If number of alleles is less than 8
kozak_score_variants<- all_kozak_score_variants[,c(1,3)]
kozak_merge <- merge(kozak_score_variants, kozak_scores, by="V1")
p1 <- hist(kozak_merge$V2,50)
p2 <- hist(kozak_merge$V3[number_alleles<10],50)
p2 <- hist(kozak_merge$V3,50)
plot(p1, col=rgb(0,0,1,1/4), xlim=c(-16,-6))
plot(p2, col=rgb(1,0,0,1/4), xlim=c(-16,-6), add=T)
wilcox.test(kozak_merge$V3, kozak_merge$V2)

multi <- duplicated(kozak_merge[,1]) | duplicated(kozak_merge[,1], fromLast=T)

#Go over the variant containing transcripts
# Some transcripts have multiple variants and should be treated separately
# Then for each transcript calculate difference in ribo-expression 
# Each transcript will contribute a delta Kozak and delta expression
kozak_diff <- kozak_merge$V2[!multi] - kozak_merge$V3[!multi]
kozak_var_ind <- kozak_var_ind[!multi,]
# grep ("10847|19240", sample_labels) gives the index of the expr value
# Go over the individuals, grep the samples and calculate difference in mean
ribo_diff <- c()
list_of_pval <- c()
for (i in 1:length(kozak_var_ind[,1])) { 
  ind_unique <- unique(grep("NA",kozak_var_ind[i,-1], value=T))
  ind_unique <- sub("NA", "GM", ind_unique)
  ribo_index <- grep(paste(ind_unique , collapse="|"), sample_labels)
  # CHECK ENST EQUIVALENT IS PRESENT IN V$E, IF NOT ADD NA
  # HERE WE CAN DO MORE WITH THE STATS --WEIGHTED MEAN, ETC
  my_index <- which ( CDS_IDs[isexpr] == enst_hgnc[grep (kozak_var_ind[i,1], enst_hgnc),2] )
  if (length(my_index) & length(ribo_index) %% 50 != 0 & length(ind_unique) > 1 & length(ind_unique) < 29) {
    ribo_diff <- c(ribo_diff, weighted.mean(v$E[my_index,ribo_index],v$weights[my_index,ribo_index] ) - weighted.mean(v$E[my_index, -ribo_index], v$weights[my_index, -ribo_index] ))
    index_factor <- rep(0,times=50)
    index_factor[ribo_index] <- 1
    list_of_pval <- c(list_of_pval, summary(lm(v$E[my_index,] ~ as.factor(index_factor), weights=v$weights[my_index,]))$coefficients[2,4])
  }
  else { 
    ribo_diff <- c(ribo_diff,NA)
    list_of_pval <- c(list_of_pval,NA)
  }
}
plot(ribo_diff, kozak_diff, cex=0.5, pch=19)
abline(v=c(0,0.5,-0.5), h=c(0,0.5, -0.5))
# fmat <- matrix(nrow=2, ncol=2)
fmat[1,] <- c(0,7)
# fmat[2,] <- c(13,13)
# fisher.test(fmat)

# When we look at the pvalues, All the significant changes have Kozak strength changes in the upper half
# The direction of the effect is inconsistent in 3/8. 2/3 are mitochondrial ribosomal proteins
# Look at a browser shot that is an average of the individuals; Look at population frequency of the variants
c1 <- p.adjust(list_of_pval)
q1 <- !is.na(c1)
q2 <- !is.na(c1) & c1 > 0 & c1 < 1
q2 <- !is.na(c1) & c1 > 0 & c1 < 0.2
quantile((kozak_diff)[q1], na.rm=T)
length((kozak_diff)[q1])
length((kozak_diff)[q2])
length(which (kozak_diff[q1] > log(2)) )
fmat[2,] <- c(5, 78)
fmat[1,] <- c(1, 277-83-1)
fisher.test(fmat)
boxplot(abs(kozak_diff[q2]), abs(kozak_diff[!q2]))
kozak_var_ind[q2,]
plot(ribo_diff[q2], kozak_diff[q2], cex=0.5, pch=19)
abline(v=c(0,0.5,-0.5), h=c(0,0.5, -0.5))

# Look at the effect of number of alleles in the significant ones
# Look at the multiple ones; the analysis of multiple ones could be similar to number of alleles
# Multiple variants always occur on different positions
# It seems like that all the significant ones are high allele frequency. We might exclude singletons from testing

# Allelic effect is interesting. For ENST00000245539.6 
my_index <- 5010
index_factor[c(1,2,3,4,5,6,7,8,9,15,22,23,26,27,28,30,31,32)] <- 2
boxplot(v$E[my_index,]~index_factor)
summary(lm(v$E[my_index,] ~ as.factor(index_factor), weights=v$weights[my_index,]))
# Two copies is better than one
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


