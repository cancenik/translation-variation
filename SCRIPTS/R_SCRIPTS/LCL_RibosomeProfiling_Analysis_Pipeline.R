library("edgeR")
library("hexbin")
library("limma")
library("qtl")
library ("sva")
library("MASS")
library("kohonen")
library("pgirmess")
library("RColorBrewer")
library("plyr")
library("apcluster")
library("RDAVIDWebService")
library("gplots")
# Use for progress bars
#library("tcltk")
source('~/project/kohonen2/R/plot.kohonen.R')
#library(nlme)
library("lme4")
library("RLRsim")
library ("VennDiagram")
library("hash")

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
select_first <- function(s) strsplit(s, "_")[[1]][1]
labels_ribo = sapply(colnames(v3$E[,85:134][,reps_ribo_joint]),select_first)
## SAVE FROM ZOOM
plot(joint_ribo_rep, labels= labels_ribo, sub = "", xlab="", main= "Clustering Ribosome Profiling Libraries")
reps_rna_joint <- duplicated(sample_labels_joint[1:84]) | duplicated(sample_labels_joint[1:84], fromLast = TRUE)
joint_dd_rna <- dist(t (v3$E[,1:84][,reps_rna_joint]))
joint_rna_hc <- hclust(joint_dd_rna)
labels_rna = sapply(colnames(v3$E[,1:84][,reps_rna_joint]) , select_first)
labels_lab_rna = sapply(colnames(v3$E[,1:84][,reps_rna_joint]), function(s) {tail(strsplit(s, "_")[[1]],n=1 )} )
plot(joint_rna_hc, labels= paste(labels_rna, labels_lab_rna, sep="_"), sub = "", xlab="", main= "Clustering RNA-Seq Libraries")
plotMDS(v3$E, labels=type)
# plotMDS (v3$E[,type=="RNA"])

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
# Test for difference in variance per gene across individuals
# Partition sum of squares to estimate between individual variance take out variance component due to rep to rep variance
# Calculate F-test p-value and compare joint significance or separate significance
sample_id_all <- paste(sample_labels_joint_common, type_common, sep = "_")
rna_ribo_mean_diff <- apply(joint_expression_common$E, 1, function(x) {mean(x[1:49] - mean(x[50:80]))} ) 


# Alternative approach is to use a random effects model
# Test for significance of the random effects with ‘RLRsim’ package

rna_cols_to_select = 1:49
ribo_cols_to_select = 50:80

# One control for this analysis is to use the subset of the data generated by Snyder
# The concern is amplified variation between individuals due to multiple sources.
#####
rna_cols_to_select = 1:36
ribo_cols_to_select = c(50:57,62:67, 74:80)
#####
random_effect_stat_rna = c()
random_effect_stat_ribo = c()
random_effect_p_val_rna = c()
random_effect_p_val_ribo = c()
# The p-value is fragile for low number of similuations so increased to 500k
# This section is quite computationally intensive
# We estimate that this will take ~12h

# # THIS WAS RUN ONCE AND RESULTS STORED
for ( i in 1:nrow(joint_expression_common$E)) { 
  # Summary keeps sum of squares as $ Sum Sq : num  2.82 3.4
  mA_ribo = lmer(joint_expression_common$E[i,ribo_cols_to_select] ~ 1 + (1| as.factor(sample_id_all[ribo_cols_to_select])), 
                         weights= joint_expression_common$weights[i,ribo_cols_to_select], REML=F)
  m0_ribo = lm(joint_expression_common$E[i,ribo_cols_to_select] ~ 1 , 
                 weights= joint_expression_common$weights[i,ribo_cols_to_select] )
  if(getME(mA_ribo, "theta") < 1e-5 ) {
    random_effect_p_val_ribo = c(random_effect_p_val_ribo,1)
    random_effect_stat_ribo = c(random_effect_stat_ribo, 0)
  }
  else {
  tmp_ribo = exactLRT(m = mA_ribo, m0 = m0_ribo, nsim=500000)
  random_effect_p_val_ribo = c(random_effect_p_val_ribo, tmp_ribo$p.value)
  random_effect_stat_ribo = c(random_effect_stat_ribo, tmp_ribo$statistic)
  }
  mA_rna =  lmer(joint_expression_common$E[i,rna_cols_to_select] ~ 1 + (1| as.factor(sample_id_all[rna_cols_to_select])) , 
                weights= joint_expression_common$weights[i,rna_cols_to_select], REML=F )  
  m0_rna =  lm(joint_expression_common$E[i,rna_cols_to_select] ~ 1 , 
                weights= joint_expression_common$weights[i,rna_cols_to_select] )  
  if(getME(mA_rna, "theta") < 1e-5 ) {
    random_effect_p_val_rna = c(random_effect_p_val_rna, 1)
    random_effect_stat_rna = c(random_effect_stat_rna, 0)      
  }
  else {
  tmp_rna = exactLRT(m = mA_rna, m0 = m0_rna, nsim=500000)
  random_effect_p_val_rna = c(random_effect_p_val_rna, tmp_rna$p.value)
  random_effect_stat_rna = c(random_effect_stat_rna, tmp_rna$statistic)      
  }
}
# save (random_effect_p_val_rna,random_effect_stat_rna, 
#       random_effect_stat_ribo, random_effect_p_val_ribo, 
#       file= paste(data_dir , "Random_Effect_Model_stats_RiboRNA", sep = "" ) )
load(file= paste(data_dir , "Random_Effect_Model_stats_RiboRNA", sep = "" ) )
hist(p.adjust(random_effect_p_val_rna, method = "holm"))
hist(p.adjust(random_effect_p_val_ribo, method = "holm"))
ribo_sig_random = which (p.adjust(random_effect_p_val_ribo, method = "holm") < 0.05)
rna_sig_random = which (p.adjust(random_effect_p_val_rna, method = "holm") < 0.05)
joint_sig_random = which (p.adjust(random_effect_p_val_rna, method = "holm") < 0.05 & 
                            p.adjust(random_effect_p_val_ribo, method = "holm") < 0.05)
length(ribo_sig_random)
length(rna_sig_random)
length(joint_sig_random)
length (random_effect_p_val_ribo )
significant_difference_counts = matrix (nrow=2, ncol=2)
significant_difference_counts[,1] <- c(length (random_effect_p_val_ribo ) -length(rna_sig_random), length(rna_sig_random) ) 
significant_difference_counts[,2] <- c(length (random_effect_p_val_ribo ) - length(ribo_sig_random), length(ribo_sig_random) )

#pdf ('~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Number_Variable.pdf', width = 5, height = 3)
barplot (significant_difference_counts, ylab = "Number of Genes", names = c("RNA Expression", "Ribosome Occupancy"),
         col = c("Yellow2", "Green2"), main = "", beside=T)
#dev.off()
random_effect_df = data.frame (ID = hgnc_to_ensg_convert(row.names(v3)), 
    RNA_Stat=random_effect_stat_rna, RNA_P = p.adjust(random_effect_p_val_rna, method = "holm"),
    RIBO_Stat=random_effect_stat_ribo, RIBO_P = p.adjust(random_effect_p_val_ribo, method = "holm"))
# write.table(random_effect_df, file = paste(data_dir , "Random_Effect_Model_stats_DF_Table.txt", sep = "" ), row.names=F )

## VENN DIAGRAM REPRESENTATION
## ADD Color
venn.mixed <- venn.diagram(
  x = list (
    Ribosome_Ocuppancy = ribo_sig_random,
    RNA_Expression = rna_sig_random
  ),
  filename = NULL,  fill = c("Yellow2", "Green2")
);

#pdf("~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Variation_Mixed_Model_VennColored.pdf");
grid.draw(venn.mixed);
#dev.off();


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

#topTable(ribo_fit2)
#topTable(rna_fit2)
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

# We can subset the set of transcripts with respect to non-zero inter-individual variation
# Variation in both RNA and Ribo is in joint_sig_random
# Variation in either
sig_ribo_rna_random = union(rna_sig_random, ribo_sig_random)
sig_ribo_rna_random_strict = intersect(rna_sig_random, ribo_sig_random)

linfeng_common <- colnames(linfeng_protein) %in% unique(sample_labels_joint)
linfeng_protein_common <- linfeng_protein[,linfeng_common]
linfeng_protein_common <- merge(linfeng_protein_common, ensg_hgnc, by.x="row.names", by.y="ENSG")
# Subset linfeng protein to only no NAs. If we want we can also include samples with missing values sqrt(#individuals)?
# For correlation we can use use="pairwise.complete.obs"
number_NAs <- 1
linfeng_protein_na <- linfeng_protein_common[apply(is.na(linfeng_protein_common), 1, sum) < number_NAs, ]
linfeng_protein_ribo_rna <- merge (v3$E, linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_nonzero_variance <- 
  merge (v3$E[sig_ribo_rna_random,], linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_zero_variance <- 
  merge (v3$E[-sig_ribo_rna_random,], linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_strict_variance <- 
  merge (v3$E[sig_ribo_rna_random_strict,], linfeng_protein_na, by.x="row.names", by.y="HGNC")

## SWITCHED TO USING THE NONZERO VARIANCE SET
# We will add the zero_variance ones on top 
## INSTEAD OF EXTRACTING THE FUNCTION WE WILL SWITCH THE DATASET
# linfeng_protein_ribo_rna = linfeng_protein_ribo_rna_zero_variance
# linfeng_protein_ribo_rna = linfeng_protein_ribo_rna_strict_variance
linfeng_protein_ribo_rna = linfeng_protein_ribo_rna_nonzero_variance

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
# 0.38; Strict 0.65; non-sig .16
median(across_ind_rna_correlation)
# 0.42; Strict 0.67; non-sig .18
ks.test (across_ind_ribo_correlation, across_ind_rna_correlation )
# p= 0.03; Strict p= 0.61
color_by_pval <- rep(0, length(ribo_replicate_mean_prot))
# FDR ~ 20% 0.0004 ; For Non-zero variation %10 FDR is 0.0002
# FDR ~10: for non-sig 0.00002 for rna ; for ribo  0.00005 - At 20% FDR -> 6 significant RNA; 13 sig Ribo
pval_cutoff <- 0.0002
color_by_pval[across_ind_ribo_correlation_pval < pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff] <- 1
color_by_pval[across_ind_ribo_correlation_pval< pval_cutoff & across_ind_rna_correlation_pval >= pval_cutoff] <- 2
color_by_pval[across_ind_ribo_correlation_pval>=pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff] <- 3
plot(across_ind_ribo_correlation, across_ind_rna_correlation, pch=19, cex=.65, tck=.02, col=c("Black", "Red", "Blue", "Gold")[as.factor(color_by_pval)], 
     xlab="Between Individual Ribosome Occupancy-Protein Expression Correlation", ylab = "Between Individual RNA-Protein Expression Correlation", 
     xlim= c(-0.5,.95), ylim = c(-0.5,.95) )
cor.test(across_ind_ribo_correlation, across_ind_rna_correlation)
# Cor 0.61; p < 2.2e-16; Doesn't change for strict
# Strict sense FDR ~25% p=0.013 from rna. 
# length(color_by_pval) = 87
# 8/87 = 9% non-significant
#### ADD FISHER'S TEST FOR ENRICHMENT -- Huge Enrichment for Both correlating significantly
fmat <- matrix(nrow=2, ncol=2)
fmat[1,1] = length(which(across_ind_ribo_correlation_pval < pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff))
fmat[1,2] =  length(which(across_ind_ribo_correlation_pval < pval_cutoff & across_ind_rna_correlation_pval > pval_cutoff))
fmat[2,1] = length(which(across_ind_ribo_correlation_pval > pval_cutoff & across_ind_rna_correlation_pval < pval_cutoff))
fmat[2,2] = length(which(across_ind_ribo_correlation_pval > pval_cutoff & across_ind_rna_correlation_pval > pval_cutoff))
fisher.test(fmat)
# Odds ratio 14.7; p< 2.2e-16 ; Doesn't change for strict 

# Remove one column which is not shared
ribo_replicate_mean_rna <- lapply(ribo_replicate_mean_prot, function(x){x <- x[-24,]})
c2 <- mapply (cbind, rna_replicate_mean_prot, ribo_replicate_mean_rna, SIMPLIFY=F)
across_ind_rna_ribo <- as.numeric(lapply(c2, function(x){ cor(x[,2], x[,4],method="spearman") }))
# Histograms of Across Ind Ribo-Prot, RNA-Prot and Ribo-RNA correlations
p1 <- hist(across_ind_ribo_correlation,40)
p2 <- hist(across_ind_rna_correlation,40)
p3 <- hist(across_ind_rna_ribo, 40)

# Calculate Density
p1$counts = p1$counts / sum(p1$counts) 
p2$counts = p2$counts / sum(p2$counts) 
p3$counts = p3$counts / sum(p3$counts) 

# Make replicate variables
# ribo_nonzero = p1
# rna_nonzero = p2
# rr_nonzero = p3
# 
# ribo_strict = p1
# rna_strict = p2
# rr_strict = p3
# 
# ribo_zero = p1
# rna_zero = p2
# rr_zero = p3

#pdf(file = "~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Across_Individual_Correlations_ALLSUPER.pdf", width=5, height=9)
par(mfrow = c(3, 1))
# USE rgb(204/255,204/255,204/255,1/8) for no-variance
# COL opaque for strict ; alpha 1/4 COL for union
plot(ribo_zero, col=rgb(204/255,204/255,204/255,1/8), xlim=c(-1,1), ylim=c(0,.15), xlab="Spearman Correlation Coefficient", main="Ribosome Occupancy-Protein Level")
#plot(ribo_nonzero, col=rgb(0,0,1,1/4), xlim=c(-1,1), ylim=c(0,.15), add=T)
plot(ribo_strict, col=rgb(0,0,1,1/2), xlim=c(-1,1), ylim=c(0,.15), add=T)

plot(rna_zero, col=rgb(204/255,204/255,204/255,1/8), xlim=c(-1,1), ylim=c(0,.15), xlab="Spearman Correlation Coefficient", main="RNA Expression-Protein Level")
#plot(rna_nonzero, col=rgb(1,0,,1/4), xlim=c(-1,1), ylim=c(0,.15), add=T)
plot(rna_strict, col=rgb(1,0,0,1/2), xlim=c(-1,1), ylim=c(0,.15), add=T)

plot(rr_zero, col=rgb(204/255,204/255,204/255,1/8), xlim=c(-1,1), ylim=c(0,.15), xlab="Spearman Correlation Coefficient", main="RNA Expression-Ribosome Occupancy")
#plot(rr_nonzero, col=rgb(0,1,0,1/4), xlim=c(-1,1), ylim=c(0,.15), add=T)
plot(rr_strict, col=rgb(0,1,0,1/2), xlim=c(-1,1), ylim=c(0,.15), add=T)

#dev.off()

# Try adding separate histograms -- We can plot all three on the same with different transparency
# For strict and variable; Non-variable can be added as grayscale. We need to plot as Density
#pdf(file = "~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Across_Individual_Correlations.pdf", width=9, height=6.5)
#pdf(file = "~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Across_Individual_Correlations_ZeroVariance.pdf", width=4, height=11)
#pdf(file = "~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Across_Individual_Correlations_Strict.pdf", width=4, height=11)
#pdf(file = "~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/Across_Individual_Correlations_Separate.pdf", width=4, height=11)
par(mfrow = c(3, 1))
plot(p1, col=rgb(0,0,1,1/4), xlim=c(-1,1), ylim=c(0,50), xlab="Spearman Correlation Coefficient", main="Ribosome Occupancy-Protein Level")
plot(p2, col=rgb(1,0,0,1/4), xlim=c(-1,1), ylim=c(0,50), xlab="Spearman Correlation Coefficient", main="RNA Expression-Protein Level")
#     add=T)
plot(p3, col=rgb(0,1,0,1/4), xlim=c(-1,1),ylim=c(0,50), xlab="Spearman Correlation Coefficient", main="RNA Expression-Ribosome Occupancy")
#     add=T)
#dev.off()
legend(.4,170,c("RNA Expression-\nProtein Expression", "Ribosome Occupancy-\nProtein Expression", "RNA Expression-\nRibosome Occupancy"), 
       yjust =0.5, x.intersp=0.2, y.intersp=1.5,bty="n", border="white", fill=c(rgb(1,0,0,1/4), rgb(0,0,1,1/4), rgb(0,1,0,1/4)), cex=.9)
#dev.off()
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
# USE i = 18 ; 20 to create example barplots showing the median to indicate what goes into SOM
exemp1 = 18
exemp2 = 20

pdf('~/Google_Drive/Manuscript Figures/Across_Gene_Comparison/Example_Genes.pdf', height=12, width=4)
par (mfrow = c(6,1), las = 1)

rna_cols = rep ("gray", length(v3$E[exemp1,type=="RNA"]))
ribo_cols = rep ("gray", length(v3$E[exemp1,type=="Ribo"]))
te_cols = rep ("gray", length(te_fit4$coefficients[exemp1,]))
rna_cols[which(v3$E[exemp1,type=="RNA"] == median(v3$E[exemp1,type=="RNA"]))] = "red"
ribo_cols[which.min(abs (v3$E[exemp1,type=="Ribo"] - median(v3$E[exemp1,type=="Ribo"] )))] = "red"
te_cols [which.min(abs (te_fit4$coefficients[exemp1,] - median(te_fit4$coefficients[exemp1,]))) ] = "red"
barplot(v3$E[exemp1,type=="RNA"], names.arg = F , col=rna_cols )
barplot(v3$E[exemp1,type=="Ribo"] , names.arg = F , col=ribo_cols)
barplot (te_fit4$coefficients[exemp1,], names.arg = F , col=te_cols)

rna_cols = rep ("gray", length(v3$E[exemp2,type=="RNA"]))
ribo_cols = rep ("gray", length(v3$E[exemp2,type=="Ribo"]))
te_cols = rep ("gray", length(te_fit4$coefficients[exemp2,]))
rna_cols[which(v3$E[exemp2,type=="RNA"] == median(v3$E[exemp2,type=="RNA"]))] = "red"
ribo_cols[which.min(abs (v3$E[exemp2,type=="Ribo"] - median(v3$E[exemp2,type=="Ribo"] )))] = "red"
te_cols [which.min(abs (te_fit4$coefficients[exemp2,] - median(te_fit4$coefficients[exemp2,]))) ] = "red"

barplot(v3$E[exemp2,type=="RNA"], names.arg = F , col=rna_cols)
barplot(v3$E[exemp2,type=="Ribo"], names.arg = F , col=ribo_cols)
barplot (te_fit4$coefficients[exemp2,], names.arg = F , col=te_cols)

dev.off()

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
ribo_rna_te_prot$ibaq.human <- log10(ribo_rna_te_prot$ibaq.human)
all_cors = cor(ribo_rna_te_prot[,-1], use="complete.obs", method="spearman")
plot(ribo_rna_te_prot$grand_mean_te, ribo_rna_te_prot$ibaq.human, pch= 19, cex =.4, tck = .02, xlim = c(-3,3), xlab="Median Translation Efficiency", ylab="log10 iBAQ protein expression")
text(par("usr")[2] - 0.75, par("usr")[4] - 0.75, labels= paste("Spearman rho", signif(all_cors[3,4], 2) , sep=" = "))
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
text(par("usr")[2]-0.5, par("usr")[4]-0.5, labels=paste ("rho=", round(rna_cor$estimate,2) , sep="=") , adj=c(1,1), cex=2)
plot(merge_ribo_rna_prot$grand_mean_ribo, log10(merge_ribo_rna_prot$ibaq.human), xlab="Ribosome Profiling Expression", ylab="log10 iBAQ protein expression", pch=19, cex=.4)
text(par("usr")[2]-0.5, par("usr")[4]-0.5, labels=paste ("rho= ", round(ribo_cor$estimate,2) , sep="=") , adj=c(1,1), cex=2)
plot(merge_ribo_rna_prot$grand_mean_ribo, merge_ribo_rna_prot$grand_mean_rna, xlab="Ribosome Profiling Expression", ylab="RNA Expression", pch=19, cex=.4)
text(par("usr")[2]-0.5, par("usr")[4]-0.5, labels=paste ("rho= ", round(ribo_rna$estimate,2) , sep="=") , adj=c(1,1), cex=2)
#

### Implement permutation scheme to test significance of difference in correlation
true_diff = ribo_cor$estimate - rna_cor$estimate

swap <- function(x) {
  if (runif(1) > 0.5) {
    return (c(x[2],x[1]))
  }
  else {
    return(c(x[1],x[2]))
  }
}

cor_dif <- c()
for (i in 1:10000) {
  perm <- apply (merge_ribo_rna_prot[,c(2,8)], 1, swap)
  cor_dif <- c(cor_dif,
        cor.test(perm[1,], log10(merge_ribo_rna_prot$ibaq.human), method="spearman")$estimate - cor.test(perm[2,], log10(merge_ribo_rna_prot$ibaq.human), method="spearman")$estimate)
}
length(which(abs(cor_dif) > true_diff))
#pdf ( "~/Google_Drive/Manuscript Figures/Across_Gene_Comparison/PermutationCorrelationDifference.pdf", width=5, height=5)
par(las=1)
hist (abs(cor_dif), xlim = c(0, 0.12), xlab= "Absolute Difference in Correlation Coefficient", main= "", breaks = 50, tck=.02)
abline(v= true_diff, col = "red")
#dev.off()

## APPLY SOMs to different versions of data. 
# APPLY BDK to ribo_rna_te_prot
# We can do a bidirectional kohonen predicting protein levels and plot the cor per cell with either or all
# We can cluster cells by which parameter correlates best with protein levels in each cell. 
redblue_cols <- function(x) {brewer.pal(x, "RdBu")}

total_cells_bdk <- floor(sqrt(dim(ribo_rna_te_prot)[2]/2) * sqrt (dim(ribo_rna_te_prot)[1]))
if (floor(sqrt(total_cells_bdk/1.3333)) %% 2 == 0) { 
  ydim = floor(sqrt(total_cells_bdk/1.3333))
} else { 
  ydim = floor(sqrt(total_cells_bdk/1.3333)) + 1
}
xdim = floor(total_cells_bdk/ydim + 0.5)

### WE NEED TO DECIDE WHETHER THIS SCALING IS USEFUL
# We might want to standardize the measurements
# One issue is that TE and RNA-RIBO does not have the same interpretation
# TE is the median of a linear model coefficient. Maybe we should switch to gene specific measure
# Another idea is to change everthing into same scale of quantiles before applying SOM
#ribo_rna_te_prot[,-1] <- scale(ribo_rna_te_prot[,-1], scale=F)
my.ecdf <- function (x) {ecdf(x)(x)}
abs.som.data <- apply(ribo_rna_te_prot, 2 , my.ecdf)
abs.som.data.noNA <- abs.som.data[!apply(is.na(abs.som.data), 1, any),]

### Add heatmap representation of the data to include as a motivator for SOM
#pdf ("~/Google_Drive/Manuscript Figures/Across_Gene_Comparison/SOM_Input.pdf", width=10, height=10)
h2 = heatmap.2 ( abs.som.data.noNA[20:47,], Colv =F, Rowv =F, scale="none", dendrogram = "none", 
      density.info="none", trace = "none",  labRow=paste("Gene",seq(1,28,length.out=28) , sep ="-"), 
      breaks = seq(0,1,length.out=50), col = function(x) {colorpanel(x, 'blue3','white', 'red2')},
      cexRow=.9,cexCol=.5,
)
#dev.off()

# Run the SOM, 1000 times and keep track of the distances pick the one with the min 75% distance
# This section was run once to determine the best seed
mean_distance <- 1
# for (i in 1:500) { 
#   set.seed(i)
#   absolute.som <- som(abs.som.data.noNA, toroidal=T, grid=somgrid(xdim, ydim, "hexagonal"))
#   if (mean(absolute.som$distances) < mean_distance) { 
#     my_seed <- i
#     mean_distance <- mean(absolute.som$distances)
#   }
# }
# best_seed is 137 for xweight .8 with 90% percentile
#set.seed(137)
# This was used before quantization of the data
# absolute.som <- bdk(abs.som.data.noNA[,1:3], Y= abs.som.data.noNA[,4] , toroidal=T, xweight=.5, contin=T, grid=somgrid(xdim, ydim, "hexagonal"))
# There is really no need to have a separate protein level when things are quantized
# Best seed for the som with quantized data is:173
set.seed(173)
absolute.som <- som(abs.som.data.noNA, toroidal=T, grid=somgrid(xdim, ydim, "hexagonal"))
# abs.som.data.noNA -> Numeric Matrix, 4th column is ibaq.human
prot_mean <-  by (data.frame(abs.som.data.noNA[,4]), absolute.som$unit.classif, FUN = colMeans)
ribo_mean <- by (data.frame(abs.som.data.noNA[,1]), absolute.som$unit.classif, FUN = colMeans)
# som.exp$unit.classif has the info about where each gene went
# We also need some visually appealing colors
## Write function to apply function to each cell of SOM.
# som.exp$unit.classif keeps track of the cell each object belongs to
# Left bottom is 1, the row above that is 1+xdim 
# Need vector of length xdim * ydim

# Unit-wise correlation is less meaningful. Might be better to try correlation after k-means

ribo_prot_cor_across_genes_som <- by (data.frame(abs.som.data.noNA[,c(1,4)]), absolute.som$unit.classif, FUN = function(x){cor(x, method="spearman")[1,2]} )
rna_prot_cor_across_genes_som <- by (data.frame(abs.som.data.noNA[,c(2,4)]), absolute.som$unit.classif, FUN = function(x){cor(x,  method="spearman")[1,2]} )
te_prot_cor_across_genes_som <- by (data.frame(abs.som.data.noNA[,c(3,4)]), absolute.som$unit.classif, FUN = function(x){cor(x)[1,2]} )
plot.kohonen(absolute.som, property=ribo_prot_cor_across_genes_som, type="property", main ="Ribosome Occupancy Protein Correlation", contin=T,zlim=c(-1,1),palette.name=redblue_cols, ncolors=11)
plot.kohonen(absolute.som, property=rna_prot_cor_across_genes_som, type="property", main = "RNA Expression Protein Correlation", contin=T,zlim=c(-1,1),palette.name=redblue_cols, ncolors=11)
plot.kohonen(absolute.som, property=te_prot_cor_across_genes_som, type="property", main="Translation Efficiency Protein Correlation",contin=T,zlim=c(-1,1), palette.name=redblue_cols, ncolors=11)

#pdf('~/Google_Drive/Manuscript Figures/Across_Gene_Comparison/Absolute_SOM_Codes.pdf', width=7, height=5)
par(mfrow = c(2, 2))
plot.kohonen(absolute.som, type = "property", property=absolute.som$codes[,1],palette.name=function(x){colorpanel(x, 'blue3','white', 'red2')}, ncolors=50, contin=T, main="Ribosome Occupancy", zlim = c(0,1))
plot.kohonen(absolute.som, type = "property", property=absolute.som$codes[,2],palette.name=function(x){colorpanel(x, 'blue3','white', 'red2')}, ncolors=50, contin=T, main="RNA Expression" , zlim = c(0,1))
plot.kohonen(absolute.som, type = "property", property=absolute.som$codes[,3],palette.name=function(x){colorpanel(x, 'blue3','white', 'red2')}, ncolors=50, contin=T, main="Translation Efficiency", zlim = c(0,1))
plot.kohonen(absolute.som, property=absolute.som$codes[,4], type = "property", palette.name=function(x){colorpanel(x, 'blue3','white', 'red2')}, ncolors=50, contin=T, main="Protein Level", zlim = c(0,1))
#dev.off()

plot.kohonen(absolute.som, property=prot_mean, type = "property", palette.name=redblue_cols, ncolors=11)

#plot.kohonen(absolute.som, type="classes", property=absolute.som$codes[,1:3], scale=T,  palette.name=function(x) {brewer.pal(x,"Dark2")}, bgcol= brewer.pal(7,"Set2")[cluster.membership])
#abs.som.which.max <- apply(absolute.som$codes[,1:3], 1, which.max)
#plot.kohonen(absolute.som, type = "property", property=abs.som.which.max,palette.name=redblue_cols, ncolors=3, contin=F, main="Which.Max" )

# We can add cluster boundaries by k-means
# kmeans(absolute.som$codes, 12)$cluster
#add.cluster.boundaries(absolute.som, abs.som.which.max)
#add.cluster.boundaries(absolute.som,kmeans(absolute.som$codes, 3)$cluster)

# Affinity propagation
ap.cluster <- apcluster(negDistMat(r=2), absolute.som$codes, q=0.25)
# ap.cluster@clusters == List of Clusters
cluster.membership <- c()
for (i in 1:length(ap.cluster@clusters)) { 
cluster.membership[ap.cluster@clusters[[i]]] <- i
}
# We can plot a pie-chart version. Convert codebook ribo, rna, te to quantile
# Use plot.kohonen(absolute.som, type="classes", property=numericmatrix(absolute.som$codes))
# Update function to introduce a more meaningful scaling for the pie chart
# Update function to change sum(codes) > 1 to > 0

# Change colors for inside the piechars and cluster colors so they separate out much more nicely
colnames(absolute.som$codes) <- c("Ribosome Occupancy", "RNA Expression", "Translation Efficiency", "Protein Level")
absolute.som$codes <- absolute.som$codes[,c(2,1,3,4)]
#pdf(file= "~/Google_Drive/Manuscript Figures/Across_Gene_Comparison/SOM_All_ClusterColored_New.pdf", width=12, height=8)
plot.kohonen(absolute.som, type="classes", property=absolute.som$codes[,c(1,3)], scale=T, 
palette.name=function(x) {c("Yellow2", "Green2")}, bgcol= gray.colors(length(ap.cluster@clusters), start=0, end=1)[cluster.membership])
# ADD CLUSTER BOUNDARIES FOR GO ENRICHED CATEGORIES
#dev.off()

## A version of the cluster membership that represents clusters only with different contrasting colors
#pdf(file= "~/Google_Drive/Manuscript Figures/Across_Gene_Comparison/SOM_All_ClusterOnlyColors.pdf", width=12, height=8)
plot.kohonen(absolute.som, type= "property", property= cluster.membership, scale = T, 
             palette.name = function(x){return (c("grey","black","brown","red","blue" ,"green","cyan","magenta","yellow" ))})
#dev.off()
## HEATMAP OF CLUSTER EXPRESSION
#pdf(file= "~/Google_Drive/Manuscript Figures/Across_Gene_Comparison/SOM_All_ClusterExpressionHeatMap.pdf", width=8, height=10)
h1 = heatmap.2 (absolute.som$codes[ap.cluster@exemplars, ], Rowv =F, scale="none", dendrogram = "none", 
density.info="none", trace = "none", keysize = 1, labRow = paste( "Cluster", seq(1, 9, by = 1), sep = " "),
breaks = 50, col = function(x) {colorpanel(x, 'blue3','white', 'red2')},cexRow=.9,cexCol=.5,
RowSideColors = c("grey","black","brown","red","blue" ,"green","cyan","magenta","yellow" ) )
#dev.off()
##

absolute.som$codes[ap.cluster@exemplars, ]
aggregate(absolute.som$codes , by = list(cluster.membership), FUN=mean)
#abs.som.which.max <- apply(absolute.som$codes[,1:3], 1, which.max)
#plot.kohonen(absolute.som, type = "property", property=abs.som.which.max,palette.name=redblue_cols, ncolors=3, contin=F, main="Which.Max" )

cluster.membership5 = cluster.membership8 = cluster.membership
cluster.membership5[cluster.membership5 != 5] = 1
cluster.membership8[cluster.membership8 != 8] = 1
add.cluster.boundaries(absolute.som, cluster.membership5, col = "blue")
add.cluster.boundaries(absolute.som, cluster.membership8, col = "red")

# corSimMat(method="spearman")
plot(ap.cluster, absolute.som$codes)
heatmap(ap.cluster)

#ap.cluster.full <- apcluster(negDistMat(r=2), absolute.som$data,q= 0)
#plot(ap.cluster.full, absolute.som$data)
#heatmap(ap.cluster.full)


# superSOM Data Structure is a list of matrices including Linfeng Proteomics with NAs -> This will use individuals in the plot.
# We can also do a general one with just absolute ribo,rna, te, and absolute SILAC amounts (Here SILAC can be coded as a classification parameter for BDF)
# Following Xie, Boyle, et al. The grid is hexagonal, toroid
# Total number of cells, sqrt(DATA_TYPE/2) x sqrt (genes x individuals )

# ACROSS INDIVIDUAL SOM SHOULD FOLLOW ACROSS INDIVIDUAL CORRELATION DATA
ribo_replicate_mean_rna_prot <- unlist(lapply(ribo_replicate_mean_prot, function(x){x <- x[-c(3,17,24),-1]}) )
rna_replicate_mean_prot_ribo <- unlist(lapply(rna_replicate_mean_prot, function(x){x <- x[-c(3,17),-1]}))
dim(ribo_replicate_mean_rna_prot) = c(27, 562)
dim(rna_replicate_mean_prot_ribo) = c(27, 562)
ribo_replicate_mean_rna_prot = t(ribo_replicate_mean_rna_prot)
rna_replicate_mean_prot_ribo = t(rna_replicate_mean_prot_ribo)
col_order = rna_replicate_mean_prot[[1]]$Group.1[-c(3,17)]
sample_labels_joint_prot[type_prot=="Prot"]
# Extra Colum is GM19139 == 22
prot_exp_rna_ribo = linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22]
prot_exp_rna_ribo= matrix(as.numeric(as.matrix(prot_exp_rna_ribo)), ncol = 27)
across_individual_supersom_nonzero_variance_data = list(
  rna = matrix(ecdf(rna_replicate_mean_prot_ribo)(rna_replicate_mean_prot_ribo), ncol=27),
  ribo= matrix(ecdf(ribo_replicate_mean_rna_prot)(ribo_replicate_mean_rna_prot), ncol=27),
  prot = matrix(ecdf(prot_exp_rna_ribo)(prot_exp_rna_ribo), ncol=27)
  )
row.names(across_individual_supersom_nonzero_variance_data$ribo) = 
row.names(linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22])
total_cells.nonzerovar_individuals <- 
floor(sqrt(length(across_individual_supersom_nonzero_variance_data)/2) * 
      sqrt (dim(across_individual_supersom_nonzero_variance_data$rna)[1] * 
      dim(across_individual_supersom_nonzero_variance_data$rna)[2]))

if (floor(sqrt(total_cells.nonzerovar_individuals/1.3333)) %% 2 == 0) { 
  ydim.ind = floor(sqrt(total_cells.nonzerovar_individuals/1.3333))
} else { 
  ydim.ind = floor(sqrt(total_cells.nonzerovar_individuals/1.3333)) + 1
}
xdim.ind = floor(total_cells.nonzerovar_individuals/ydim.ind + 0.5)

# mean_distance <- 1
# for (i in 1:100) { 
#    set.seed(i)
#    supersom.nonzero.ind = supersom (data = across_individual_supersom_nonzero_variance_data, 
#                                     grid=somgrid ( xdim.ind, ydim.ind, "hexagonal"), 
#                                     toroidal=T, contin = T)
#    
#   if (mean(supersom.nonzero.ind$distances) < mean_distance) { 
#      my_seed <- i
#      mean_distance <- mean(supersom.nonzero.ind$distances)
#    }
# }
set.seed(61)
supersom.nonzero.ind = supersom (data = across_individual_supersom_nonzero_variance_data, 
                                 grid=somgrid ( xdim.ind, ydim.ind, "hexagonal"), 
                                 toroidal=T, contin = T)
plot.kohonen (supersom.nonzero.ind, type = "changes")
########

# Calculate cor.estimate for each level grouped by unit.classif
# X is passed as a dataframe to function
gene_wise_cors <- function (mt1 , mt2) { 
  if ( !identical( nrow(mt1),nrow(mt2)) ) { 
    stop ("Unequal rows")
  } 
  sapply(seq.int(nrow(mt1) ), 
         function(k) {cor (mt1[k,], mt2[k,], use="pairwise.complete.obs", method="spearman") } )
}

# Change pval calculation testing to one-sided "greater
gene_wise_cor_pvals <- function (mt1 , mt2) { 
  if (!identical( nrow(mt1),nrow(mt2)) ) { 
    stop ("Unequal rows")
  } 
  sapply(seq.int(nrow(mt1) ), 
         function(k) { cor.test(mt1[k,], mt2[k,], 
                      use="pairwise.complete.obs", method="spearman", alternative="g")$p.value } )
}

across_ind_ribo_rna  = gene_wise_cors(supersom.nonzero.ind$data$ribo, supersom.nonzero.ind$data$rna)
across_ind_ribo_prot  = gene_wise_cors(supersom.nonzero.ind$data$ribo, supersom.nonzero.ind$data$prot)
across_ind_rna_prot  = gene_wise_cors(supersom.nonzero.ind$data$rna, supersom.nonzero.ind$data$prot)

# across_ind_ribo_rna_pval  = gene_wise_cor_pvals(supersom.nonzero.ind$data$ribo, supersom.nonzero.ind$data$rna)
# across_ind_ribo_prot_pval  = gene_wise_cor_pvals(supersom.nonzero.ind$data$ribo, supersom.nonzero.ind$data$prot)
# across_ind_rna_prot_pval  = gene_wise_cor_pvals(supersom.nonzero.ind$data$rna, supersom.nonzero.ind$data$prot)


#SOM, Unit wise correlations -- We are looking at unit wise median correlation
max_cor_unit = c()
min_cor_pval_unit_type = c()
min_cor_pval_unit = c()
cor_diff = c()
float_comparison_with_epsilon = function (x, y, eps)  { 
  return (x + eps > y & x-eps < y)
}
# Adapted from the wiki page
Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- 1-pchisq(Xsq, df = 2*length(p))
  return(p.val)
}
Correlation_Meta <- function (v) { 
 z.r = log((1+v) / (1-v) )  / 2
 ave.z.r = mean (z.r)
 return ( (exp(2*ave.z.r) -1) / (exp(2*ave.z.r) +1) )
}
# First determine the min-pval then update the rest
# Finding the Fisher aggregate p-val and then deciding on that has the problem that 
# some negative correlations are highly significant and skew the highest correlation

## NOT CLEAR WHAT to DO WITH P_VALs
## THE MEDIAN OF PVAL DOESNT MAKE SENSE. 
# Current strategy is to calculate pvals first
# Then we find the values closes to the median pval 
# Then we find the corresponding units and take the median of their correlations
# Maybe better is to combine p-values with Fisher's method?
tol = 1e-5
for ( i in 1: (supersom.nonzero.ind$grid$xdim * supersom.nonzero.ind$grid$ydim)) { 
  if ( identical(nrow(supersom.nonzero.ind$data$ribo[supersom.nonzero.ind$unit.classif == i, ]), NULL ) ) { 
    max_cor_unit = c(max_cor_unit, 0)
    min_cor_pval_unit_type = c(min_cor_pval_unit_type,  -1) 
    min_cor_pval_unit = c(min_cor_pval_unit , 1)
    cor_diff = c(cor_diff , 0 ) 
    next
  }
  rna_prot_unit = gene_wise_cors( supersom.nonzero.ind$data$rna[supersom.nonzero.ind$unit.classif == i, ], 
                                  supersom.nonzero.ind$data$prot[supersom.nonzero.ind$unit.classif == i, ])
  ribo_prot_unit = gene_wise_cors( supersom.nonzero.ind$data$ribo[supersom.nonzero.ind$unit.classif == i, ], 
                                   supersom.nonzero.ind$data$prot[supersom.nonzero.ind$unit.classif == i, ])
  ribo_prot_unit_pval = gene_wise_cor_pvals( supersom.nonzero.ind$data$ribo[supersom.nonzero.ind$unit.classif == i, ], 
                                             supersom.nonzero.ind$data$prot[supersom.nonzero.ind$unit.classif == i, ])
  rna_prot_unit_pval  = gene_wise_cor_pvals( supersom.nonzero.ind$data$rna[supersom.nonzero.ind$unit.classif == i, ], 
                                             supersom.nonzero.ind$data$prot[supersom.nonzero.ind$unit.classif == i, ])
  min_cor_pval_unit = c(min_cor_pval_unit, 
                        min( Fisher.test (rna_prot_unit_pval), Fisher.test (ribo_prot_unit_pval)))
  min_cor_pval_unit_type = c(min_cor_pval_unit_type,
                             which.min( c(Fisher.test (rna_prot_unit_pval), Fisher.test (ribo_prot_unit_pval))))

  # The pval closest to median can be of opposite sign
#   min_cor_pval_unit = c(min_cor_pval_unit, 
#                         min( median(rna_prot_unit_pval), median(ribo_prot_unit_pval),  median(te_prot_unit_pval)))
#   min_cor_pval_unit_type = c(min_cor_pval_unit_type, 
#                              which.min( c(median(rna_prot_unit_pval), median(ribo_prot_unit_pval),  median(te_prot_unit_pval) )))
  if (min_cor_pval_unit_type[i] == 1 ) { 
#     abs.diff = abs(rna_prot_unit_pval -  min_cor_pval_unit[i])
#     closest.pval = which (float_comparison_with_epsilon (abs.diff, min(abs.diff), tol) )
#     max_cor_unit = c(max_cor_unit , median(rna_prot_unit[closest.pval]) )
    max_cor_unit = c(max_cor_unit , Correlation_Meta(rna_prot_unit) )
    cor_diff = c(cor_diff , Correlation_Meta(rna_prot_unit) - Correlation_Meta(ribo_prot_unit) )
  }
  else if (min_cor_pval_unit_type[i] == 2 ) { 
#     abs.diff = abs(ribo_prot_unit_pval -  min_cor_pval_unit[i])    
#     closest.pval = which (float_comparison_with_epsilon (abs.diff, min(abs.diff), tol) )
#     max_cor_unit = c(max_cor_unit , median(ribo_prot_unit[closest.pval]) )
    max_cor_unit = c(max_cor_unit , Correlation_Meta(ribo_prot_unit) )
    cor_diff = c(cor_diff , Correlation_Meta(ribo_prot_unit) - Correlation_Meta(rna_prot_unit) )
    
  }
  else { 
    stop("Incorrect Type")
  }
#  max_cor_unit = c(max_cor_unit , median(rna_prot_unit) )
}
#pdf(file="~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/SuperSOM.pdf", width=4, height=10.5)
#par(mfrow = c(3, 1))
pal <- function (x) {return( c("white", "#B300FF", "#F7FF00"))}
plot.kohonen(supersom.nonzero.ind, property=min_cor_pval_unit_type, type = "property", 
             palette.name=pal, ncolors=3, main = "Highest Correlated Expression Value")
plot.kohonen(supersom.nonzero.ind, property=max_cor_unit, type = "property", 
             palette.name=redblue, ncolors=50, contin=T, 
             zlim=c(-1,1), main = "Spearman Correlation")
plot.kohonen(supersom.nonzero.ind, property=p.adjust(min_cor_pval_unit, method="holm"), type = "property", 
             palette.name=function(x){colorpanel(x,'red', 'white')}, ncolors=50,
             contin=T, main = "Holm's Adjusted Meta P-value")

cor_diff[cor_diff < 0] = 0
# Grey > .35; Can modify by changing bgcolors in property plot
plot.kohonen(supersom.nonzero.ind, property=cor_diff, type = "property", 
             palette.name=function(x){colorpanel(x,'white', 'red')}, ncolors=50, contin=T, 
             zlim=c(0,.35), main = "Spearman Correlation Difference")
#dev.off()

# Plot spearman cor diff > .2 ; adjusted p < .05
high_diff = which(cor_diff > .2)
high_diff = high_diff[p.adjust ( min_cor_pval_unit, method = "holm")[high_diff] < .05]
cor_diff_sigs = rep (0,times=length(cor_diff))
cor_diff_sigs[high_diff] = 1
#pdf(file="~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/SuperSOM_CorDiff.2_Pval05.pdf", width=4, height=3.5)
plot.kohonen(supersom.nonzero.ind, property=cor_diff_sigs, type = "property", 
             palette.name=function(x){colorpanel(x,'white', 'red')}, ncolors=2, 
             main= "Significant Correlation Difference")
#dev.off()

cor_diff_sig_units = which ( as.logical(cor_diff_sigs) == T)
cor_diff_sig_units_rna = cor_diff_sig_units[which (min_cor_pval_unit_type[cor_diff_sig_units] == 1 )]
cor_diff_sig_units_ribo = cor_diff_sig_units[which (min_cor_pval_unit_type[cor_diff_sig_units] == 2 )]
cor_diff_sig_units_rna_genes_index = supersom.nonzero.ind$unit.classif %in% cor_diff_sig_units_rna
cor_diff_sig_units_ribo_genes_index = supersom.nonzero.ind$unit.classif %in% cor_diff_sig_units_ribo

plot.kohonen(supersom.nonzero.ind, type="counts", palette.name=function(x){colorpanel(x, 'white', 'red')}, ncolors=20)
#pdf(file="~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/SuperSOM_Codes.pdf", width=4, height=10.5)
par(mfrow = c(3, 1))
plot.kohonen(supersom.nonzero.ind)
#dev.off()
# Compared to random median correlation per unit is not much higher
# However, there are many more units where the correlations are higher. 
# > quantile(max_cor_unit)
# 0%         25%         50%         75%        100% 
# -0.04895105  0.20279720  0.29195804  0.37762238  0.70279720 
# > quantile(tr)
# 0%         25%         50%         75%        100% 
# 0.006993007 0.209790210 0.307692308 0.461538462 0.842657343 
# True data gives 
# > table(max_cor_unit_type)
# max_cor_unit_type
# 1   2   3 
# 92 148  26 
# # Random Classification gives
# > table(max_cor_unit_type)
# max_cor_unit_type
# 1   2   3 
# 114 143   9 


### GO ENRICHMENT 
# Visualization idea for FuncAssociate GO. 
# Generate term-term kappa statistic similarity matrix. 
# Use cytoscape to visualize and color or size by LOD/p-value
# kappa statistic for GO to Gene matrix

# We use funcassociate to generate results
# We then use kappa_statistic_GO_network.pl to process the results
# The output from R is going to formatted and opened in Cytoscape 
# WE will use the node attributes as the LOD enrichment
# We will define edge similarity with kappa statistic
go_dag_joint = '~/project/CORE_DATAFILES/FUNCASSOCIATE/Mixed_Effect_FuncAssociate/RESULTS/JointIDs_funcassociate_results.tsv_Kappa_Network.sif'
go_dag_ribo = '~/project/CORE_DATAFILES/FUNCASSOCIATE/Mixed_Effect_FuncAssociate/RESULTS/Ribo_Only_funcassociate_results.tsv_Kappa_Network.sif'
go_dag_rna = '~/project/CORE_DATAFILES/FUNCASSOCIATE/Mixed_Effect_FuncAssociate/RESULTS/RNA_Only_funcassociate_results.tsv_Kappa_Network.sif'

go_dag_abs5 = read.table('~/project/CORE_DATAFILES/FUNCASSOCIATE/ABS_SOM/RESULTS/funcassociate_results_Cluster5_3_foldenriched_selected.tsv_Kappa_Network.sif', header = T)
go_dag_abs8 = read.table('~/project/CORE_DATAFILES/FUNCASSOCIATE/ABS_SOM/RESULTS/funcassociate_results_Cluster8_3_foldenriched.tsv_Kappa_Network.sif', header = T)

go_dag = read.table(go_dag_joint, header=T)
go_dag = read.table(go_dag_ribo, header=T)
go_dag = read.table(go_dag_rna, header=T)
go_dag = go_dag_abs5
go_dag = go_dag_abs8

# Calculates kappa similarity between two binary vectors 
calculate_kappa <- function (a1, a2) { 
  Pr_a = sum (!xor(a1,a2)) / length(a1)
  a1_1 = sum(a1) / length(a1)
  a2_1 = sum(a2) /length(a2)
  Pr_e = a1_1 * a2_1 + (1-a1_1) * (1-a2_1)
  kappa = (Pr_a - Pr_e) / (1-Pr_e)
  return (signif(kappa, 2) ) 
}


first_kappas <- c()
for ( i in 1:dim(go_dag)[1]) { 
  for ( j in (i+1):dim(go_dag)[1]) { 
    first_kappas = c(first_kappas, 
                     paste(go_dag[i, 1], go_dag[j, 1],
                           calculate_kappa (go_dag[i, 2:dim(go_dag)[2]], go_dag[j, 2:dim(go_dag)[2]] ), sep="\t"  ) 
    ) 
  }
}
#write(first_kappas, file = paste(go_dag_joint, "modified", sep="_"))
#write(first_kappas, file = paste(go_dag_ribo, "modified", sep="_"))
#write(first_kappas, file = paste(go_dag_rna, "modified", sep="_"))
#write(first_kappas, file = '~/project/CORE_DATAFILES/FUNCASSOCIATE/ABS_SOM/RESULTS/funcassociate_results_Cluster5_3_foldenriched_selected.tsv_Kappa_Network.sif_modified' )

### PATHWAY ENRICHMENT ANALYSIS WITH RDAVID
david<-DAVIDWebService$new(email="ccenik@stanford.edu")
# We will use ENSEMBL_GENE_ID but need to crop the -001 and select species to human
# GENE_ID 92% coverage
# TRANSRIPT_ID 87$ == We will use GENE_IDS

hgnc_to_ensg_convert <- function(x) {
  xdf <- data.frame(ID=x)
  list = unlist(strsplit(merge(xdf, hgnc_to_ensg, by.x="ID", by.y="V5" )$V2, split=".", fixed=T))
  list = list[seq(1, length(list), 2)]
  return (list)
}

filter_by_fdr_fold_enrichment <- function (annot.chart, fdr.threshold, fold.enrichment) { 
 filtered = annot.chart[annot.chart$Fold.Enrichment > fold.enrichment & annot.chart$FDR < fdr.threshold , ]
 return (filtered)
}

# We need a better representation of the results
# We can only set a p-val threshold; update the results with FDR cutoff
# Then we should sort the results by fold enrichment

#setCurrentSpecies(david, )
getAllAnnotationCategoryNames(david)
#getIdTypes(david)
setAnnotationCategories (david, c("GOTERM_CC_ALL", "GOTERM_BP_ALL", "GOTERM_MF_ALL", "KEGG_PATHWAY", "REACTOME_PATHWAY"))

background_list = hgnc_to_ensg_convert(row.names(v3))
addList(david, background_list, idType="ENSEMBL_GENE_ID", listName="V3", listType="Background")

# # Across Ind RNA-Prot; RNA-RIBO Correlation; Nothing significant Commented_out
# across_individual_background_list = hgnc_to_ensg_convert(row.names(linfeng_protein_ribo_rna))
# addList(david, across_individual_background_list, idType="ENSEMBL_GENE_ID", listName="Across_Ind_Background", listType="Background")
# 
# high_ribo_cor_list = hgnc_to_ensg_convert(row.names(linfeng_protein_ribo_rna)[across_ind_ribo_correlation_pval<= pval_cutoff & across_ind_rna_correlation_pval >  pval_cutoff])
# high_rna_cor_list = hgnc_to_ensg_convert(row.names(linfeng_protein_ribo_rna)[across_ind_rna_correlation_pval<= pval_cutoff & across_ind_ribo_correlation_pval >  pval_cutoff ])
# high_ribo_rna_list = hgnc_to_ensg_convert(row.names(linfeng_protein_ribo_rna)[across_ind_rna_correlation_pval<= pval_cutoff & across_ind_ribo_correlation_pval <= pval_cutoff])
# 
# addList(david, high_ribo_cor_list, idType="ENSEMBL_GENE_ID", listName="HighBetIndRiboCorProt", listType="Gene")
# addList(david, high_rna_cor_list, idType="ENSEMBL_GENE_ID", listName="HighBetIndRNACorProt", listType="Gene")
# addList(david, high_ribo_rna_list, idType="ENSEMBL_GENE_ID", listName="HighBetIndRiboRNACorProt", listType="Gene")
# 
# setCurrentBackgroundPosition(david,2)
# setCurrentGeneListPosition(david,1)
# AnnotCluster = getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
# filter_by_fdr_fold_enrichment(AnnotCluster, .1, 2)$Term

# ## Non-zero variance Analysis-- DONE in FuncAssociate in Ordered mode
# ribo_sig_names <- hgnc_to_ensg_convert( row.names(v3)[ribo_sig_random] ) 
# rna_sig_names = hgnc_to_ensg_convert( row.names(v3)[rna_sig_random] )
# joint_sig_names = hgnc_to_ensg_convert( row.names(v3)[joint_sig_random] )
# addList(david, ribo_sig_names, idType="ENSEMBL_GENE_ID", listName="ribo_sig_names", listType="Gene")
# addList(david, rna_sig_names, idType="ENSEMBL_GENE_ID", listName="rna_sig_names", listType="Gene")
# addList(david, joint_sig_names, idType="ENSEMBL_GENE_ID", listName="joint_sig_names", listType="Gene")
# addList(david, setdiff(rna_sig_names, joint_sig_names), idType="ENSEMBL_GENE_ID", listName="rnaonly_sig_names", listType="Gene")
# addList(david, setdiff(ribo_sig_names, joint_sig_names), idType="ENSEMBL_GENE_ID", listName="riboonly_sig_names", listType="Gene")
# setAnnotationCategories (david, c("GOTERM_CC_ALL", "GOTERM_BP_ALL", "GOTERM_MF_ALL", "KEGG_PATHWAY", "REACTOME_PATHWAY"))
# setCurrentBackgroundPosition(david,2)
# setCurrentGeneListPosition(david, 2)
# AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.001, count=2L)
# filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)$Term

## Differential Expression -- RNA/Ribo/TE => This should go with Bilal's RADIAL SETS
# te_fit3 is the across individual difference in translation efficiency, ribo_fit2, rna_fit2 equivalent ones for ribo and rna
# First set is genes differentially expressed at any level
# te.diff.results -- lfc=1
#results.ribo <- decideTests(ribo_fit2, p.value=0.01, lfc=log2(1.5))
#results.rna <- decideTests(rna_fit2, p.value=0.01, lfc=log2(1.5))

### Both up/down have interesting categories for enrichment. We can use the radial sets as a visualization
te_down_list = hgnc_to_ensg_convert(names(apply(te.diff.results == -1 , 1 , any))[apply(te.diff.results == -1 , 1 , any)])
te_up_list = hgnc_to_ensg_convert(names(apply(te.diff.results == 1 , 1 , any))[apply(te.diff.results == 1 , 1 , any)])
te_any  = hgnc_to_ensg_convert(names(apply(te.diff.results != 0 , 1 , any))[apply(te.diff.results != 0  , 1 , any)])
addList(david, te_down_list, idType="ENSEMBL_GENE_ID", listName="te_down_list", listType="Gene")
addList(david, te_up_list, idType="ENSEMBL_GENE_ID", listName="te_up_list", listType="Gene")
addList(david, te_any, idType="ENSEMBL_GENE_ID", listName="te_any", listType="Gene")
AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)

# Way too many genes are up/down in rna/ribo so we might have to do something more specific
# One idea is to change to back to rna_any/ribo_any and ask if different degrees are enriched in something
rna_down_list = hgnc_to_ensg_convert(names(apply(results.rna == -1 , 1 , any))[apply(results.rna == -1 , 1 , any)])
rna_up_list = hgnc_to_ensg_convert(names(apply(results.rna == 1 , 1 , any))[apply(results.rna == 1 , 1 , any)])
rna_any  = hgnc_to_ensg_convert(names(apply(results.rna != 0 , 1 , any))[apply(results.rna != 0  , 1 , any)])
addList(david, rna_down_list, idType="ENSEMBL_GENE_ID", listName="rna_down_list", listType="Gene")
addList(david, rna_up_list, idType="ENSEMBL_GENE_ID", listName="rna_up_list", listType="Gene")
addList(david, rna_any, idType="ENSEMBL_GENE_ID", listName="rna_any", listType="Gene")
AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)

# Ribo-Down/ribo-any has slight enrichment for immune-related functions
ribo_down_list = hgnc_to_ensg_convert(names(apply(results.ribo == -1 , 1 , any))[apply(results.ribo == -1 , 1 , any)])
ribo_up_list = hgnc_to_ensg_convert(names(apply(results.ribo == 1 , 1 , any))[apply(results.ribo == 1 , 1 , any)])
ribo_any  = hgnc_to_ensg_convert(names(apply(results.ribo != 0 , 1 , any))[apply(results.ribo != 0  , 1 , any)])
addList(david, ribo_down_list, idType="ENSEMBL_GENE_ID", listName="ribo_down_list", listType="Gene")
addList(david, ribo_up_list, idType="ENSEMBL_GENE_ID", listName="ribo_up_list", listType="Gene")
addList(david, ribo_any, idType="ENSEMBL_GENE_ID", listName="ribo_any", listType="Gene")
AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)

## SOM enrichments -- We can do enrichment on each of the cells of the SOM or we can cluster the cells and enrichment on the cluster
# Absolute.som Enrichments
# absolute.som$unit.classif -> Where each gene goes
# cluster.membership -> Code clustering based on affinity propagation
# row.names(ribo_rna_te_prot)[!apply(is.na(abs.som.data), 1, any)] -> Gene_Ids
# Each cell enrichment
abs.som.genes = row.names(ribo_rna_te_prot)[!apply(is.na(abs.som.data), 1, any)]
background.list.abs.som = hgnc_to_ensg_convert(abs.som.genes )
addList(david, background.list.abs.som, idType="ENSEMBL_GENE_ID", listName="background.list.abs.som", listType="Background")

setAnnotationCategories (david, c("GOTERM_CC_ALL", "GOTERM_BP_ALL", "GOTERM_MF_ALL", "KEGG_PATHWAY", "REACTOME_PATHWAY"))

# # No results except for 4 units
# for (i in 1:(absolute.som$grid$xdim *absolute.som$grid$ydim) ) { 
# unit_list = hgnc_to_ensg_convert(abs.som.genes[absolute.som$unit.classif == i ])
# addList(david, unit_list, idType="ENSEMBL_GENE_ID", listName=paste("UnitList", i, sep="_" ), listType="Gene")
# setCurrentBackgroundPosition(david,3)
# AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
# FilteredChart = filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)
#   if (length(FilteredChart$Term) != 0L) { 
#     out.file = paste("Absolute.SOM.Unit", i, sep="_")
#     out.df = data.frame(Term=FilteredChart$Term, FE = FilteredChart$Fold.Enrichment, FDR= FilteredChart$FDR )  
#     write.table(out.df, file = paste('~/project/CORE_DATAFILES/GO_RESULTS/', out.file, sep=""),row.names=F)
#   }
# }



# Enrichment by clustered Affinity Propagation
# Cluster 4 -> Chromosome enrichhed high rna not too low te but low protein. Any reason proteomics?
# Cluster 8 => Super high RNA levels not very highly translated but reaches high protein levels
# Cluster 5 = > Much more highly translated RNAs in intersting categories.
for (i in 1:max(cluster.membership) ) { 
units_in_cluster = c(1:(absolute.som$grid$xdim *absolute.som$grid$ydim))[cluster.membership == i]
cluster_list = hgnc_to_ensg_convert(abs.som.genes[absolute.som$unit.classif %in% units_in_cluster])
out.file = paste("Absolute.SOM.Cluster", i, sep="_")
write.table(cluster_list, row.names=F, file = paste('~/project/CORE_DATAFILES/FUNCASSOCIATE/ABS_SOM/',out.file, sep="") )

addList(david, cluster_list, idType="ENSEMBL_GENE_ID", listName=paste("cluster_list", i, sep="_" ), listType="Gene")
setCurrentBackgroundPosition(david,2)
AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
FilteredChart = filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)
  if (length(FilteredChart$Term) != 0L) { 
    out.df = data.frame(Term=FilteredChart$Term, FE = FilteredChart$Fold.Enrichment, FDR= FilteredChart$FDR )  
    write.table(out.df, file = paste('~/project/CORE_DATAFILES/GO_RESULTS/', out.file, sep=""),row.names=F)    
  }
}

## Repeat these analysis with FUNCASSOCIATE
write.table(hgnc_to_ensg_convert (abs.som.genes), '~/project/CORE_DATAFILES/FUNCASSOCIATE/ABS_SOM/AbsSOMBackground', row.names=F)

  
### RELATIVE SOM ENRICHMENT
# Possibly three classes, high difference in spearmant correlation difference
# Those that have overall high correlation

# som.exp.prot$unit.classif
# min_cor_pval_unit_type
# min_cor_pval_unit
# max_cor_unit
# row.names(som.exp.prot$data$ribo)

## RELATIVE SOM IDS ARE EQUAL TO row.names(linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22])


relative.som.background = hgnc_to_ensg_convert(row.names(linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22]))
addList(david, relative.som.background, idType="ENSEMBL_GENE_ID", listName="relative.som.background", listType="Background")
# # Relative Som Enrichment for significant correlation difference
# # cor_diff_sig_units_rna_genes_index ; No significant enrichment for either RNA or Ribo
# cor_diff_sig_rna_ids = hgnc_to_ensg_convert (
#   row.names(linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22])[cor_diff_sig_units_rna_genes_index]
#   )
# addList(david, cor_diff_sig_rna_ids, idType="ENSEMBL_GENE_ID", listName="cor_diff_sig_rna_ids", listType="Gene")
# setCurrentBackgroundPosition(david,2)
# AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
# FilteredChart = filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)

# Relative SOM UnitWise Enrichment -- 3 cells with enrichment
for (i in 1:(supersom.nonzero.ind$grid$xdim *supersom.nonzero.ind$grid$ydim) ) { 
  unit_list = hgnc_to_ensg_convert(row.names(linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22])[supersom.nonzero.ind$unit.classif == i ] )
  addList(david, unit_list, idType="ENSEMBL_GENE_ID", listName=paste("RelativeSOMUnitList", i, sep="_" ), listType="Gene")
  setCurrentBackgroundPosition(david,2)
  AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
  FilteredChart = filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)
  if (length(FilteredChart$Term) != 0L) { 
    out.file = paste("Relative.SOM.Unit", i, sep="_")
    out.df = data.frame(Term=FilteredChart$Term, FE = FilteredChart$Fold.Enrichment, FDR= FilteredChart$FDR )  
    write.table(out.df, file = paste('~/project/CORE_DATAFILES/GO_RESULTS/', out.file, sep=""),row.names=F)
  }
}

## Relative SOM Enrichment grouped by best correlating feature
## Here we took a cumulative approach. We can also do a clustering or unit-wise approach
## There is a no significant enrichment with the cumulative approach
## No dependence on whether we threshold on the spearman

# pval_threshold = 0.05
# spearman_threhold = .5
# table(min_cor_pval_unit_type[which(p.adjust(min_cor_pval_unit, method = "holm") < pval_threshold)])
# rna_cor = c()
# ribo_cor = c()
# for ( i in which(p.adjust(min_cor_pval_unit, method = "holm") < pval_threshold & max_cor_unit > spearman_threhold)) { 
#  if (min_cor_pval_unit_type[i] == 1) { 
#    rna_cor = c(rna_cor, 
#                hgnc_to_ensg_convert(
#     row.names(linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22])[supersom.nonzero.ind$unit.classif == i ] ) ) 
#  }
#  else if (min_cor_pval_unit_type[i] == 2) { 
#    ribo_cor = c(ribo_cor, 
#                hgnc_to_ensg_convert(
#     row.names(linfeng_protein_ribo_rna[,type_prot=="Prot"][,-22])[supersom.nonzero.ind$unit.classif == i ] ) )             
#  }
# }
# length(ribo_cor)
# length(rna_cor)
# 
# addList(david, rna_cor, idType="ENSEMBL_GENE_ID", listName="rna_cor", listType="Gene")
# addList(david, ribo_cor, idType="ENSEMBL_GENE_ID", listName="ribo_cor", listType="Gene")
# setCurrentBackgroundPosition(david,2)
# AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.01, count=2L)
# FilteredChart = filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)
# FilteredChart$Term
# 
# significantly_correlated_rnacells <- cbind(som.exp.prot$codes$rna,som.exp.prot$codes$ribo, som.exp.prot$codes$prot)[which(p.adjust(min_cor_pval_unit) < pval_threshold & min_cor_pval_unit_type == 1),] 
# significantly_correlated_ribocells <- cbind(som.exp.prot$codes$rna,som.exp.prot$codes$ribo, som.exp.prot$codes$prot)[which(p.adjust(min_cor_pval_unit) < pval_threshold & min_cor_pval_unit_type == 2),] 
# 

##############END OF GO ANALYSIS

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
#pdf(file="~/Google_Drive/Manuscript Figures/Kozak_Analysis/Difference_in_Kozak_Scores.pdf", width=3, height=4, bg= "transparent")
boxplot(kozak_merge$V2 -kozak_merge$V4, notch=T)
abline(h = 0)
#dev.off()
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
legend("topright", paste("p-val = ",  signif(summary.lm(lm(ribo_only$E[my_index,] ~ index_factor, weights=ribo_only$weights[my_index,]))$coefficients[2,4], digits=2 ), sep="" ), inset=0.05, bty= "n" )
  boxplot(rna_only$E[my_index,]~rna_index_factor, ylab="RNA Expression", xlab= "Kozak Strength", names=unique(sort(round(index_factor, digits=2))))
legend("topright", paste("p-val = ",  signif(summary.lm(lm(rna_only$E[my_index,] ~ rna_index_factor, weights=rna_only$weights[my_index,]))$coefficients[2,4], digits=2 ), sep="" ), inset=0.05, bty= "n" )
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
    #pdf(file=paste("~/Google_Drive/Manuscript Figures/Kozak_Analysis/Ribo_", row.names(ribo_only)[my_index],".pdf", sep=""), width=5, height=5   )
    boxplot(ribo_only$E[my_index,]~ index_factor, ylab= "Ribosome Occupancy", xlab="Allele Number" , names=unique(sort(index_factor)), main=row.names(ribo_only)[my_index] )
    legend("bottomright", paste("p-val = ",  signif(my_pval, digits=2 ), sep="" ), inset=0.05, bty= "n" )
    #dev.off()
    rna_index_factor <- rep (0, times=length(sample_labels_rna))
    rna_index <-  grep(paste(ind_unique , collapse="|"), sample_labels_rna)
    rna_index_values <- grep(paste(ind_unique , collapse="|"), sample_labels_rna, value=T)
    rna_index_factor[rna_index] <-  allele_num[match(rna_index_values, ind_unique)] 
    rna_pval <- summary(lm(rna_only$E[my_index,]~ rna_index_factor, weights=rna_only$weights[my_index,]))$coefficients[2,4]
    #pdf(file=paste("~/Google_Drive/Manuscript Figures/Kozak_Analysis/RNA_", row.names(rna_only)[my_index], ".pdf",  sep=""), width=5, height=5   )
    boxplot(rna_only$E[my_index,]~ rna_index_factor, ylab= "RNA Expression", xlab="Allele Number" , names=unique(sort(index_factor)), main= row.names(rna_only)[my_index] )
    legend("topright", paste("p-val = ",  signif(rna_pval, digits=2 ), sep="" ), inset=0.05, bty= "n" )  
    #dev.off()
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
length(which(p.adjust(list_of_pval, method="hommel") < .05))
significant_ribo_diff <- which(p.adjust(list_of_pval, method="hommel") < .05)
color_by_pval <- rep(0, length(list_of_pval))
pval_cutoff <- max(list_of_pval[significant_ribo_diff])
color_by_pval[list_of_pval <= pval_cutoff] <- 1
plot(ribo_diff, kozak_diff, cex=0.65, pch=19, xlim=c(-1.2, 1.2), tck = .02, col=c("Black", "Red")[as.factor(color_by_pval)])
abline(v=c(0), h=c(0,log(2), -log(2)))

# Some figure to show that the significantly different Ribo Diff Ones have significant effect on Kozak
# We can just state this; when we take transcripts with ribo difference significant at 5% FDR
# Wilcox.test p-value is 0.02
boxplot(abs(kozak_diff[significant_ribo_diff]), abs(kozak_diff[!is.na(p.adjust(list_of_pval, method="hommel"))]))
wilcox.test(abs(kozak_diff[significant_ribo_diff]), abs(kozak_diff[!is.na(p.adjust(list_of_pval, method="hommel"))]))

#Experimental Results -- RATIOS
t1 <- c(3.96542345327789,
        3.44234121174092,
        3.04438144930569,
        2.93967231943929,
        4.11347731000547)
t2 = c(2.82217506138878,
       2.59405879757123,
       2.35815546097475,
       2.963013925323,
       2.28901207104044) 
t3= c(3.60185503263483,
      3.62021809898122,
      3.28314793141455,
      2.99207787282458,
      3.50478602453603)
t4 = c(18.2902479279538,
       20.6546019318151,
       18.8193020249892,
       16.6166007581944,
       21.3674390573499)
t5= c(11.0387066961044,
      11.8931431823988,
      11.6826039464806,
      10.1285659460151,
      12.0310466442286)
t6= c(17.0039823022692,
      17.6715785418472,
      17.4222866730181,
      17.3597159803328,
      21.5216655304032)
t7 = c(5.54029043612037,
       3.79875677177386,
       4.34192690011633,
       5.46527638585086,
       5.13069328479285)
t8= c(2.71925726758724,
      2.01505284250671,
      2.26591664908055,
      2.3069549139525,
      3.65752843380776)
t9 = c(3.17805617165845,
       3.05691918246317,
       4.87943869472912,
       2.98700902762487,
       6.17002355323692)
t10 = c(12.7397566904909,
        13.0141868294551,
        15.0902961231795,
        15.2084415898006,
        20.5442642085684)
t11 = c(11.6546081495051,
        10.5407431725019,
        11.5617412714154,
        9.09057275966018,
        13.2427496147232)
# Based on http://nar.oxfordjournals.org/content/32/20/e160.full#disp-formula-2
# we define outliers as 
outlier_detect<- function(r) {
  range = c(median(r) + 1.5*IQR(r), median(r) - 1.5*IQR(r))
  outliers= !(r > range[1] | r < range[2])
  return(outliers)
}
# 4-6-10
#1-7 ; 2-8, 3-9, 5-11
# Welch two sample t-test
t.test(t1, t7)
t.test(t1[outlier_detect(t1)], t7[outlier_detect(t7)])
# t6, and t10 have clear outliers that can be removed
t.test(t4, t10)
t.test(t4[outlier_detect(t4)], t10[outlier_detect(t10)])
t.test(t6, t10)
t.test(t6[outlier_detect(t6)], t10[outlier_detect(t10)])
t.test(t4, t6)
t.test(t4[outlier_detect(t4)], t6[outlier_detect(t6)])

df = matrix( c(t1, t2 , t3 , t4, t5, t6, t7, t8, t9, t10, t11), ncol=11)
outs <- apply(df, 2, outlier_detect)
df.summary  = data.frame(MEAN =apply(df, 2, mean), 
                         SE = apply(df, 2, function(x){sd(x)/ sqrt(length(x))}), 
                         PAIRING = factor (c(1, 2, 3, 4,5,4,1,2,3,4,5 ) )  ) 

barplot2(df.summary$MEAN[c(1,7)], space = 0, col=c("blue", "salmon"), ylim=c(0,5.5), ylab="Rluc/Fluc Ratio",
         plot.ci=T, ci.l= df.summary$MEAN[c(1,7)] - df.summary$SE[c(1,7)], ci.u = df.summary$MEAN[c(1,7)] + df.summary$SE[c(1,7)])

barplot2(df.summary$MEAN[c(6,4,10)], space = 0, col=c("blue", "salmon", "green"), ylim=c(0,20), ylab="Rluc/Fluc Ratio",
         plot.ci=T, ci.l= df.summary$MEAN[c(6,4,10)] - df.summary$SE[c(6,4,10)], ci.u = df.summary$MEAN[c(6,4,10)] + df.summary$SE[c(6,4,10)])


#ggplot(df.summary, aes( PAIRING, MEAN)) + 
#  geom_bar( stat="identity")

# Take the translation efficiency table and for each position look for significant diff with Kruskal
# Generate this table by merging kozak data with translation efficiency table (grand_mean_te)
kozak_seq_score_table <- read.table('~/project/CORE_DATAFILES/Kozak_IDs_PWM_Strand_Seq_HGNC_TE_Ribo.bed')

kruskal_pvals = c()
# Translation Efficiency and Ribosome Occupancy are different somewhat
f <- function(s, letter) strsplit(s, "")[[1]][letter]
g <- function(s) strsplit(s, "")[[1]][c(4,10)]
for ( j in c(1:6, 10:11)) { 
seq_factor <- sapply(as.character(kozak_seq_score_table$V7), f, letter=j)
k1 <- kruskal.test(kozak_seq_score_table$V9 ~ as.factor(seq_factor))
kruskal_pvals = c(kruskal_pvals, k1$p.value)
print (k1$p.value * 8)
print(kruskalmc(kozak_seq_score_table$V9 ~ as.factor(seq_factor), probs=.05/8))
if (j < 7) { 
  pos = j -7
}
else { 
  pos = j -6
}
#pdf (paste('~/Google Drive/Manuscript Figures/Kozak_Analysis/Translation_Efficiency_by_Position',pos ,sep="_"), width=5, height=5 )
boxplot(kozak_seq_score_table$V9 ~ as.factor(seq_factor), 
        varwidth=T, ylim=c(-.75,.75), notch=T, cex=.2, whisklty=0, staplelty=0, ylab="Translation Efficiency", main = paste("Position", pos , sep = ": ") )

# boxplot(kozak_seq_score_table$V9 ~ as.factor(seq_factor), 
#         varwidth=T, cex=.2, outline = F, ylab="Translation Efficiency", main = paste("Position", pos , sep = ": ") )
# # max(boxplot.stats(kozak_seq_score_table$V9)$out[< 0])
# # [1] -2.238098
# # > min(boxplot.stats(kozak_seq_score_table$V9)$out[a1])
# # [1] 2.367177
# collapse_out = kozak_seq_score_table$V9
# collapse_out[collapse_out > 2.367177 ] = 2.367177
# collapse_out[collapse_out < -2.238098 ] = -2.238098
# bwplot(collapse_out ~ as.factor(seq_factor), panel=function(...){panel.violin(...); panel.bwplot(do.out = F, ...)})
# #dev.off()

#boxplot(kozak_seq_score_table$V10 ~ as.factor(seq_factor), varwidth=T, ylim=c(4,6), notch=T, range=.001, cex=.2)
}

# Add Heatmap of p-values from the kruskal vallis
# Use rect for drawing; Use -log10p 
adjusted_p_kruskal = log10(p.adjust (kruskal_pvals, method = "bonferroni"))
# We can also draw 6*7 matrix with p-values for each change
# Fill color = colorpanel(1024,'blue', 'red')
i = rep(1,1+abs(min(adjusted_p_kruskal)))
j = seq(1,1+abs(min(adjusted_p_kruskal)))
# 20 color heatmap legend
colors = c( colorpanel(17,'red', "#FFEFEF"), colorpanel(3,'white', 'blue') )
#pdf ("~/Google_Drive/Manuscript Figures/Kozak_Analysis/Color_Key_Kozak_Significance.pdf", width=3, height=3)
plot(NA, type = "n", ann=FALSE, xlim=c(1,2),ylim=c(1,2+abs(min(adjusted_p_kruskal))),xaxt="n",yaxt="n",bty="n")
rect(i,0+j, i+1, 1+j, col =colors)
#dev.off()
# Log-pvalue of the effect on translation efficiency summarized by position
p_corresponding_colors = rev(colors)[round(abs(adjusted_p_kruskal))+1]
i_cor = rep(1, length(p_corresponding_colors))
j_cor = seq(1, length(p_corresponding_colors))

#pdf ("~/Google_Drive/Manuscript Figures/Kozak_Analysis/KozakPosition_Significance.pdf", width=3, height=3)
plot(NA, type = "n", ann=FALSE, xlim=c(1,2),ylim=c(1,2+abs(min(adjusted_p_kruskal))),xaxt="n",yaxt="n",bty="n")
rect(i_cor,0+j_cor, i_cor+1, 1+j_cor, col = p_corresponding_colors)
#dev.off()

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
# 
# # GO over each gene. Calculate F-value for ribo and rna separately
# # One issue is that mean diff is highly correlated with p-value.
# # The higher the mean diff, the higher the p-value
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
# ####################################
# #END OF OLD VARIATION ANALYSIS


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

### UNUSED SOM PREDICTION 
# We can do a version of this where there are more classes based on differences
# Test prediction -- Not doing much better than rna alone
# set.seed(173)
# training <- sample(c(T,F), nrow(abs.som.data.noNA),  prob=c(.95,.05), replace=T)
# absolute.som.train <- som(abs.som.data.noNA[training,1:3], toroidal=T, grid=somgrid(xdim, ydim, "hexagonal"))
# Xtest= abs.som.data.noNA[!training,1:3]
# Ytrain=abs.som.data.noNA[training,4]
# YTrue= abs.som.data.noNA[!training,4]
# absolute.predictions <- predict.kohonen(absolute.som.train, newdata=Xtest, trainY=Ytrain)
# diff_true_pred = abs(absolute.predictions$prediction - YTrue )
# diff_rna_true = abs( abs.som.data.noNA[!training,2] - YTrue)
# diff_random_true = abs( runif(length(YTrue)) - YTrue)
# hist(diff_random_true, 25, xlim = c(0,1))
# hist(diff_rna_true, 25, xlim = c(0,1))
# hist(diff_true_pred, 25, xlim = c(0,1))
# length(which(diff_rna_true < .2))
# length(which(diff_true_pred < .2))

#plot(absolute.som)
#plot(absolute.som, type="quality")
#plot(absolute.som, type="count")
#plot(absolute.som, type="changes")


# PREDICTION DOESN'T WORK THAT MUCH BETTER THAN RNA ALONE HERE. 
# The training set and test set don't have the same distribution of prot quantiles
# Xtest <- list ( ribo= ribo.quantiles[singleNARows, -dropCols] , 
#                 rna = rna.quantiles[singleNARows, -dropCols] ,
#                 te = te.quantiles[singleNARows, -dropCols],
#                 prot=prot.quantiles[singleNARows, -dropCols])
# 
# prot.predictions <- predict.kohonen(som.exp.prot, newdata=Xtest, whatmap=c(1,2,3) )
# diff_in_pred <- prot.predictions$prediction - prot.quantiles[singleNARows,-dropCols]
# 
# THRESHOLD <- .1
# sum( abs(apply(diff_in_pred, 1, median, na.rm=T)) < THRESHOLD) 
# random.prediction <- matrix(runif(3708), ncol= ncol(prot.predictions$prediction))
# diff_wrandom <- random.prediction - prot.quantiles[singleNARows,-dropCols]
# sum( abs(apply(diff_wrandom, 1, median, na.rm=T)) < THRESHOLD) 
# diff_wrna <- rna.quantiles[singleNARows, -dropCols] - prot.quantiles[singleNARows,-dropCols]
# diff_betweenpred_and_random <-abs(diff_in_pred) - abs(diff_wrandom)
# diff_betweenpred_and_rna <- abs(diff_in_pred) - abs(diff_wrna)



############UNUSED ACROSS INDIVIDUAL VARIATION ANALYSIS
# # We should test the effect of unbalanced design: 
# # Idea is subsample libraries so that there are 2 of each 
# # Test robustness to sampling
# # The p-values are highly unstable based on the sampling
# # table(as.factor(sample_id_all[1:49]))
# # table(as.factor(sample_id_all[50:80]))
# rna_cols_to_select = c()
# ribo_cols_to_select = c()
# for (l in levels(as.factor(sample_id_all[1:49])) ) { 
#   rna_cols_to_select= c( rna_cols_to_select, 
#                          sample (which(as.factor(sample_id_all[1:49]) == l) , 2, replace=F) )
# }
# for (l in levels(as.factor(sample_id_all[50:80])) ) { 
#   ribo_cols_to_select= c( ribo_cols_to_select, 
#                          49 + sample (which(as.factor(sample_id_all[50:80]) == l) , 2, replace=F) )
# }
# table(sample_id_all[rna_cols_to_select])
# table(sample_id_all[ribo_cols_to_select])
# 
# # We can use eta2 as measure of effect size SS-W/SS-TOTAL
# ribo_F = c()
# ribo_eta2 = c()
# rna_F = c()
# rna_eta2 = c()
# rna_cols_to_select = 1:49
# ribo_cols_to_select = 50:80
# for ( i in 1:nrow(joint_expression_common$E)) { 
#   # Summary keeps sum of squares as $ Sum Sq : num  2.82 3.4
#   tmp_ribo = summary(aov(joint_expression_common$E[i,ribo_cols_to_select] ~ as.factor(sample_id_all[ribo_cols_to_select]), 
#                          weights= joint_expression_common$weights[i,ribo_cols_to_select] ))[[1]]
#   ribo_F = c(ribo_F, tmp_ribo$Pr[1] ) 
#   ribo_eta2 = c(ribo_eta2, tmp_ribo[1,2]/ sum(tmp_ribo[,2]) )
#   tmp_rna =  summary(aov(joint_expression_common$E[i,rna_cols_to_select] ~ as.factor(sample_id_all[rna_cols_to_select]), 
#                          weights= joint_expression_common$weights[i,rna_cols_to_select] ))[[1]]   
#   rna_F = c(rna_F, tmp_rna$Pr[1] ) 
#   rna_eta2 = c(rna_eta2,tmp_rna[1,2]/ sum(tmp_rna[,2]) )                 
# }
# 
# ribo_F_corrected = p.adjust(ribo_F, method="holm")
# rna_F_corrected = p.adjust(rna_F, method = "holm")
# ribo_sig = which (ribo_F_corrected < 0.05)
# rna_sig = which (rna_F_corrected < 0.05)
# joint_sig = which (rna_F_corrected < 0.05 & ribo_F_corrected < 0.05)
# length(ribo_sig)
# length(rna_sig)
# length(joint_sig)
# plot(rna_eta2, rna_F_corrected, cex = .2, pch = 19)
# plot(ribo_eta2, ribo_F_corrected, cex = .2, pch = 19)

### THIS APPROACH BASED ON CVs or SIMPLE RATIOS IS COMMENTED OUT
# # Improved method using weighted coefficient of variation
# 
# # weighted_coef_var is a function to be applied to limma object with E and weights
# weighted_coef_var <- function(x)  { 
# #  wt_cv = c()
# #  pb <- tkProgressBar(title="Progress Bar", min = 0, max = dim(x$E)[1], width=300)
#   wt_cv_unbiased = c()
#   for (i in 1:dim(x$E)[1]) { 
#     wt_mean = weighted.mean(x$E[i,] , x$weights[i,] )
#     # Weighted sample variance --biased estimator
# #    wt_var = sum(x$weights[i,] * (x$E[i,] - wt_mean)^2) / sum(x$weights[i,])
#     # Weighted sample variance -- unbiased estimator
#     wt_var_unbiased = (sum(x$weights[i,]) * (sum(x$weights[i,] * (x$E[i,] - wt_mean)^2) ) ) / 
#       (sum(x$weights[i,])^2 - sum(x$weights[i,]^2) ) 
# #   wt_cv = c(wt_cv, sqrt(wt_var) / wt_mean)
#     wt_cv_unbiased = c(wt_cv_unbiased, sqrt(wt_var_unbiased) / wt_mean)    
# #  setTkProgressBar(pb, i, label=paste( round(i/dim(x$E)[1]*100, 0),"% done"))
#   }
# #  return (wt_cv)
# #  close (pb)
#   return (wt_cv_unbiased)
# }
# 
# # Given limma object calculates weighted mean and weighted weights
# 
# weighted_mean_limma <- function(x)  { 
#   wt_mean = c()
#   for (i in 1:dim(x$E)[1]) { 
#     wt_mean = c(wt_mean, weighted.mean(x$E[i,] , x$weights[i,] ))
#   }
#   return (wt_mean)
# }
# 
# weighted_weight_limma <- function(x)  { 
#   wt_weights =c()
#   for (i in 1:dim(x$E)[1]) { 
#     wt_weights = c(wt_weights, weighted.mean(x$weights[i,], x$weights[i,]))
#   }
#   return (wt_weights)
# }
# 
# #coef_var_median
# wt_cvs_ribo = data.frame(matrix(0, ncol = length(unique(sample_labels_joint_common[type_common=="Ribo"])),
#                      nrow = nrow(joint_expression_common[,type_common=="Ribo"]) ) )
# colnames(wt_cvs_ribo) <- unique(sample_labels_joint_common[type_common=="Ribo"])
# 
# wt_means_ribo = wt_cvs_ribo
# wt_weights_ribo = wt_cvs_ribo
# 
# for (j in unique(sample_labels_joint_common[type_common=="Ribo"]) ){ 
#   cols <- which(sample_labels_joint_common[type_common=="Ribo"] == j)
#   wt_cvs_ribo[[j]] = weighted_coef_var(joint_expression_common[,type_common=="Ribo"][,cols])
#   wt_means_ribo[[j]] = weighted_mean_limma(joint_expression_common[,type_common=="Ribo"][,cols])
#   wt_weights_ribo[[j]] = weighted_weight_limma(joint_expression_common[,type_common=="Ribo"][,cols])
# }
# 
# wt_cvs_rna = data.frame(matrix(0, ncol = length(unique(sample_labels_joint_common[type_common=="RNA"])),
#                                 nrow = nrow(joint_expression_common[,type_common=="RNA"]) ) )
# colnames(wt_cvs_rna) <- unique(sample_labels_joint_common[type_common=="RNA"])
# wt_means_rna = wt_cvs_rna
# wt_weights_rna = wt_cvs_rna
# 
# for (j in unique(sample_labels_joint_common[type_common=="RNA"]) ){ 
#   cols <- which(sample_labels_joint_common[type_common=="RNA"] == j)
#   wt_cvs_rna[[j]] = weighted_coef_var(joint_expression_common[,type_common=="RNA"][,cols])
#   wt_means_rna[[j]] = weighted_mean_limma(joint_expression_common[,type_common=="RNA"][,cols])
#   wt_weights_rna[[j]] = weighted_weight_limma(joint_expression_common[,type_common=="RNA"][,cols])
#   
# }
# 
# median_ribo_cv <- apply(wt_cvs_ribo, 1, median)
# median_rna_cv <- apply(wt_cvs_rna, 1, median)
# hist(median_ribo_cv, 50)
# hist(median_rna_cv, 50)
# 
# # NEED TO MAKE SURE THAT THESE REMAINING SECTIONS WORK: 
# per_ind_ribo_cv = weighted_coef_var(list(E =  wt_means_ribo, weights = wt_weights_ribo))
# per_ind_rna_cv = weighted_coef_var(list(E = wt_means_rna, weights= wt_weights_rna))
# hist(per_ind_ribo_cv[per_ind_ribo_cv < 1], 50)
# hist(per_ind_rna_cv[per_ind_rna_cv>0 & per_ind_rna_cv < 1], 50)
# 
# inter_ind_cv_to_technical_ribo = per_ind_ribo_cv / median_ribo_cv
# inter_ind_cv_to_technical_rna = per_ind_rna_cv / median_rna_cv
# hist(inter_ind_cv_to_technical_ribo[inter_ind_cv_to_technical_ribo < 10], 50)
# hist(inter_ind_cv_to_technical_rna[inter_ind_cv_to_technical_rna > 0 &inter_ind_cv_to_technical_rna < 10], 50)
# 
# ratio_of_inter_ind_cv <- inter_ind_cv_to_technical_rna/inter_ind_cv_to_technical_ribo
# # which(ratio_of_inter_ind_cv < 0 )
# # [1] 3251 3297 3790 5269 8037
# # which(ratio_of_inter_ind_cv> 5)
# # [1]  932 1067 1442 1571 2752 3347 3523 4289 5005 5551 6233 7818 8616 9185
# # which(ratio_of_inter_ind_cv> 10)
# # [1]  932 1067 5551 8616
# 
# hist(ratio_of_inter_ind_cv[ratio_of_inter_ind_cv > 0 & ratio_of_inter_ind_cv < 5], 50 )
# 
# # Median_[ribo/rna]_cv is our estimated technical noise
# # For interindividual cv: Take weighed mean of expression and weights for each individual
# # Then use these to calculated weighted cv
# 
#   
# # There doesn't seem to be any significant association between expression and the cv ratios as expected
# # We can add a permutation scheme here to get p-values on these differences
# # The easiest permutation scheme is to permute type_common to generate new RNA, Ribo Classification
# # Run through existing code and compare the actual difference in CV witht the permutation p-value
# # We might rely on nominal p-value threshold for selecting significant ones
# # type_common will be modified making sure that each individual is preserved as in sample_labels_joint_common
# # Loop over unique(sample_labels_joint_common)
# weighted_sd <- function (y) { 
#   return ( sd(y[,2]*y[,4]) )
# }
# 
# rna_replicate_mean <- apply (joint_expression_common$E[,type_common=="RNA"], 1, function(x) {
# aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="RNA"])), mean)  
# } )
# ribo_replicate_mean <- apply (joint_expression_common$E[,type_common=="Ribo"], 1, function(x) {
#   aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="Ribo"])), mean)  
# } )
# rna_replicate_weight_mean <- apply (joint_expression_common$weights[,type_common=="RNA"], 1, function(x) {
#   aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="RNA"])), mean)  
# } )
# ribo_replicate_weight_mean <- apply (joint_expression_common$weights[,type_common=="Ribo"], 1, function(x) {
#   aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="Ribo"])), mean)  
# } )
# rna_replicate_mean_weights <- mapply(cbind, rna_replicate_mean, rna_replicate_weight_mean, SIMPLIFY=F)
# ribo_replicate_mean_weights <- mapply(cbind, ribo_replicate_mean, rna_replicate_weight_mean, SIMPLIFY=F)
# 
# rna_cv_between_individuals <- as.numeric(lapply(rna_replicate_mean_weights, weighted_sd))
# ribo_cv_between_individuals <- as.numeric(lapply(ribo_replicate_mean_weights, weighted_sd))
# 
# # Calculate median CV of replicate CVs
# # THIS SHOULD BE SWITCHED WITH WEIGHTED.SD/ WEIGHTED.MEAN AND THEN TAKE THE MEDIAN
# rna_replicatecvs <- apply (joint_expression_common$E[,type_common=="RNA"]*joint_expression_common$weights[,type_common=="RNA"], 1, function(x) {
#   aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="RNA"])), sd)  
# } )
# ribo_replicatecvs <- apply (joint_expression_common$E[,type_common=="Ribo"]*joint_expression_common$weights[,type_common=="Ribo"], 1, function(x) {
#   aggregate(x, by= list(as.factor(sample_labels_joint_common[type_common=="Ribo"])), sd)  
# } )
# 
# rna_repcv_median <- as.numeric(lapply(rna_replicatecvs, function(z){median(z$x)}))
# ribo_repcv_median <- as.numeric(lapply(ribo_replicatecvs, function(z){median(z$x)}))
# 
# rnacv <- (rna_cv_between_individuals/rna_repcv_median) 
# ribocv <- (ribo_cv_between_individuals/ribo_repcv_median)
# 
# #pdf(file = "~/Google_Drive/Manuscript Figures/RNA_Between_Individual_Variance.pdf", width=5, height=5)
# #99% of the dataset is less than 8; so limit xlim to 0_to_8
# # We changed number of breaks so that the number of breaks in x-axis is similar
# hist(rnacv, 100, main= "RNA Expression", xlab = "Between/ Within Individual Coefficient of Variation", xlim=c(0,8))
# #dev.off()
# #pdf(file = "~/Google_Drive/Manuscript Figures/Ribo_Between_Individual_Variance.pdf", width=5, height=5)
# hist(ribocv, 200 , main= "Ribosome Occupancy", xlab = "Between/ Within Individual Coefficient of Variation", xlim=c(0,8))
# #dev.off()
# 
# #low_sig_to_noise <- which( rnacv < 1 & ribocv < 1)
# #names(rna_replicate_mean_weights)
# length(which(rnacv/ribocv > 2))
# length(which(rnacv/ribocv < .5))
# sorted_cvs <- sort(rnacv/ribocv , index.return=T)
# #write.table(row.names(joint_expression_common)[sorted_cvs$ix], file = "~/project/CORE_DATAFILES/Sorted_InterIndividualCV_RNA_to_Ribo.txt", row.names=F)
# 
# p1 <- hist(rnacv/ribocv, 100 )
# plot(p1, col=rgb(0,0,1,1/4),  tck=.02, xlab="RNA Expression to Ribosome Occupancy Ratio", main="Between/ Within Individual CV", ylim=c(0,500))
# abline(v=1, lwd=3)

# plot(p.adjust(random_effect_p_val_ribo, method = "holm") , ribo_F_corrected, pch = 19, cex=.2)
# plot(p.adjust(random_effect_p_val_rna, method = "holm") , rna_F_corrected, pch = 19, cex=.2)


# # To test whether sources of data increases across individual variance
# # Simple test is to select only polyA_RNA individuals and equivalents in Ribo
# # The results are consistent
# ribo_polyA = c()
# rna_polyA = c()
# for ( i in 1:nrow(joint_expression_common$E)) { 
#   ribo_polyA = c(ribo_polyA, 
#              summary(aov(joint_expression_common$E[i,c(50,51,54:57, 74:80)] ~ as.factor(sample_id_all[c(50,51,54:57, 74:80)]), 
#                          weights= joint_expression_common$weights[i,c(50,51,54:57, 74:80)] ))[[1]]$Pr[1] )
#   
#   rna_polyA = c(rna_polyA, 
#             summary(aov(joint_expression_common$E[i,1:18] ~ as.factor(sample_id_all[1:18]), 
#                         weights= joint_expression_common$weights[i,1:18] ))[[1]]$Pr[1] )
#   
# }

### OLD ACROSS INDIVIDUAL SUPERSOM ANALYSIS
linfeng_prot_common_with_te <-  colnames(linfeng_protein) %in% colnames(te_fit3$coefficients)
linfeng_te_columns <- linfeng_protein[,linfeng_prot_common_with_te]
class(linfeng_te_columns) <- "numeric"
linfeng_te_match <- cbind ( linfeng_te_columns[,1], rep(NA, dim(linfeng_te_columns)[1]), linfeng_te_columns[,2:5], rep(NA, dim(linfeng_te_columns)[1]), linfeng_te_columns[,6:12])
colnames(linfeng_te_match) <- sort(colnames(te_fit3$coefficients))
linfeng_te_match <- merge(linfeng_te_match, ensg_hgnc, by.x="row.names", by.y="ENSG")
# joint_count_ids is the same as row.names(v3); Pad the data with NAs when there is no proteomics
linfeng_te_match <- merge (data.frame(HGNC=joint_count_ids), linfeng_te_match, by="HGNC", all.x=T)
row.names(linfeng_te_match) <- linfeng_te_match[,1]
linfeng_te_match <- linfeng_te_match[,-c(1,2)]

# Simplifies the colnames; no order problems
colnames(ribo_fit2$coefficients) <- sort(colnames(te_fit3$coefficients))
colnames(rna_fit2$coefficients) <- sort(colnames(te_fit3$coefficients))
te_matrix = te_fit3$coefficients[,sort(colnames(te_fit3$coefficients), index.return=T)$ix]
# The easiest way for interpretation is using quantized data
ribo.quantiles = matrix(ecdf(ribo_fit2$coefficients)(ribo_fit2$coefficients), ncol = 14, dimnames = list(rownames(ribo_fit2$coefficients), colnames(ribo_fit2$coefficients)))
rna.quantiles = matrix(ecdf(rna_fit2$coefficients)(rna_fit2$coefficients), ncol = 14, dimnames = list(rownames(rna_fit2$coefficients), colnames(rna_fit2$coefficients)))
te.quantiles = matrix(ecdf(te_matrix)(te_matrix), ncol = 14, dimnames = list(rownames(te_matrix), colnames(te_matrix)))
prot.quantiles = matrix(ecdf(as.matrix(linfeng_te_match))(as.matrix(linfeng_te_match)), ncol = 14, dimnames = list(rownames(linfeng_te_match), colnames(linfeng_te_match)))
dropCols <- c(2,7)
dropRows <- !apply(is.na(prot.quantiles[,-dropCols]), 1, any)
#singleNARows <- which(apply(is.na(prot.quantiles[,-dropCols]),1, sum) == 1)

som.data = list (  ribo= ribo.quantiles[dropRows, -dropCols] , rna = rna.quantiles[dropRows, -dropCols] ,te = te.quantiles[dropRows, -dropCols])
som.data.prot = list ( ribo= ribo.quantiles , rna = rna.quantiles ,te = te.quantiles, prot=prot.quantiles)
som.data.prot.noNA = list ( Ribosome_Occupancy= ribo.quantiles[dropRows, -dropCols] , 
                            RNA_Expression = rna.quantiles[dropRows, -dropCols] ,
                            Translation_Efficiency = te.quantiles[dropRows, -dropCols], 
                            Protein_Expression=prot.quantiles[dropRows, -dropCols])

# Another version of the SuperSOM is using linfeng_protein_ribo_rna
# Sample IDs => sample_labels_joint_prot
# Types_of_Data => Type_Prot
# rna_replicate_mean_prot; ribo_replicate_mean_prot: Lists 
# rna_in_prot <- as.character(rna_replicate_mean_prot[[1]]$Group.1) %in% sample_labels_joint_prot[type_prot=="Prot"]
# ribo_in_prot <- as.character(ribo_replicate_mean_prot[[1]]$Group.1) %in% sample_labels_joint_prot[type_prot=="Prot"]

# Create num mat with each individual specific data
rna.num.mat = matrix(nrow=length(rna_replicate_mean_prot) , ncol = 27)
ribo.num.mat = matrix ( nrow = length(ribo_replicate_mean_prot), ncol = 28)
for (i in 1:length(rna_replicate_mean_prot)) { 
  rna.num.mat[i, ] = rna_replicate_mean_prot[[i]]$x[rna_in_prot]
}
for (i in 1:length(rna_replicate_mean_prot)) { 
  ribo.num.mat[i, ] = ribo_replicate_mean_prot[[i]]$x[ribo_in_prot]
}
colnames( rna.num.mat) = as.character (rna_replicate_mean_prot[[i]]$Group.1[rna_in_prot] ) 
colnames( ribo.num.mat) = as.character (ribo_replicate_mean_prot[[i]]$Group.1[ribo_in_prot] ) 

# Ribo and Prot has one extra col
dropProt = c(23)
ribo.num.mat  = ribo.num.mat[,-dropProt]
prot.num.mat = matrix(as.numeric(as.matrix(linfeng_protein_ribo_rna[,type_prot=="Prot"])),ncol=28)
prot.num.mat = prot.num.mat[, -dropProt]
rna.full.quantiles = matrix(ecdf(rna.num.mat)(rna.num.mat), ncol = 27 )
ribo.full.quantiles = matrix(ecdf(ribo.num.mat)(ribo.num.mat), ncol = 27 )
prot.full.quantile = matrix( ecdf (prot.num.mat)(prot.num.mat), ncol = 27)

supersom.fullrna.ribo.prot = list ( 
  ribo= ribo.full.quantiles , rna = rna.full.quantiles , prot=prot.full.quantile)



#total_cells <- floor(sqrt(length(som.data.prot)/2) * sqrt (dim(som.data.prot$ribo)[1] * dim(som.data.prot$ribo)[2]))
total_cells.noNA <- floor(sqrt(length(som.data.prot.noNA)/2) * 
                            sqrt (dim(som.data.prot.noNA$Ribosome_Occupancy)[1] * 
                                    dim(som.data.prot.noNA$Ribosome_Occupancy)[2]))
# if (floor(sqrt(total_cells/1.3333)) %% 2 == 0) { 
#   ydim.total = floor(sqrt(total_cells/1.3333))
# } else { 
#   ydim.total = floor(sqrt(total_cells/1.3333)) + 1
# }
# xdim.total = floor(total_cells/ydim.total + 0.5)

if (floor(sqrt(total_cells.noNA/1.3333)) %% 2 == 0) { 
  ydim.total.noNA = floor(sqrt(total_cells.noNA/1.3333))
} else { 
  ydim.total.noNA = floor(sqrt(total_cells.noNA/1.3333)) + 1
}
xdim.total.noNA = floor(total_cells.noNA/ydim.total.noNA + 0.5)

#som.exp = supersom(data =som.data, grid=somgrid(xdim.total.noNA, ydim.total.noNA, "hexagonal"), toroidal=T, contin=T)
# plot(som.exp, type="codes")
# plot(som.exp, type="quality")
# plot(som.exp, type="mapping", pch=19, cex=.3)
# plot(som.exp, type="changes")
# plot(som.exp, type="counts")
# plot.kohonen(som.exp, type="property", property=ribo_prot_cor_in_som, main= "Between Individual Ribosome Occupancy Protein Level Correlation", palette.name=redblue_cols, contin=T,zlim=c(-1,1), ncolors=11)
# plot.kohonen(som.exp, type="property", property=rna_prot_cor_in_som, main = "Between Individual RNA Occupancy Protein Level Correlation", palette.name=redblue_cols,contin=T, zlim=c(-1,1),ncolors=11)
# plot.kohonen(som.exp, type="property", property=te_prot_cor_in_som, main = "Between Individual Translation Efficiency Protein Level Correlation", palette.name=redblue_cols,contin=T, zlim=c(-1,1),ncolors=11)
# plot.kohonen(som.exp, type="counts" )

# We can update the xdim - ydim
total_cells.full <- floor(sqrt(length(supersom.fullrna.ribo.prot)/2) * 
                            sqrt (dim(supersom.fullrna.ribo.prot$prot)[1] * 
                                    dim(supersom.fullrna.ribo.prot$prot)[2]))

if (floor(sqrt(total_cells.full/1.3333)) %% 2 == 0) { 
  ydim.total.full = floor(sqrt(total_cells.full/1.3333))
} else { 
  ydim.total.full = floor(sqrt(total_cells.full/1.3333)) + 1
}
xdim.total.full = floor(total_cells.full/ydim.total.full + 0.5)

supersom.full_sample = supersom (data = supersom.fullrna.ribo.prot, 
                                 grid=somgrid ( xdim.total.full, ydim.total.full, "hexagonal"), toroidal=T, contin = T)
plot.kohonen (supersom.full_sample, type = "changes")

som.exp.prot.noRibo = supersom(data =som.data.prot.noNA, whatmap = 2:4, grid=somgrid(xdim.total.noNA, ydim.total.noNA, "hexagonal"), toroidal=T, contin=T)
# Give proteins higher weight for tighter clustering -> Another rationale is colinearity between the expression measures
# Given RNA and Ribo => TE is fixed. So, we should all gie the independent components 1/3 weight
# mean_distance <- 1
# iter = 500
# # pb <- tkProgressBar(title="Progress Bar", min = 0, max = iter, width=300)
# for (i in 11:iter) { 
#   set.seed(i)
#   relative.som <- supersom(data =som.data.prot.noNA, grid=somgrid(xdim.total.noNA, ydim.total.noNA, "hexagonal"), toroidal=T, contin=T, weights=c(2/9,2/9,2/9, 1/3))
#   if (mean(relative.som$distances) < mean_distance) { 
#     my_seed <- i
#     mean_distance <- mean(relative.som$distances)
#   }
# #   setTkProgressBar(pb, i, label=paste( round(i/iter*100, 0),"% done"))
#  }
# close (pb)
# Best in 500 seeds is 378
set.seed(378)
som.exp.prot = supersom(data = som.data.prot.noNA, grid=somgrid(xdim.total.noNA, ydim.total.noNA, "hexagonal"), toroidal=T, contin=T, weights=c(2/9,2/9,2/9, 1/3))

plot.kohonen(som.exp.prot, type="counts")
plot.kohonen(som.exp.prot, type="changes")
plot.kohonen(som.exp.prot, type="quality")
mean(som.exp.prot$distances)
pdf(file= "~/Google_Drive/Manuscript Figures/Across_Individual_Comparison/SuperSomCodes.pdf", width=8.5, height=11)
par(mfrow = c(2, 2))
plot.kohonen(som.exp.prot, type="codes")
dev.off()

# plot.kohonen(som.exp.prot, property=som.exp.prot$codes$ribo, type = "property", palette.name=redblue_cols, ncolors=11, contin=T)


## UNUSED GO ANALYSIS
# RNA/RIBO
high_rna_ribo_variation = names(rna_replicate_mean_weights)[which(rnacv/ribocv > 2)]
low_rna_ribo_variation = names(rna_replicate_mean_weights)[which(rnacv/ribocv < .5)]
high_rna_ribo_variation = hgnc_to_ensg_convert(high_rna_ribo_variation)
low_rna_ribo_variation = hgnc_to_ensg_convert(low_rna_ribo_variation)
addList(david, high_rna_ribo_variation, idType="ENSEMBL_GENE_ID", listName="HighRNARiboVariation", listType="Gene")
addList(david, low_rna_ribo_variation, idType="ENSEMBL_GENE_ID", listName="LowRNARiboVariation", listType="Gene")
setAnnotationCategories (david, c("GOTERM_CC_ALL", "GOTERM_BP_ALL", "GOTERM_MF_ALL", "KEGG_PATHWAY", "REACTOME_PATHWAY"))
AnnotCHART <- getFunctionalAnnotationChart(david, threshold=0.001, count=2L)
filter_by_fdr_fold_enrichment(AnnotCHART, .05,2)

# setCurrentGeneListPosition(david, 1)

# MAKE A HASHTABLE OF GO TERMS
# Comparative graph doesn't look very pretty 
# DAVID RESULTS and manual fisher calculation does not agree well
# Better strategy is output all cluster ids into separate files and run funcassociate
GO = hash()
fmat5 = matrix(nrow=2, ncol =2)
fmat8 = matrix(nrow=2, ncol =2)
dag = readLines('~/project/CORE_DATAFILES/FUNCASSOCIATE/Mixed_Effect_FuncAssociate/RESULTS/funcassociate_go_associations.txt')
dag = dag[-c(1:23)]
for (line in dag) {
  lineelements = unlist(strsplit(line,split="\t")[[1]])
  GOterm = lineelements[1]
  Genes  = unlist(strsplit(lineelements[3],split=" ")[[1]])
  fmat5[1,] = c( length (intersect (cluster5, Genes)) , length(setdiff(cluster5, Genes)))
  fmat5[2,] = c (length(setdiff( Genes, cluster5) ), 9651-length(union (cluster5, Genes)) ) 
  fmat8[1,] = c( length (intersect (cluster8, Genes)) , length(setdiff(cluster8, Genes)))
  fmat8[2,] = c (length(setdiff( Genes, cluster8) ), 9651-length(union (cluster8, Genes)) ) 
  GO[[GOterm]] = as.numeric(c ( fisher.test(fmat5)$estimate, fisher.test(fmat8)$estimate,  fisher.test(fmat5, conf.level = 0.9)$conf.int,  fisher.test(fmat8, conf.level = 0.9)$conf.int)   )
}

c5 = read.table('~/project/CORE_DATAFILES/GO_RESULTS/Absolute.SOM.Cluster_GOID_Only_5', stringsAsFactors=F)
c8 = read.table('~/project/CORE_DATAFILES/GO_RESULTS/Absolute.SOM.Cluster_GOID_Only_8', stringsAsFactors=F)

listofkeys = union( c5[,1], c8[,1])
barplotdata = matrix(nrow = 2, ncol = length(listofkeys))
colnames(barplotdata) <- rep("GO", length(listofkeys))
i = 1
for (key in listofkeys) {
  if (!is.null (GO[[key]])) {
    if(GO[[key]][1] / GO[[key]][2] > 5 | GO[[key]][1] / GO[[key]][2] < .2) {
      barplotdata[,i] = c(GO[[key]][1], GO[[key]][2])
      colnames(barplotdata)[i] = key
      i= i+1
    }
    else {
      barplotdata = barplotdata[,-i]    
    }
  }
  else {
    barplotdata = barplotdata[,-i]
  }
}
barplot(barplotdata, beside=T, horiz = T, cex.names = .5, col = c("blue", "magenta"))

