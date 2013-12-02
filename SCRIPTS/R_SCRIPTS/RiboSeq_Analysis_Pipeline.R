### Invocation
# R --no-save --no-restore --silent --args COL1 COL2 COL3 COL4 EX/CDS COND1 COND2 < RiboSeq_Analysis_Pipeline.R > OUTFILE
####
library("DEXSeq")
setwd('~/EMI_RIBOSEQ/')

arg <- commandArgs(TRUE)
cols_to_use <- as.numeric(c(arg[1], arg[2], arg[3], arg[4]))
if (arg[5] == "EX" ) {
   dat <- read.table('all_exon_counts')	
} else if ( arg[5] == "CDS" ){
   dat <- read.table('all_cds_counts')
} else {
     arg[5]
}
colnames(dat) <- c("Name", "C1_ribo", "C2_ribo", "C1_RNA", "C2_RNA", "O1_RNA", "O2_RNA","S1_RNA", "S2_RNA", "O1_ribo", "O2_ribo", "S1_ribo","S2_ribo")

# There are duplicates to be cleaned
data <- dat[!duplicated(dat[,c("Name")]), ]
names <- unlist(strsplit(as.character(data$Name), ":"))
gene_ids <- names[seq(1, length(names),2)]
exon_ids <- names[seq(2, length(names),2)]

count_data <- as.matrix(data[,cols_to_use])
condition <- c(arg[6], arg[6], arg[7], arg[7])
replicate <- c(1,2,1,2)
rownames(count_data) <- data$Name
design <- data.frame(condition=condition, replicate=replicate)
analysis <- newExonCountSet(
countData = count_data,
design=design,
geneIDs=gene_ids,
exonIDs=exon_ids)

file <- paste("10%FDR", arg[5], arg[6], arg[7], ".txt" , sep="_")
#file
normalized <- estimateSizeFactors(analysis)
sizeFactors(normalized)
#sizeFactors(normalized) <- c(1,1,1,1)
normalized <- estimateDispersions(normalized)
normalized <- fitDispersionFunction(normalized)
normalized <- testForDEU( normalized )
normalized <- estimatelog2FoldChanges( normalized )
res1 <- DEUresultTable(normalized)
head( res1 )
table ( res1$padjust < 0.1 )
table ( tapply( res1$padjust < 0.1, geneIDs(normalized), any ) )
changes <- which(res1$padjust < 0.1)
res1[changes,]
write.table(file = file , res1[changes,])
 
# #Cont to Over_ Translation_Efficiency
# # Need to sum over genes and get rid of 0s. 
# library("anota")
# dataT = as.matrix(data[,c(9,10,11,12,13,14)])
# dataP = as.matrix(data[,c(7,8,15,16,17,18)])
# grouping <- as.factor(gene_ids)
# summed_countsT <- matrix (nrow= nlevels(grouping), ncol=1)
# for (i in 1:dim(dataT)[2]) {
#     tmp <- tapply(dataT[,i],grouping, sum)
#     summed_countsT <- cbind(summed_countsT, tmp)
# }
# summed_countsT <- summed_countsT[,-1]

# summed_countsP <- matrix (nrow= nlevels(grouping), ncol=1)
# for (i in 1:dim(dataP)[2]) {
#     tmp <- tapply(dataP[,i],grouping, sum)
#     summed_countsP <- cbind(summed_countsP, tmp)
# }
# summed_countsP <- summed_countsP[,-1]

# all_counts <- cbind(summed_countsP, summed_countsT )
# logical <- all_counts > 10
# k <- c()
# for (i in 1:dim(logical)[1]) {
#     if(all(logical[i,])) {
#     			 k <- c(k, i)
# 			 }
# }

# PhenoVec <- c("cont", "cont", "over", "over", "siRNA", "siRNA")
# anotaQcOut <- anotaPerformQc (dataT=log10(summed_countsT[k,]), dataP= log10(summed_countsP[k,]), phenoVec=PhenoVec, onlyGroup=TRUE, useDfb=F)
#  anotaResidOut <- anotaResidOutlierTest(anotaQcObj=anotaQcOut)

# anotaSigOut <- anotaGetSigGenes(dataT=log10(summed_countsT[k,]), dataP= log10(summed_countsP[k,]), phenoVec=PhenoVec, anotaQcObj=anotaQcOut, useRVM=FALSE, useProgBar=FALSE)

# anotaSelected <- anotaPlotSigGenes(anotaSigObj=anotaSigOut, selContr=1, maxP=0.5, minSlope=(-0.5), maxSlope=1.5, selDeltaPT=0.5)



####### OLD CODE ##########
#cols_to_use
#is.vector(cols_to_use)
#exon_counts <- read.table('all_exon_counts')
#cds_counts <- read.table('all_cds_counts')
#colnames(exon_counts) <- colnames(cds_counts) <- c("Name", "C1_ribo", "C2_ribo", "C1_RNA", "C2_RNA", "O1_RNA", "O2_RNA","S1_RNA", "S2_RNA", "O1_ribo", "O2_ribo", "S1_ribo","S2_ribo")

# Coordinates if needed
#cds_coor <- read.table('cds_hg18_renamed.bed')
#exon_coor <- read.table('exons_hg18_renamed.bed')



# # RIBOSEQ
# # Control to OVER
# count_data <- as.matrix(data[,c(7,8,15,16)])
# # Control TO SIRNA
# count_data <- as.matrix(data[,c(7,8,17,18)])
# # OVER TO SIRNA
# count_data <- as.matrix(data[,c(15:18)])
# # RNASEQ
# # Control to OVER
# count_data <- as.matrix(data[,c(9,10,11,12)])
# # Control TO SIRNA
# count_data <- as.matrix(data[,c(9,10,13,14)])
# # OVER TO SIRNA
# count_data <- as.matrix(data[,c(11:14)])
# condition <- c("control", "control", "over", "over")
# condition <- c("control", "control", "siRNA", "siRNA")
# condition <- c( "over", "over", "siRNA", "siRNA")

