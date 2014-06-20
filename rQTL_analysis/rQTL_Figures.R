### CC rQTL FINAL FIGURES

## PLOT OF BETAS RNA_RIBO; RIBO_PROT
# Q1: Should we be plotting p < .05 in both or only ribo
# Q2: Should we be using the EMP2? Or P? 
p <- read.table("~/project/rQTL_analysis/qtlMapping.protein.21YRI.txt",header=T,as.is=T)
r2 <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.21YRI.txt",header=T,as.is=T)
rn1 <- read.table("~/project/rQTL_analysis/qtlMapping.rna1.21YRI.txt",header=T,as.is=T)
te <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.CovarRNA1.21YRI.txt",header=T,as.is=T)

rna_ribo =  merge(rn1,r2,by="SNP")
cor(rna_ribo$BETA.x,rna_ribo$BETA.y,method="spearman") #0.35
plot(rna_ribo$BETA.x,rna_ribo$BETA.y)
# Nominal Ribo QTL P < .05
rna_ribo_ribosig <- rna_ribo[rna_ribo$EMP2.y<0.01,]
cor(rna_ribo_ribosig$BETA.x,rna_ribo_ribosig$BETA.y,method="spearman") #0.81
pdf (file = "Comparison_of_Betas", width= 4, height = 4)
plot(rna_ribo_ribosig$BETA.x,rna_ribo_ribosig$BETA.y, xlab="RNA Beta", ylab= "Ribo Beta",
     tck=.02, pch=19, cex=.65)
abline(0,1)
abline(v=0,h=0)
dev.off()

# prot_ribo = merge( p, r2, by ="SNP")
# cor(prot_ribo$BETA.x,prot_ribo$BETA.y,method="spearman") #0.29
# plot(prot_ribo$BETA.x,prot_ribo$BETA.y)
# # Nominal Ribo QTL P < .05
# prot_ribo_ribosig <- prot_ribo[prot_ribo$EMP2.y<0.01,]
# cor(prot_ribo_ribosig$BETA.x,prot_ribo_ribosig$BETA.y,method="spearman") #0.64
# plot(prot_ribo_ribosig$BETA.x,prot_ribo_ribosig$BETA.y,  ylab= "Ribo Beta", xlab= "Protein Beta",
#      tck=.02, pch=19, cex=.65,  xlim = c(-1, 1), ylim= c(-1,1))
# abline(0,1)
# abline(v=0,h=0)

# PLOT OF RNA P-VALs for SIGNIFICANT RIBO QTLS at FDR 30% -- BEST HITS
# PLOT OF PROTEIN P_VALs
rb <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.21YRI.bestHits.txt",header=T,as.is=T)
# Use p and rn1
sign <- rb$SNP[rb$EMP2<0.01]
hist(rn1$P[is.element(rn1$SNP,sign)],nclass=50,main="nom.pValue as eQTL \nof SNP with p.emp2<0.01 (n=159, fdr=45%) for rQTL",xlab="nominal p as eQTL ")
hist(p$P[is.element(p$SNP,sign)],nclass=50,main="nom.pValue as pQTL \nof SNP with p.emp2<0.01 (n=159, fdr=45%) for rQTL",xlab="nominal p as pQTL ")

#results:  fdr10      fdr30    fdr50
# protein    4         29        48    (logRatio)

# ribo1      0          55        227  ()
# ribo2      0          67        256  (qn)
# ribo3      0          48        156  (sv20)
# ribo4      15         60        164  (sv3) 

# rna1       22         128       488  (qn)
# rna2       15         100       513   (sv20)
# rna3       31         164       398   (sv3)
# CHECK P-VAL FDR CORRESPONDENCE

# QQPLOT FOR VARIANT POSITION AND TYPE FOR RNA AND RIBO
# USE all; p, rn1, r2
annot <- read.table("~/project/rQTL_analysis/ALL_SNPs_Annotation_v1",header=F,as.is=T)
annot$snpID <- paste(annot$V1,annot$V2,sep=";")

makeSNPid <- function(x){
  k <- unlist(strsplit(x,";"))
  res <- paste(k[1],k[2],sep=";")
  return(res)
}


p$snpID <- unlist(lapply(p$SNP,makeSNPid))
r2$snpID <- unlist(lapply(r2$SNP,makeSNPid))
rn1$snpID <- unlist(lapply(rn1$SNP,makeSNPid))
te$snpID <- unlist(lapply(te$SNP,makeSNPid))

p <- merge(p,annot,by="snpID")
r <- merge(r2,annot,by="snpID")
rna <- merge(rn1,annot,by="snpID")
te <- merge(te,annot,by="snpID")

r <- r[order(r$P,decreasing=T),]
p <- p[order(p$P,decreasing=T),]
rna <- rna[order(rna$P,decreasing=T),]
te <- te[order(te$P, decreasing=T), ]

annotated_start = annot[grep("Annotated_Start",annot$V3),]$snpID
nonsyn = annot[grep("Nonsynonymous",annot$V3),]$snpID
syn = annot[grep("Synonymous",annot$V3),]$snpID
utr3 = annot[grep("3UTR",annot$V3),]$snpID
utr5 = annot[grep("5UTR",annot$V3),]$snpID

plot_qq = function (x) { 
obs_non <- sort(-log10(x$P[x$snpID %in% nonsyn]))
exp_non <- sort(-log10(seq(1:length(x$P[x$snpID %in% nonsyn]))/length(x$P[x$snpID %in% nonsyn])))
obs_syn <- sort(-log10(x$P[x$snpID %in% syn]))
exp_syn <- sort(-log10(seq(1:length(x$P[x$snpID %in% syn]))/length(x$P[x$snpID %in% syn])))
obs_utr3 <- sort(-log10(x$P[x$snpID %in% utr3]))
exp_utr3 <- sort(-log10(seq(1:length(x$P[x$snpID %in% utr3]))/length(x$P[x$snpID %in% utr3])))
obs_utr5 <- sort(-log10(x$P[x$snpID %in% utr5]))
exp_utr5 <- sort(-log10(seq(1:length(x$P[x$snpID %in% utr5]))/length(x$P[x$snpID %in% utr5])))

plot(exp_non,obs_non,pch=".",cex=4,col="red", main = deparse(substitute(x)))
points(exp_syn,obs_syn,pch=".",cex=4,col="blue")
points(exp_utr3,obs_utr3,pch=".",cex=4,col="green")
points(exp_utr5,obs_utr5,pch=".",cex=4,col="purple")
abline(0,1)
}

plot_qq(r)
plot_qq(p)
plot_qq(rna)
plot_qq(te)

