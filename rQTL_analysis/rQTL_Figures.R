### CC rQTL FINAL FIGURES
# ece1_qtldata
# [1] 6.191203 5.776160 3.664089 5.431239 6.733963 3.894626 6.664624 4.975776 6.320706 5.368113 6.378847
# [12] 6.239615 5.565489 6.185235 5.713562 3.907457 5.531433 4.164887 4.991089 5.183971 5.020389 5.174953
# [23] 6.463338 3.955402 6.195927 5.787483 6.559498 4.862363 4.375954 6.993334
# > ece1_genotypes
# [1] 1 0 0 0 0 1 0 1 0 1 1 1 1 1 1 2 0 2 0 1 1 1 1 1 1 1 1 1 1 1
# > ece1_ids
# [1] NA10847 NA12878 NA12890 NA12891 NA12892 NA18489 NA18499 NA18501 NA18502 NA18504 NA18505 NA18507
# [13] NA18511 NA18519 NA18522 NA18523 NA18526 NA18870 NA18951 NA19098 NA19099 NA19137 NA19138 NA19139
# [25] NA19193 NA19200 NA19201 NA19238 NA19239 NA19240
# yri = c(6:16, 18, 20:30)
# > ece1_rna
# ENST  GM10847  GM12878  GM12890  GM12891  GM12892  GM18489  GM18499  GM18501  GM18502
# 2401 ENST00000264205.6 4.614093 5.267574 2.500083 4.538748 5.491529 2.573254 5.816071 2.749624 5.601134
# GM18504  GM18505  GM18507 GM18511  GM18519  GM18522  GM18523  GM18526  GM18870  GM18951  GM19098
# 2401 6.026954 5.230803 4.533472 3.91226 5.313835 4.583976 4.091149 4.132438 3.787842 4.585178 4.065722
# GM19099  GM19137  GM19138  GM19193  GM19200  GM19201  GM19238 GM19239  GM19240
# 2401 4.241511 4.271718 5.417058 4.763614 4.928664 5.288601 2.214951 3.23084 5.583783
# boxplot(ece1_qtldata[yri] ~ece1_genotypes[yri])
# summary(lm(ece1_qtldata[yri] ~ece1_genotypes[yri]))

# > summary(lm(ece1_exp_rna[yri_rna] ~ ece1_genotypes_rna[yri_rna]))
# 
# Call:
#   lm(formula = ece1_exp_rna[yri_rna] ~ ece1_genotypes_rna[yri_rna])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.2499 -0.3552  0.2296  0.7022  1.5621 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   5.3494     0.5590    9.57 6.59e-09 ***
#   ece1_genotypes_rna[yri_rna]  -0.8846     0.5142   -1.72    0.101    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.028 on 20 degrees of freedom
# Multiple R-squared:  0.1289,  Adjusted R-squared:  0.08535 
# F-statistic:  2.96 on 1 and 20 DF,  p-value: 0.1008

### rQTL Information for NTPCR is robust to YRI only
# dn52empp:rQTL_analysis cancenik$ grep 'ENST00000366628' qtlMapping.ribo2.21YRI.txt 
# ENST00000366628.4;87;T;C  ENST00000366628.4	1	87	C	ADD	21	0.4362	3.348	0.003383	0.003	0.003	C	T	0.3571
# dn52empp:rQTL_analysis cancenik$ grep 'ENST00000366628' qtlMapping.rna1.21YRI.txt 
# ENST00000366628.4;87;T;C	ENST00000366628.4	1	87	C	ADD	21	0.1836	1.767	0.09336	0.09299	0.09299	C	T	0.3571
# Similarly it is a protein QTL
# dn52empp:rQTL_analysis cancenik$ grep 'ENST00000366628' qtlMapping.protein.21YRI.txt 
# ENST00000366628.4;87;T;C  ENST00000366628.4	1	87	C	ADD	21	0.1595	2.492	0.0221	0.0219	0.0219	C	T	0.3571


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
# Ribo QTL P < .05
rna_ribo_ribosig <- rna_ribo[rna_ribo$EMP2.y<0.003,]
cor(rna_ribo_ribosig$BETA.x,rna_ribo_ribosig$BETA.y,method="spearman") #0.86

# Color the points by rna EMP2
pval_colors <- rep("Black", times = length(rna_ribo_ribosig$BETA.x ))
pval_colors[rna_ribo_ribosig$P.x < .05] <- "Red"

pdf (file = "~/Google_Drive/Manuscript Figures/rQTLs/Comparison_of_Betas.pdf",
     width= 4, height = 4)
plot(rna_ribo_ribosig$BETA.x,rna_ribo_ribosig$BETA.y, xlab="RNA Beta", ylab= "Ribo Beta",
     tck=.02, pch=19, cex=.65, col = pval_colors)
abline(0,1)
abline(v=0,h=0)
dev.off()

# We will probably not going to compare protein Betas
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
# rb$EMP2<0.0029 ~ fdr 30
sign <- rb$SNP[rb$EMP2<0.0029]
pdf (file = "~/Google_Drive/Manuscript Figures/roQTLs//PVals_as_eQTL.pdf",
     width= 6, height = 5)
hist(rn1$P[is.element(rn1$SNP,sign)],xlim= c(0,1), nclass=50,main="P as eQTL of SNP with\n rb$EMP2<0.0029 (n=67, fdr=30%)",xlab="P as eQTL ")
dev.off()
pdf (file = "~/Google_Drive/Manuscript Figures/roQTLs/PVals_as_pQTL.pdf",
     width= 6, height = 5)
hist(p$P[is.element(p$SNP,sign)],nclass=50, xlim= c(0,1), main="P as pQTL of SNP with\n rb$EMP2<0.0029 (n=16, fdr=30%) for rQTL",xlab="P as pQTL ")
dev.off()

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

r <- r[order(r$EMP2,decreasing=T),]
p <- p[order(p$EMP2,decreasing=T),]
rna <- rna[order(rna$EMP2,decreasing=T),]
te <- te[order(te$EMP2, decreasing=T), ]

annotated_start = annot[grep("Annotated_Start",annot$V3),]$snpID
nonsyn = annot[grep("Nonsynonymous",annot$V3),]$snpID
syn = annot[grep("Synonymous",annot$V3),]$snpID
utr3 = annot[grep("3UTR",annot$V3),]$snpID
utr5 = annot[grep("5UTR",annot$V3),]$snpID

plot_qq = function (x) { 
obs_non <- sort(-log10(x$EMP2[x$snpID %in% nonsyn]))
exp_non <- sort(-log10(seq(1:length(x$EMP2[x$snpID %in% nonsyn]))/length(x$EMP2[x$snpID %in% nonsyn])))
obs_syn <- sort(-log10(x$EMP2[x$snpID %in% syn]))
exp_syn <- sort(-log10(seq(1:length(x$EMP2[x$snpID %in% syn]))/length(x$EMP2[x$snpID %in% syn])))
obs_utr3 <- sort(-log10(x$EMP2[x$snpID %in% utr3]))
exp_utr3 <- sort(-log10(seq(1:length(x$EMP2[x$snpID %in% utr3]))/length(x$EMP2[x$snpID %in% utr3])))
obs_utr5 <- sort(-log10(x$EMP2[x$snpID %in% utr5]))
exp_utr5 <- sort(-log10(seq(1:length(x$EMP2[x$snpID %in% utr5]))/length(x$EMP2[x$snpID %in% utr5])))

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

