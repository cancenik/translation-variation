options(scipen=100)  
# 675 line cont
annot <- read.table("~/project/rQTL_analysis/ALL_SNPs_Annotation_v1",header=F,as.is=T)
dim(annot)
annot[1:5,]

annot$snpID <- paste(annot$V1,annot$V2,sep=";")
p <- read.table("~/project/rQTL_analysis/qtlMapping.protein.21YRI.txt",header=T,as.is=T)
r <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.21YRI.txt",header=T,as.is=T)
rna <- read.table("~/project/rQTL_analysis/qtlMapping.rna1.21YRI.txt",header=T,as.is=T)
te <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.CovarRNA1.21YRI.txt",header=T,as.is=T)

p[1:5,]
table(is.na(r$P))
table(is.na(p$P))
table(is.na(rna$P))

#make the same SNP ID:

makeSNPid <- function(x){
	k <- unlist(strsplit(x,";"))
	res <- paste(k[1],k[2],sep=";")
	return(res)
}

p$snpID <- unlist(lapply(p$SNP,makeSNPid))
r$snpID <- unlist(lapply(r$SNP,makeSNPid))
rna$snpID <- unlist(lapply(rna$SNP,makeSNPid))
te$snpID <- unlist(lapply(te$SNP,makeSNPid))

p <- merge(p,annot,by="snpID")
dim(annot)
dim(p)
p[1:5,]
r <- merge(r,annot,by="snpID")
rna <- merge(rna,annot,by="snpID")
te <- merge(te,annot,by="snpID")

#redo annotation:
splitAnnot <- function(x){
	k <- unlist(strsplit(x,";"))
	return(k[length(k)])
}

p$tag <-  unlist(lapply(p$V3,splitAnnot))
r$tag <-  unlist(lapply(r$V3,splitAnnot))
rna$tag <-  unlist(lapply(rna$V3,splitAnnot))
te$tag <-  unlist(lapply(te$V3,splitAnnot))

table(p$tag)
p[p$tag=="A_G",]
annot[grep("Annotated_Start",annot$V3),]
p[grep("Annotated_Start",p$V3),]
annot[grep("ENST00000390556.2;1452",annot$snpID),]
p[1:5,]
boxplot(abs(p$STAT)~p$tag,ylim=c(0,1))



table(r$tag)
#simplify:
r <- r[is.element(r$tag,c("3UTR","5UTR","Nonsynonymous","Synonymous")),]
rna <- rna[is.element(rna$tag,c("3UTR","5UTR","Nonsynonymous","Synonymous")),]
p <- p[is.element(p$tag,c("3UTR","5UTR","Nonsynonymous","Synonymous")),]
te <- te[is.element(te$tag,c("3UTR","5UTR","Nonsynonymous","Synonymous")),]


boxplot(abs(r$STAT)~r$tag,ylim=c(0,1))
boxplot(abs(r$STAT)~r$tag)

summary(abs(r$STAT[r$tag=="Nonsynonymous"]))
summary(abs(r$STAT[r$tag=="Synonymous"]))
summary(abs(r$STAT[r$tag=="3UTR"]))
summary(abs(r$STAT[r$tag=="5UTR"]))
# > summary(abs(r$STAT[r$tag=="Nonsynonymous"]))
     # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
 # 0.000922  0.350600  0.747400  0.932500  1.316000 10.970000 
# > summary(abs(r$STAT[r$tag=="Synonymous"]))
     # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
 # 0.000028  0.343100  0.739800  0.904600  1.274000 10.970000 
# > summary(abs(r$STAT[r$tag=="3UTR"]))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000039 0.343400 0.724000 0.884500 1.246000 9.534000 
# > summary(abs(r$STAT[r$tag=="5UTR"]))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000776 0.344200 0.746100 0.929000 1.336000 5.445000 

summary(abs(p$STAT[p$tag=="Nonsynonymous"]))
summary(abs(p$STAT[p$tag=="Synonymous"]))
summary(abs(p$STAT[p$tag=="3UTR"]))
summary(abs(p$STAT[p$tag=="5UTR"]))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000642 0.353400 0.770600 0.955900 1.270000 7.173000 
# > summary(abs(p$STAT[p$tag=="Synonymous"]))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000324 0.333000 0.710400 0.884000 1.202000 7.286000 
# > summary(abs(p$STAT[p$tag=="3UTR"]))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000716 0.328900 0.724200 0.882300 1.222000 7.286000 
# > summary(abs(p$STAT[p$tag=="5UTR"]))
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000158 0.334000 0.717000 0.873300 1.230000 5.027000 

summary(abs(r$STAT[grep("Annotated_Start",r$V3)]))
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 0.6537  0.9862  1.3860  1.3460  1.6480  2.0630 
#these are very significant ->but only 6...
length(r$STAT[grep("Annotated_Start",r$V3)]) #6

summary(abs(r$STAT[grep("Annotated_End",r$V3)]))
length(r$STAT[grep("Annotated_End",r$V3)]) #0


summary(abs(p$STAT[grep("Annotated_Start",p$V3)]))
length(p$STAT[grep("Annotated_Start",p$V3)]) #1
p[grep("Annotated_Start",p$V3),]
r[r$snpID=="ENST00000359678.5;297",]

summary(abs(r$STAT[r$tag=="Nonsynonymous"]))
summary(abs(r$STAT[r$tag=="Nonsynonymous"]))



########################try qqplots:
r <- r[order(r$P,decreasing=T),]
p <- p[order(p$P,decreasing=T),]
rna <- rna[order(rna$P,decreasing=T),]

obsAll <- sort(-log10(r$P))
obsAll[1:100]
table(obsAll==0)
hist(r$P)
obsAll <- sort(-log10(r$P))
exp <- sort(-log10(seq(1:length(r$P))/length(r$P)))
hist(seq(1:length(r$P))/length(r$P))

min(r$P)
min(rna$P)
dim(rna)
dim(r)
dim(p)



par(mfrow=c(2,2))
obsAll <- sort(-log10(r$P))
exp <- sort(-log10(seq(1:length(r$P))/length(r$P)))

plot(exp,obsAll,pch=".",cex=4,main="riboSeq")
abline(0,1)

obs <- r$P[r$tag=="Nonsynonymous"]
exp <- seq(1:length(r$P[r$tag=="Nonsynonymous"]))/length(r$P[r$tag=="Nonsynonymous"])
obs <- sort(-log10(r$P[r$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="Nonsynonymous"]))/length(r$P[r$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(r$P[r$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="Synonymous"]))/length(r$P[r$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(r$P[r$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="5UTR"]))/length(r$P[r$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(r$P[r$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="3UTR"]))/length(r$P[r$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")




#####rna
obsAll <- sort(-log10(rna$P))
exp <- sort(-log10(seq(1:length(rna$P))/length(rna$P)))
plot(exp,obsAll,pch=".",cex=4,main="rna")
abline(0,1)
obs <- sort(-log10(rna$P[rna$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(rna$P[rna$tag=="Nonsynonymous"]))/length(rna$P[rna$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(rna$P[rna$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(rna$P[rna$tag=="Synonymous"]))/length(rna$P[rna$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(rna$P[rna$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(rna$P[rna$tag=="5UTR"]))/length(rna$P[rna$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(rna$P[rna$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(rna$P[rna$tag=="3UTR"]))/length(rna$P[rna$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")



#protein 
obsAll <- sort(-log10(p$P))
exp <- sort(-log10(seq(1:length(p$P))/length(p$P)))
plot(exp,obsAll,pch=".",cex=4,main="protein")
abline(0,1)
obs <- sort(-log10(p$P[p$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(p$P[p$tag=="Nonsynonymous"]))/length(p$P[p$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(p$P[p$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(p$P[p$tag=="Synonymous"]))/length(p$P[p$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(p$P[p$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(p$P[p$tag=="5UTR"]))/length(p$P[p$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(p$P[p$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(p$P[p$tag=="3UTR"]))/length(p$P[p$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")


obsAll <- sort(-log10(te$P))
exp <- sort(-log10(seq(1:length(te$P))/length(te$P)))
plot(exp,obsAll,pch=".",cex=4,main="TE")
abline(0,1)
obs <- sort(-log10(te$P[te$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(te$P[te$tag=="Nonsynonymous"]))/length(te$P[te$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(te$P[te$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(te$P[te$tag=="Synonymous"]))/length(te$P[te$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(te$P[te$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(te$P[te$tag=="5UTR"]))/length(te$P[te$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(te$P[te$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(te$P[te$tag=="3UTR"]))/length(te$P[te$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")




#protein not 6

p6 <- p[p$CHR!=6,]
dim(p6)

obsAll <- sort(-log10(p6$P))
exp <- sort(-log10(seq(1:length(p6$P))/length(p6$P)))
plot(exp,obsAll,pch=".",cex=4,main="Protein, no chr6")
abline(0,1)
obs <- sort(-log10(p6$P[p6$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="Nonsynonymous"]))/length(p6$P[p6$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(p6$P[p6$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="Synonymous"]))/length(p6$P[p6$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(p6$P[p6$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="5UTR"]))/length(p6$P[p$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(p6$P[p6$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="3UTR"]))/length(p6$P[p6$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")

legend("bottomright",fill=c("red","blue","green","purple"),c("NonSyn","Syn","5UTR","3UTR"))

p[1:5,]
ns <- p[p$tag=="Nonsynonymous",]
ns[1:5,]
dim(ns) #947
length(unique(ns$ENST)) #566
hist(ns$P)
s <- p[p$tag=="Synonymous",]
dim(s) #2007
hist(s$P)

#less than 0.01
ns <- ns[ns$P<0.1,]
dim(ns) #130
length(unique(ns$ENST)) #566 -> 82ENST
ns[1:5,]
table(ns$CHR)
 # 1  2  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 22 23 25 
# 15 13  9  2  7 28  5  2  8  2  3  6  9  2  3  6  3  4  1  1  1 
#clearl an excess on chr 6! hla issue
table(p$tag)
length(unique(p$ENST[p$tag=="Nonsynonymous"])) #612



p6 <- p[p$CHR!=6,]
dim(p6)

obsAll <- sort(-log10(p6$P))
exp <- sort(-log10(seq(1:length(p6$P))/length(p6$P)))
plot(exp,obsAll,pch=".",cex=4)
abline(0,1)
obs <- sort(-log10(p6$P[p6$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="Nonsynonymous"]))/length(p6$P[p6$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(p6$P[p6$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="Synonymous"]))/length(p6$P[p6$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(p6$P[p6$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="5UTR"]))/length(p6$P[p$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(p6$P[p6$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(p6$P[p6$tag=="3UTR"]))/length(p6$P[p6$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")



#would have to check that the peptides did not overlap with the 1000 SNPs
#take the top examle
pns <- p[p$tag=="Nonsynonymous",]
pns <- pns[order(pns$P),]
pns[1:20,]
###make it for the riboseq 2050 genes
obsAll <- sort(-log10(r$P))
exp <- sort(-log10(seq(1:length(r$P))/length(r$P)))

par(mfrow=c(1,2))
plot(exp,obsAll,pch=".",cex=4)
abline(0,1)
obs <- sort(-log10(r$P[r$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="Nonsynonymous"]))/length(r$P[r$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(r$P[r$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="Synonymous"]))/length(r$P[r$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(r$P[r$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="5UTR"]))/length(r$P[r$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(r$P[r$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="3UTR"]))/length(r$P[r$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")


r <- r[is.element(r$ENST,p$ENST),]
obsAll <- sort(-log10(r$P))
exp <- sort(-log10(seq(1:length(r$P))/length(r$P)))

plot(exp,obsAll,pch=".",cex=4)
abline(0,1)
obs <- sort(-log10(r$P[r$tag=="Nonsynonymous"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="Nonsynonymous"]))/length(r$P[r$tag=="Nonsynonymous"])))
points(exp,obs,pch=".",cex=4,col="red")
obs <- sort(-log10(r$P[r$tag=="Synonymous"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="Synonymous"]))/length(r$P[r$tag=="Synonymous"])))
points(exp,obs,pch=".",cex=4,col="blue")

obs <- sort(-log10(r$P[r$tag=="5UTR"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="5UTR"]))/length(r$P[r$tag=="5UTR"])))
points(exp,obs,pch=".",cex=4,col="green")

obs <- sort(-log10(r$P[r$tag=="3UTR"]))
exp <- sort(-log10(seq(1:length(r$P[r$tag=="3UTR"]))/length(r$P[r$tag=="3UTR"])))
points(exp,obs,pch=".",cex=4,col="purple")




#################################
# COMPARING BETAS FOR THE SIGNIFICANT OR LOW P_VALUE VARIANTS
p <- read.table("~/project/rQTL_analysis/qtlMapping.protein.21YRI.txt",header=T,as.is=T)
r <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.21YRI.txt",header=T,as.is=T)
p[1:5,]
p <- merge(p,r,by="SNP")

plot(p$BETA.x,p$BETA.y)
abline(v=0,h=0)
abline(0,1)
plot(p$STAT.x,p$STAT.y)
abline(v=0,h=0)
abline(0,1)
cor(p$BETA.x,p$BETA.y,method="spearman")
cor(p$BETA.x,p$BETA.y,method="spearman")
cor(p$STAT.x,p$STAT.y,method="spearman")


#rm chr6:
class(p$CHR)
p <- p[p$CHR.x!=6,]
dim(p)
plot(p$BETA.x,p$BETA.y,xlim=c(-2,2),ylim=c(-2,2))
abline(v=0,h=0)
abline(0,1)
abline(0,1)


par(mfrow=c(2,1))
hist(p$BETA.x,xlim=c(-2,2),nclass=100)
abline(v=c(-0.25,0.25))
hist(p$BETA.y,xlim=c(-2,2),nclass=100)
abline(v=c(-0.25,0.25))


#look at it for p<0.05
pProt0.05 <- p[p$P.x<0.05,]
dim(pProt0.05) #462
plot(pProt0.05$BETA.x,pProt0.05$BETA.y,xlim=c(-2,2),ylim=c(-2,2))
abline(v=0,h=0)
abline(0,1)
summary(lm(pProt0.05$BETA.y~pProt0.05$BETA.x))
#est 0.911791 
summary(lm(pProt0.05$BETA.x~pProt0.05$BETA.y))

summary(pProt0.05$BETA.y) #ribo
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.62700 -0.15210 -0.03011 -0.02454  0.10830  1.98000 

summary(pProt0.05$BETA.x) #protein
    # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.477500 -0.126200 -0.050080 -0.008039  0.120300  0.608600 



pProt0.05 <- p[p$P.y<0.05,]
dim(pProt0.05) #442
plot(pProt0.05$BETA.x,pProt0.05$BETA.y,xlim=c(-2,2),ylim=c(-2,2))
abline(v=0,h=0)
abline(0,1)
abline(0,1)
summary(lm(pProt0.05$BETA.y~pProt0.05$BETA.x))
#est2.41410
summary(lm(pProt0.05$BETA.x~pProt0.05$BETA.y))
#0.1468178
summary(pProt0.05$BETA.y) #ribo
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.62700 -0.31200 -0.16880 -0.02367  0.28560  1.98000 

summary(pProt0.05$BETA.x) #protein
     # Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.4376000 -0.0482400 -0.0006661 -0.0028960  0.0464700  0.5379000 


pProt0.05 <- p[p$P.y<0.05&p$P.x<0.05,]
dim(pProt0.05) #63
par(mfrow=c(1,1))
plot(pProt0.05$BETA.x,pProt0.05$BETA.y,xlim=c(-2,2),ylim=c(-2,2),xlab="beta protein log2 ratio", ylab="riboseq (1) beta - log counts scale",main="comparing genetic effect magnitude\n nom.p <0.05 in riboseq and protein")
abline(v=0,h=0)
abline(0,1)





summary(lm(pProt0.05$BETA.y~pProt0.05$BETA.x))
summary(lm(pProt0.05$BETA.x~pProt0.05$BETA.y))
summary(pProt0.05$BETA.y) #ribo
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.62700 -0.39360 -0.23260 -0.03702  0.25770  1.98000 

summary(pProt0.05$BETA.x) #protein
    # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.43760 -0.16400 -0.07704 -0.02098  0.13500  0.53790 


pProt0.05[1:5,]
length(unique(pProt0.05$ENST.x)) #37


#################################################################
best <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.21YRI.bestHits.txt",header=T,as.is=T)
best[1:5,]
	
source("http://www.bioconductor.org/biocLite.R")
#biocLite("qvalue")
citation("qvalue")
library(qvalue)

p <- best$EMP2
hist(p,nclass=1000)
min(p)
qobj <- qvalue(p,fdr.level=0.1)
qplot(qobj)
table(qobj$significant) #0 at fdr 10%, 12 at 20%, 55 at 30% , 227 at FDR50%

par(mfrow=c(1,2))
hist(p,nclass=1000)
hist(p[best$CHR!=6],nclass=1000)

best <- read.table("~/project/rQTL_analysis/qtlMapping.protein.21YRI.bestHits.txt",header=T,as.is=T)
best[1:5,]

p <- best$EMP2
hist(p,nclass=1000)
min(p)
qobj <- qvalue(p,fdr.level=0.5)
qplot(qobj)
table(qobj$significant) #48 at FDR50%

par(mfrow=c(1,2))
hist(p,nclass=1000)
hist(p[best$CHR!=6],nclass=1000)

plot(-log10(best$P),-log10(best$EMP1))
plot(-log10(best$EMP2),-log10(best$EMP1))

#############################
#FDR
###############################
best <- read.table("~/project/rQTL_analysis/qtlMapping.protein.21YRI.bestHits.txt",header=T,as.is=T)
best[1:5,]
library(qvalue)
p <- best$EMP2
hist(p,nclass=1000)
min(p)
qobj <- qvalue(p,fdr.level=0.5)
qplot(qobj)
table(qobj$significant) #48 at FDR50%
qobj <- qvalue(p,fdr.level=0.3)
qplot(qobj)
table(qobj$significant) #29 at FDR30%
qobj <- qvalue(p,fdr.level=0.1)
table(qobj$significant) #4 at FDR10%

par(mfrow=c(1,2))
hist(p,nclass=1000)
hist(p[best$CHR!=6],nclass=1000)


bestR2 <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.21YRI.bestHits.txt",header=T,as.is=T)
bestR2[1:5,]
pR2 <- bestR2$EMP2
hist(pR2,nclass=1000)
min(pR2)
qobj <- qvalue(pR2,fdr.level=0.5)
qplot(qobj)
table(qobj$significant) #256 at FDR50%
qobj <- qvalue(pR2,fdr.level=0.3)
qplot(qobj)
table(qobj$significant) #67 at FDR30%
qobj <- qvalue(pR2,fdr.level=0.1)
qplot(qobj)
table(qobj$significant) #0 at FDR10%

par(mfrow=c(1,2))
hist(pR2,nclass=1000)
hist(pR2[best$CHR!=6],nclass=1000)

bestRNA1 <- read.table("~/project/rQTL_analysis/qtlMapping.rna1.21YRI.bestHits.txt",header=T,as.is=T)
bestRNA1[1:5,]
pRNA1 <- bestRNA1$EMP2
hist(pRNA1,nclass=1000)
min(pRNA1)
qobj <- qvalue(pRNA1,fdr.level=0.5)
qplot(qobj)
table(qobj$significant) #488 at FDR50%
qobj <- qvalue(pRNA1,fdr.level=0.3)
qplot(qobj)
table(qobj$significant) #128 at FDR30%
qobj <- qvalue(pRNA1,fdr.level=0.1)
table(qobj$significant) #22 at FDR10%

par(mfrow=c(1,2))
hist(pRNA1,nclass=1000)
hist(pRNA1[bestRNA1$CHR!=6],nclass=1000)


#this is ribo2~geno+rna1
bestR5 <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.CovarRNA1.21YRI.bestHits.txt",header=T,as.is=T)
bestR5[1:5,]
pR5 <- bestR5$EMP2
hist(pR5,nclass=1000)
min(pR5)
qobj <- qvalue(pR5,fdr.level=0.5)
qplot(qobj)
table(qobj$significant) #10 at FDR50%
qobj <- qvalue(pR5,fdr.level=0.3)
qplot(qobj)
table(qobj$significant) #3 at FDR30%
qobj <- qvalue(pR5,fdr.level=0.1)
qplot(qobj)
table(qobj$significant) #0 at FDR10%


#results:  fdr10      fdr30    fdr50
# protein    4         29        48    (logRatio)

# ribo1      0          55        227  ()
# ribo2      0          67        256  (qn)
# ribo3      0          48        156  (sv20)
# ribo4      15         60        164  (sv3) 


# rna1       22         128       488  (qn)
# rna2       15         100       513   (sv20)
# rna3       31         164       398   (sv3)


#ribo2~g+rna1  0         3         10



############
#looking at beta values, protein rna riboseq
###########

p <- read.table("~/project/rQTL_analysis/qtlMapping.protein.21YRI.txt",header=T,as.is=T)
r2 <- read.table("~/project/rQTL_analysis/qtlMapping.ribo2.21YRI.txt",header=T,as.is=T)
rn1 <- read.table("~/project/rQTL_analysis/qtlMapping.rna1.21YRI.txt",header=T,as.is=T)
hist(rn1$P,ylim=c(0,3000))

test <- merge(rn1,r2,by="SNP")
dim(test)
cor(test$BETA.x,test$BETA.y,method="spearman") #0.35
plot(test$BETA.x,test$BETA.y)
#pretty correlated!
plot(-log10(test$P.x),-log10(test$P.y))

test <- merge(rn1,p,by="SNP")
cor(test$BETA.x,test$BETA.y,method="spearman") #0.30
plot(test$BETA.x,test$BETA.y)
#rm the chr6:
test <- test[test$CHR.x!=6,]
cor(test$BETA.x,test$BETA.y,method="spearman") #0.28
plot(test$BETA.x,test$BETA.y)
#pretty correlated!
plot(-log10(test$P.x),-log10(test$P.y))

test <- merge(rn1,r2,by="SNP")
cor(test$BETA.x,test$BETA.y,method="spearman") #0.35
plot(test$BETA.x,test$BETA.y)
testS <- test[test$P.x<0.05&test$P.y<0.05,]
dim(testS)
cor(testS$BETA.x,testS$BETA.y,method="spearman") #0.92

par(mfrow=c(1,1))
plot(testS$BETA.x,testS$BETA.y,main="beta for ribo and rna\n for associations with p<0.05 in both",xlab="beta rna",ylab="beta ribo")
abline(0,1)
abline(v=0,h=0)

#look at them for the proteins:
testS[1:5,]
testS$sign <- testS$BETA.x*testS$BETA.y
opp <- testS$SNP[testS$sign<0]


test <- merge(rn1,p,by="SNP")
cor(test$BETA.x,test$BETA.y,method="spearman") #0.35
plot(test$BETA.x,test$BETA.y)
points(test$BETA.x[is.element(test$SNP,opp)],test$BETA.y[is.element(test$SNP,opp)],col="red") 

testS <- test[test$P.x<0.05&test$P.y<0.05,]
testS <- testS[testS$CHR.x!=6,]

cor(testS$BETA.x,testS$BETA.y,method="spearman") #0.92

plot(testS$BETA.x,testS$BETA.y,main="beta for protein and rna\n for associations with p<0.05 in both",xlab="beta rna",ylab="beta protein")
abline(0,1)
abline(v=0,h=0)
points(testS$BETA.x[is.element(testS$SNP,opp)],testS$BETA.y[is.element(testS$SNP,opp)],col="red") 



test <- merge(r2,p,by="SNP")
dim(test)
cor(test$BETA.x,test$BETA.y,method="spearman") #0.35
plot(test$BETA.x,test$BETA.y)
points(test$BETA.x[is.element(test$SNP,opp)],test$BETA.y[is.element(test$SNP,opp)],col="red") 

testS <- test[test$P.x<0.05&test$P.y<0.05,]
cor(testS$BETA.x,testS$BETA.y,method="spearman") #0.92
plot(testS$BETA.x,testS$BETA.y,main="beta for protein and rna\n for associations with p<0.05 in both",xlab="beta rna",ylab="beta protein")
abline(0,1)
abline(v=0,h=0)
testS <- testS[testS$CHR.x!=6,]
cor(testS$BETA.x,testS$BETA.y,method="spearman") #0.92
plot(testS$BETA.x,testS$BETA.y,main="beta for protein and rna\n for associations with p<0.05 in both",xlab="beta rna",ylab="beta protein")
abline(0,1)
abline(v=0,h=0)

############
#looking at p value for rna assoc of ribo QTLs
###########

rb <- read.table("qtlMapping.ribo2.21YRI.bestHits.txt",header=T,as.is=T)
rna <- read.table("qtlMapping.rna1.21YRI.txt",header=T,as.is=T)
dim(rb)
dim(rna)
rb[1:5,]
hist(rb$EMP2,nclass=50)
hist(rna$EMP2,nclass=50)

par(mfrow=c(2,2))
sign <- rb$SNP[rb$EMP2<0.05]
length(sign) #603
hist(rna$P[is.element(rna$SNP,sign)],nclass=50,main="nom.pValue as eQTL \nof SNP with p.emp2<0.05 (n=603, fdr=60%) for rQTL",xlab="nominal p as eQTL ")
qobj <- qvalue(rna$P[is.element(rna$SNP,sign)],fdr.level=0.5)
qplot(qobj) #pie0 is 60%

sign <- rb$SNP[rb$EMP2<0.01]
length(sign) #159

hist(rna$P[is.element(rna$SNP,sign)],nclass=50,main="nom.pValue as eQTL \nof SNP with p.emp2<0.01 (n=159, fdr=45%) for rQTL",xlab="nominal p as eQTL ")
qobj <- qvalue(rna$P[is.element(rna$SNP,sign)],fdr.level=0.5)
qplot(qobj) #pie0 is 44%

sign <- rb$SNP[rb$EMP2<0.005]
length(sign) #98
hist(rna$P[is.element(rna$SNP,sign)],nclass=50,main="nom.pValue as eQTL \nof SNP with p.emp2<0.005 (n=98, fdr=36%) for rQTL",xlab="nominal p as eQTL ")
qobj <- qvalue(rna$P[is.element(rna$SNP,sign)],fdr.level=0.5)
qplot(qobj) #pie0 is 20%


sign <- rb$SNP[rb$EMP2<0.001]
length(sign) #35
hist(rna$P[is.element(rna$SNP,sign)],nclass=50,main="nom.pValue as eQTL \nof SNP with p.emp2<0.001 (n=35; fdr=21%) for rQTL",xlab="nominal p as eQTL ")
qobj <- qvalue(rna$P[is.element(rna$SNP,sign)],fdr.level=0.5)
qplot(qobj) #pie0 is 19%


par(mfrow=c(1,2))
hist(pRNA1,nclass=1000)
hist(pRNA1[bestRNA1$CHR!=6],nclass=1000)

###########plot the hist of emp2
bestR2 <- read.table("qtlMapping.ribo2.21YRI.bestHits.txt",header=T,as.is=T)
pR2 <- bestR2$EMP2
bestRNA1 <- read.table("qtlMapping.rna1.21YRI.bestHits.txt",header=T,as.is=T)
pRNA1 <- bestRNA1$EMP2
best <- read.table("qtlMapping.protein.21YRI.bestHits.txt",header=T,as.is=T)
p <- best$EMP2
best[1:6,]

par(mfrow=c(2,2))
hist(pR2,nclass=50,main= "p.emp2 Riboseq2\n(one per gene ~10000 genes)\n21 unrelated YRI",ylim=c(0,350))
hist(pRNA1,nclass=50,main="p.emp2 RNA1\n(one per gene ~10000 genes)",ylim=c(0,350))
hist(p,nclass=50,,main="p.emp2 protein\n(one per gene ~2000 genes)",ylim=c(0,100))




#are the rQTL, eQTL?
bestR2 <- bestR2[order(bestR2$EMP2,decreasing=F),]
head(bestR2)

top67 <- bestR2[1:67,]
rna <- read.table("qtlMapping.rna1.21YRI.txt",header=T,as.is=T)
head(rna)

top67 <- merge(top67,rna,by="SNP")
dim(top67)
head(top67)
hist(top67$EMP2.y)
plot(top67$BETA.x,top67$BETA.y)
abline(0,1)
abline(v=0,h=0)
#results:  fdr10      fdr30    fdr50
# protein    4         29        48    (logRatio)

# ribo1      0          55        227  ()
# ribo2      0          67        256  (qn)
# ribo3      0          48        156  (sv20)
# ribo4      15         60        164  (sv3) 
 
# rna1       22         128       488  (qn)
# rna2       15         100       513   (sv20)
# rna3       31         164       398   (sv3)