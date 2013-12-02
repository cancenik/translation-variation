library("edgeR")
library("goseq")
library(GO.db)

### Calculates Ribo to RNA differences in Kotaro's data
### NEED TO AUTOMATE TABLE GENERATION FOR A GIVEN SET OF CONTRASTS. 
### WE SHOULD ALSO DEFINE WHAT WE SHOULD USE FOR MIN_COUNTS

## DATA INPUT _Split by Region
#READ COUNTS transcript
data <- read.table ("Reformatted_Transcript_Counts_All_Libraries.tsv", header=T)
X <- split (data, data$REGION)
CDS <- X[[1]]
UTR3 <- X[[2]]
UTR5 <- X[[3]]

summed_cds_counts <- as.matrix(CDS[,c(16,18,4,6,8,10,12,14,15,17,3,5,7,9,11,13)])
rownames(summed_cds_counts) <- CDS$ID
keep <- rowSums(summed_cds_counts) > 50
# isexpr <- rowSums(cpm(summed_cds_counts) > 1) >= 4
summed_cds_counts <- summed_cds_counts[keep,]
cor (summed_cds_counts)

d <- DGEList(counts=summed_cds_counts )
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors (d, method= "TMM")
d$samples
#apply(d$counts, 2, quantile)

pdf ("All_Samples_CDSCounts_MDS.pdf")
plotMDS(d)
dev.off()

s1 <- unlist(strsplit(colnames(summed_cds_counts), "_"))
s2 <- unlist(strsplit(sub(".","_",s1[seq(1, length(s1), 2)], fixed=TRUE), ".", fixed=TRUE))

Stage <- s2[seq(1, length(s2), 2)]
Sample <- s2[seq(2, length(s2), 2)]
Treatment <- s1[seq(2, length(s1), 2)]

experimental_design <- data.frame (Sample=Sample, Treatment=Treatment, Stage=Stage)
Treatment <- factor(experimental_design$Treatment, levels=c("RNA", "RIBO") )
Sample <- factor(experimental_design$Sample, levels=c("1", "2", "3") )
Stage <- factor(experimental_design$Stage, levels=c("E9_5", "E10_5", "E11_5"))
Group <- factor(paste(experimental_design$Treatment, experimental_design$Stage, sep="."))

design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)

my.contrasts <- makeContrasts(
E95.RibovsRNA=RIBO.E9_5-RNA.E9_5,
E105.RibovsRNA=RIBO.E10_5-RNA.E10_5, 
E115.RibovsRNA=RIBO.E11_5-RNA.E11_5,
RNA_105vs95=RNA.E10_5-RNA.E9_5,
RNA_115vs95=RNA.E11_5-RNA.E9_5,
RNA_115vs105=RNA.E11_5-RNA.E10_5,
RIBO_105vs95=RIBO.E10_5-RIBO.E9_5,
RIBO_115vs95=RIBO.E11_5-RIBO.E9_5,
RIBO_115vs105=RIBO.E11_5-RIBO.E10_5,
TE.105vs95=(RIBO.E10_5-RNA.E10_5)-(RIBO.E9_5-RNA.E9_5),
TE.115vs95=(RIBO.E11_5-RNA.E11_5)-(RIBO.E9_5-RNA.E9_5),
TE.115vs105=(RIBO.E11_5-RNA.E11_5)-(RIBO.E10_5-RNA.E10_5),
levels=design)

d <- estimateGLMCommonDisp(d, design, verbose=TRUE)
d <- estimateGLMTrendedDisp(d, design)
fit <- glmFit(d, design)

te115lrt <- glmLRT(d, fit, contrast=my.contrasts[,"TE.115vs95"])
E95_RibovsRNA <- glmLRT(d, fit, contrast=my.contrasts[,"E95.RibovsRNA"])
E105_RibovsRNA <- glmLRT(d, fit, contrast=my.contrasts[,"E105.RibovsRNA"])
E115_RibovsRNA <- glmLRT(d, fit, contrast=my.contrasts[,"E115.RibovsRNA"])

rna105 <- glmLRT(d, fit, contrast=my.contrasts[,"RNA_105vs95"])
rna115 <- glmLRT(d, fit, contrast=my.contrasts[,"RNA_115vs95"])
rna115v105 <- glmLRT(d, fit, contrast=my.contrasts[,"RNA_115vs105"])

topTags(E95_RibovsRNA)
summary(de <- decideTestsDGE(E95_RibovsRNA))
detags <- rownames(d)[as.logical(de)]
detags.over <- rownames(d)[which(de==1)]
detags.under <- rownames(d)[which(de==-1)]
o <- order(E95_RibovsRNA$table$PValue)
cpm(d)[o[1:100],]
write.table(topTags(E95_RibovsRNA, 14352), file="CDS_RiboSeq_vs_RNASeq_95_CountGT50.txt")
write.table(d[o[1:2005],]$counts, file="CDS_RiboSeq_vs_RNASeq_95_RawCounts.txt")

topTags(E105_RibovsRNA)
summary(de <- decideTestsDGE(E105_RibovsRNA))
detags <- rownames(d)[as.logical(de)]
detags.over <- rownames(d)[which(de==1)]
detags.under <- rownames(d)[which(de==-1)]
o <- order(E105_RibovsRNA$table$PValue)
cpm(d)[o[1:100],]
write.table(topTags(E105_RibovsRNA, 14352), file="CDS_RiboSeq_vs_RNASeq_105_CountGT50.txt")
write.table(d[o[1:3232],]$counts, file="CDS_RiboSeq_vs_RNASeq_105_RawCounts.txt")

topTags(E115_RibovsRNA)
summary(de <- decideTestsDGE(E115_RibovsRNA))
detags <- rownames(d)[as.logical(de)]
detags.over <- rownames(d)[which(de==1)]
detags.under <- rownames(d)[which(de==-1)]
o <- order(E115_RibovsRNA$table$PValue)
cpm(d)[o[1:100],]
write.table(topTags(E115_RibovsRNA, 14352), file="CDS_RiboSeq_vs_RNASeq_115_CountGT50.txt")
write.table(d[o[1:4858],]$counts, file="CDS_RiboSeq_vs_RNASeq_115_RawCounts.txt")

topTags(te115lrt)
summary(de <- decideTestsDGE(te115lrt))
detags <- rownames(d)[as.logical(de)]
detags.over <- rownames(d)[which(de==1)]
detags.under <- rownames(d)[which(de==-1)]
o <- order(te115lrt$table$PValue)
cpm(d)[o[1:100],]
write.table(topTags(te115lrt, 14352), file="CDS_TE_115_vs_95_CountGT50.txt")
write.table(d[o[1:100],]$counts, file="CDS_TE_115_vs_95_RawCounts.txt")

topTags(rna105)
summary(de <- decideTestsDGE(rna105))
detags <- rownames(d)[as.logical(de)]
o <- order(rna105$table$PValue)
cpm(d)[o[1:100],]
write.table(topTags(rna105, 14149), file="CDSOnly_Differential_RNA_Expression_105_95_rowCount50.txt")
write.table(d[o[1:615],]$counts, file="CDS_Only_Differential_RNA_Expression_105_95_RawCounts.txt")

## GOSEQ_ANALYSIS

## Need to modify goseq to output data to calculate effect size
## Current modification works only for Wallenius
goseq_modification <- 
function (pwf, genome, id, gene2cat = NULL, test.cats = c("GO:CC", 
    "GO:BP", "GO:MF"), method = "Wallenius", repcnt = 2000) 
{
    if (any(!test.cats %in% c("GO:CC", "GO:BP", "GO:MF", "KEGG"))) {
        stop("Invaled category specified.  Valid categories are GO:CC, GO:BP, GO:MF or KEGG")
    }
    if ((missing(genome) | missing(id))) {
        if (is.null(gene2cat)) {
            stop("You must specify the genome and gene ID format when automatically fetching gene to GO category mappings.")
        }
        genome = "dummy"
        id = "dummy"
    }
    if (!any(method %in% c("Wallenius", "Sampling", "Hypergeometric"))) {
        stop("Invalid calculation method selected.  Valid options are Wallenius, Sampling & Hypergeometric.")
    }
    if (!is.null(gene2cat) && (!is.data.frame(gene2cat) & !is.list(gene2cat))) {
        stop("Was expecting a dataframe or a list mapping categories to genes.  Check gene2cat input and try again.")
    }
    pwf = unfactor(pwf)
    gene2cat = unfactor(gene2cat)
    if (is.null(gene2cat)) {
        message("Fetching GO annotations...")
        gene2cat = getgo(rownames(pwf), genome, id, fetch.cats = test.cats)
        names(gene2cat) = rownames(pwf)
        cat2gene = reversemapping(gene2cat)
        gene2cat = reversemapping(cat2gene)
    }
    else {
        message("Using manually entered categories.")
        if (class(gene2cat) != "list") {
            genecol_sum = as.numeric(apply(gene2cat, 2, function(u) {
                sum(u %in% rownames(pwf))
            }))
            genecol = which(genecol_sum != 0)
            if (length(genecol) > 1) {
                genecol = genecol[order(-genecol_sum)[1]]
                warning(paste("More than one possible gene column found in gene2cat, using the one headed", 
                  colnames(gene2cat)[genecol]))
            }
            if (length(genecol) == 0) {
                genecol = 1
                warning(paste("Gene column could not be identified in gene2cat conclusively, using the one headed", 
                  colnames(gene2cat)[genecol]))
            }
            othercol = 1
            if (genecol == 1) {
                othercol = 2
            }
            gene2cat = split(gene2cat[, othercol], gene2cat[, 
                genecol])
            cat2gene = reversemapping(gene2cat)
            gene2cat = reversemapping(cat2gene)
        }
        if (sum(unique(unlist(gene2cat, use.names = FALSE)) %in% 
            rownames(pwf)) > sum(unique(names(gene2cat)) %in% 
            rownames(pwf))) {
            gene2cat = reversemapping(gene2cat)
        }
        gene2cat = gene2cat[names(gene2cat) %in% rownames(pwf)]
        cat2gene = reversemapping(gene2cat)
        gene2cat = reversemapping(cat2gene)
    }
    nafrac = (sum(is.na(pwf$pwf))/nrow(pwf)) * 100
    if (nafrac > 50) {
        warning(paste("Missing length data for ", round(nafrac), 
            "% of genes.  Accuarcy of GO test will be reduced.", 
            sep = ""))
    }
    pwf$pwf[is.na(pwf$pwf)] = pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)], 
        pwf$bias.data)]
    cats = names(cat2gene)
    DE = rownames(pwf)[pwf$DEgenes == 1]
    num_de = length(DE)
    num_genes = nrow(pwf)
    pvals = data.frame(category = cats, over_represented_pvalue = NA, 
        under_represented_pvalue = NA, num_in_category=NA, num_de = NA, stringsAsFactors = FALSE)
    if (method == "Sampling") {
        num_DE_mask = rep(0, length(cats))
        a = table(unlist(gene2cat[DE], FALSE, FALSE))
        num_DE_mask[match(names(a), cats)] = as.numeric(a)
        num_DE_mask = as.integer(num_DE_mask)
        gene2cat = gene2cat[rownames(pwf)]
        names(gene2cat) = rownames(pwf)
        message("Running the simulation...")
        lookup = matrix(0, nrow = repcnt, ncol = length(cats))
        for (i in 1:repcnt) {
            a = table(as.character(unlist(gene2cat[order(runif(num_genes)^(1/pwf$pwf), 
                decreasing = TRUE)[1:num_de]], FALSE, FALSE)))
            lookup[i, match(names(a), cats)] = a
            pp(repcnt)
        }
        message("Calculating the p-values...")
        pvals[, 2] = (colSums(lookup >= outer(rep(1, repcnt), 
            num_DE_mask)) + 1)/(repcnt + 1)
        pvals[, 3] = (colSums(lookup <= outer(rep(1, repcnt), 
            num_DE_mask)) + 1)/(repcnt + 1)
    }
    if (method == "Wallenius") {
        message("Calculating the p-values...")
        degenesnum = which(pwf$DEgenes == 1)
        cat2genenum = relist(match(unlist(cat2gene), rownames(pwf)), 
            cat2gene)
        alpha = sum(pwf$pwf)
        pvals[, 2:5] = t(sapply(cat2genenum, function(u) {
            num_de_incat = sum(degenesnum %in% u)
            num_incat = length(u)
            avg_weight = mean(pwf$pwf[u])
            weight = (avg_weight * (num_genes - num_incat))/(alpha - 
                num_incat * avg_weight)
            c(dWNCHypergeo(num_de_incat, num_incat, num_genes - 
                num_incat, num_de, weight) + pWNCHypergeo(num_de_incat, 
                num_incat, num_genes - num_incat, num_de, weight, 
                lower.tail = FALSE), pWNCHypergeo(num_de_incat, 
                num_incat, num_genes - num_incat, num_de, weight), num_incat, num_de_incat)
        }))
    }
    if (method == "Hypergeometric") {
        message("Calculating the p-values...")
        degenesnum = which(pwf$DEgenes == 1)
        cat2genenum = relist(match(unlist(cat2gene), rownames(pwf)), 
            cat2gene)
        pvals[, 2:3] = t(sapply(cat2genenum, function(u) {
            num_de_incat = sum(degenesnum %in% u)
            num_incat = length(u)
            c(dhyper(num_de_incat, num_incat, num_genes - num_incat, 
                num_de) + phyper(num_de_incat, num_incat, num_genes - 
                num_incat, num_de, lower.tail = FALSE), phyper(num_de_incat, 
                num_incat, num_genes - num_incat, num_de))
        }))
    }
    pvals = pvals[order(pvals$over_represented_pvalue), ]
    return(pvals)
}

reversemapping <- function (map)
{
    tmp <- unlist(map, use.names = FALSE)
    names(tmp) <- rep(names(map), times = as.numeric(summary(map)[,
        1]))
    return(split(names(tmp), as.vector(tmp)))
}

go_ontology <- read.table("goseq_input_gene_go_associations_enhanced.txt", header =F)
assayed.tags <-	rownames(d)
gene.vector <- as.integer(assayed.tags%in%detags.over)
#gene.vector <- as.integer(assayed.tags%in%detags.under)
names(gene.vector)=assayed.tags
# normalize by average cpm as opposed to length
normalization.factor <- log10(apply (cpm(d), 1, mean))
#normalization.factor <- rowSums(cpm(d))
# Used mm9 (43), gene ID code -> knownGene
pwf <- nullp(gene.vector, "mm9", "knownGene", bias.data=normalization.factor, plot.fit=F)
pdf ("DE_vs_Log10_MeanCPM_CDSOnly.pdf")
plotPWF(pwf, binsize=250)
dev.off()
# Should use my own GO category table
GO.wall=goseq(pwf,"mm9", "knownGene", gene2cat= go_ontology)
# Permutation p_value
# GO.perm <- goseq(pwf,"mm9", "knownGene", gene2cat= go_ontology, method="Sampling", repcnt=1000)

head(GO.wall)
enriched.GO <- GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")< 0.05]	
for(go in enriched.GO[1:10]){
       print(GOTERM[[go]])	
       cat("--------------------------------------\n")
}

# Modified GOSEQ

GO.wall=goseq_modification(pwf,"mm9", "knownGene", gene2cat= go_ontology)
## Calculate odds ratio for significant at 5%
fmat <- matrix(nrow=2, ncol =2)
enriched.GO.table <- GO.wall[p.adjust(GO.wall$over_represented_pvalue,method="BH")< 0.05,]
enriched.GO.table$log10odds <- apply(enriched.GO.table, 1, function (x) 
{
fmat[1,] <-  c(as.numeric(x[5]), length(detags.over) - as.numeric(x[5]))
fmat[2,] <- c(as.numeric(x[4])- as.numeric(x[5]), length(gene.vector) - length(detags.over) - (as.numeric(x[4])- as.numeric(x[5]))  )
log10(fisher.test(fmat)$estimate )
} ) 

for(go in enriched.GO.table$category[enriched.GO.table$log10odds > 0.5]){
       print(GOTERM[[go]])	
       cat("--------------------------------------\n")
}


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

pdf("RNA_Pairs_Log10.pdf")
pairs(log10(summed_counts[,c(1,2,3,4,11,5, 12)]), pch=19, cex=0.25,
diag.panel= panel.hist, lower.panel=panel.smoothScatter)
dev.off()

pdf("Ribo_Pairs_Log10.pdf")
pairs(log10(summed_counts[,c(6,7,8,9,13,10,14)]), pch=19, cex=0.25,
diag.panel= panel.hist, lower.panel=panel.smoothScatter)
dev.off()


