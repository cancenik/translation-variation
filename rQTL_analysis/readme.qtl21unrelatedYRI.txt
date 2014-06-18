#Results file for QTL mapping in the 21 YRI unrelated:

#there are 2 types of result file per phenotype:
1.Plink results for each SNP tested:
-> eg: qtlMapping.protein.21YRI.txt
2.The most significant test per gene was extracted:
The EMP2  p value in this file can be used for determining FDR.
-> eg:qtlMapping.protein.21YRI.bestHits.txt


Key to the phenotype file used:   model : pheno ~ geno
 log2ratio (sup Table 1)                                          -> protein
TMM_VarianceMeanDetrended_CPM_GT1_in36_RiboProfiling_Expression    -> ribo1   
TMM_VarianceMeanDetrended_QN_FullModel_RiboProfiling_Expression    -> ribo2
Top20_SVA_Removed_QN_FullModel_RiboProfiling_Expression            -> ribo3
Top3_SVA_Removed_QN_FullModel_RiboProfiling_Expression             -> ribo4

TMM_VarianceMeanDetrended_QN_FullModel_RNA_Expression              -> rna1
Top20_SVA_Removed_QN_FullModel_RNA_Expression                      -> rna2
Top3_SVA_Removed_QN_FullModel_RNA_Expression                       -> rna3


Other tests:
qtlMapping.ribo2.CovarRNA1   ->   ribo2 ~ geno + rna1



Headers in the results files are:
SNP	ENST	CHR	BP	A1	TEST	NMISS	BETA	STAT	P	EMP1	EMP2	ALT	REF	ALT.freq
ENST00000000233.5;646;A;C	ENST00000000233.5	7	646	C	ADD	21	-0.03966	-1.454	0.1622	0.1615	0.1615	C	A	0.2381

SNP= unique SNP id =ENST;bp in transcript;ref allele;alt allele
BP: is the bp in the transcript
A1 is the alternate allele
NMISS: number of individuals with non missing geno and pheno: this is mainly 21 for all (except those on chr X)
       this is not 22 (we have 22 unrelated YRI in the geno file, but only 21 for the phenotypes ) 
BETA: increase in phenotype for each increasing number of alt allele (geno is coded as 0,1,2 copies of the alternate allele)
P: p value from a linear model pheno~geno
EMP1: p value determined by permutations (they are pretty similar, except EMP1 plateaus at 10-5 because I only ran 10000 permutations)
EMP2: p value correcting for the number of tests within a gene, based on the maxT method, and 10,000 permutations

Using qvalue to get the number of QTL at FDR 10% (based on the distribution of EMP2, 1 per gene), I get:




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




