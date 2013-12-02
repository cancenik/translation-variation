library(scatterplot3d)
library(rgl)
library (psych)
library(hash)

# NOTE THAT SOME SECTIONS OF CODE ARE COMPUTATIONALLY INTENSIVE

emi <- read.table('~/Desktop/Emi_Ingolia/Transcript_Counts_EMI')
ingolia <- read.table('~/Desktop/Emi_Ingolia/Transcript_Counts_Ingolia')
 
#### Test for the effect of deeper transcriptome coverage. Randomly sample from EMI
set.seed(403)
s1 <- sort(sample(1:sum(emi$V2), sum(ingolia$V2)))

emi_sampled <- emi
emi_sampled[,2] <- 0
for (i in 1:nrow(emi)) {
	emi_sampled$V2[i] <- sum(c(1+sum(emi$V2[1:i-1]):sum(emi$V2[1:i])) %in% s1)
}
sum (emi_sampled$V2)
sum(ingolia$V2)
###### END OF COVERAGE TEST -- No Difference in results so I used the non-subsampled data

#Remove those with less than 10 reads 
emi10 <- emi[emi$V2 > 10, ]
emi10 <- emi_sampled[emi_sampled$V2 > 10, ]
ingolia10 <- ingolia[ingolia$V2 > 10, ]

#emi_ingolia_merge <- merge(emi, ingolia, by="V1")
emi_ingolia_merge <- merge(emi10, ingolia10, by="V1")
colnames(emi_ingolia_merge) <- c("ID", "EMI", "Ingolia")

plot(log10(emi_ingolia_merge$EMI), log10(emi_ingolia_merge$Ingolia), pch=20, col="blue", cex=0.2, main="HEK293 cell Ribosome Profiling\nlog10 Read Numbers Per Transcript", xlab= "RNase A + RNase T1", ylab = "RNase I (Ingolia)")
legend("topleft", "Pearson Correlation: 0.88", bty= "n")

prot <- read.table('~/Desktop/Emi_Ingolia/Geiger_Mann_HEK293_Data_Quantification_in_2_SplitbyGene.csv')
colnames(prot)[1] <- c("ID")
all<- merge (emi_ingolia_merge, prot, by="ID")

for (i in 4:9) {
	all[,i][is.nan(all[,i])] <- NA
}

## ALL
scatterplot3d(log10(all$EMI), log10(all$Ingolia), prot_ibaq_means , pch=20, cex.symbols = 0.2, color = "blue", box = F, xlab= "RNase A + RNase T1", ylab = "RNase I", zlab = "Quantitative Mass Spectrometry", angle=45, ylim= c(1,5), xlim = c(1,5))
legend("topleft", "Pearson Correlation: 0.66 vs 0.62\n p-value < 0.0001", bty= "n", cex = 1.5 )

plot3d(log10(all$EMI), prot_ibaq_means, log10(all$Ingolia), xlab = "RNase A + T1", ylab = "Protein Quantification", zlab= "RNase I")

prot_ibaq_means <- apply(all[,4:6],1, function(x) {mean(x, na.rm=T)} )
plot(log10(all$EMI),prot_ibaq_means, pch=20, col="blue", cex=0.2 )
c_i_e <- cor.test(log10(emi_ingolia_merge$EMI), log10(emi_ingolia_merge$Ingolia))
c_e_p <- cor.test(log10(all$EMI),prot_ibaq_means)
# 0.659 -- emi10
c_i_p <- cor.test(log10(all$Ingolia),prot_ibaq_means)
# 0.6249- ing10

# Permutation Test for Significance of Correlation Coefficient Differences
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
	perm <- apply (all[,2:3], 1, swap)
	cor_dif <- c(cor_dif,cor.test(log10(perm[1,]), prot_ibaq_means)$estimate - cor.test(log10(perm[2,]), prot_ibaq_means)$estimate)
}
quantile(cor_dif)
length(which (abs(cor_dif) > (c_e_p$estimate - c_i_p$estimate))) 
hist (abs(cor_dif), xlim = c(0, 0.035), xlab= "Absolute Difference in Correlation Coefficient", main= "", breaks = 50)
abline(v= c_e_p$estimate - c_i_p$estimate, col = "red")

#### TRIPLICITY MAXIMIZATION FOR RNase A+ T1
# Use the last two nucleotides and length to reassign to frame. 
# We will maximize the triplicity
emi_modulo_counts <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/ALL_NUCLEOTIDES/appris_aln_Q2.sam_modulo_Counts')
# Simplest approach -- Take max value and assign to 0 
# Use only 24nt to 31nt for simplicity
mod_counts_emi_24_31 <- emi_modulo_counts[284:692,]
# REMOVE Positions with N or NA
mod_counts_emi_24_31 <- mod_counts_emi_24_31[-grep ("N", mod_counts_emi_24_31$V2), ]
mod_counts_emi_24_31 <- mod_counts_emi_24_31[!is.na(mod_counts_emi_24_31$V2),]
mod_counts_emi_24_31$stratify <- paste(mod_counts_emi_24_31$V1, mod_counts_emi_24_31$V2, sep="_")
aggregate(V4 ~ V3, data=mod_counts_emi_24_31, FUN=sum )

mod_max <- aggregate (V4 ~ stratify, data= mod_counts_emi_24_31, FUN=max )
sum (mod_max$V4) / sum(mod_counts_emi_24_31$V4)
m1<- merge (mod_max , mod_counts_emi_24_31)[,c("stratify", "V3")]
h1 <- hash(m1$stratify, m1$V3)
for (i in 1:nrow(mod_counts_emi_24_31)) {
	if (h1[[mod_counts_emi_24_31$stratify[i]]] == 0 ) {
		next;
	}
	else if (h1[[mod_counts_emi_24_31$stratify[i]]] == 1 ) {
		mod_counts_emi_24_31$V3[i] <- (mod_counts_emi_24_31$V3[i] +2) %% 3
	}
	else {
		mod_counts_emi_24_31$V3[i] <- (mod_counts_emi_24_31$V3[i] +1) %% 3

	}
}

a1 <- aggregate(V4 ~ V3, data=mod_counts_emi_24_31, FUN=sum )
#  V3      V4
# 1  0 1483548
# 2  1  904436
# 3  2  807502
barplot (cbind(a1$V4/sum(mod_counts_emi_24_31$V4)), names= c("RNase A + T1"), beside= T, ylim = c(0,0.5)) 

b1 <- c()
n1 <- c() 
c1 <- aggregate(V4 ~ V1+V3, data=mod_counts_emi_24_31, FUN=sum )
for (i in 24:31) {
	b1 <- cbind (b1, c1[c1$V1 ==i, 3] )
	n1 <- c(n1, i)
}
barplot (b1, names= n1, beside= T, xlab = "RNase A + T1", ylab= "Number of Reads") 

a1$V4 <- a1$V4/sum(mod_counts_emi_24_31$V4)

### TRIPLICITY OF MAPPED READS USING HEATMAPS
library(gplots)
generate_heatmatrix <- function (x) {
	tmp <- rep (0, 71) 
    ht_mat_ret <- matrix (nrow = 8, ncol= 71 )
	for (i in 24:31) {
	 for (j in -35: 35) { 
 	  if (j %in% x[x$V1 == i, 2]){
  	   tmp[j+36] <- x[x$V1 == i & x$V2 == j, 3 ]  	
	  }	
	 } 	
  	 ht_mat_ret[i-23, ]	<- tmp
 	 tmp <- rep(0,71)
	}
ht_mat2 <- ht_mat_ret[,c(-1,-71)]
rownames(ht_mat2) <- c("24","25","26","27","28", "29", "30","31")
colnames (ht_mat2) <- as.character(seq(-34, 34, 1))
return (ht_mat2)
}

# This function sums a vector to three values using mod 3
calculate_mod_sum <- function (a) { 
	mod0 <- sum(a[seq(1, length(a), 3)])
	mod1 <- sum(a[seq(2, length(a), 3)])
	mod2 <- sum(a[seq(3, length(a), 3)])
	return(c(mod0, mod1, mod2))
}

ing_readspileup <- read.table('~/Desktop/Emi_Ingolia/INGOLIA_HEK_SAM_COUNTS/appris_aln_Q2.sam_readpileup_tocds')
emi_readspileupall <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/ALL_NUCLEOTIDES/appris_aln_Q2.sam_readpileup_tocds')
emi_readspileupG <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_G/appris_aln_Q2.sam_readpileup_tocds')
emi_readspileupC <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_C/appris_aln_Q2.sam_readpileup_tocds')
emi_readspileupA <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_A/appris_aln_Q2.sam_readpileup_tocds')
emi_readspileupT <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_T/appris_aln_Q2.sam_readpileup_tocds')
 
ing_mat <- generate_heatmatrix(ing_readspileup) 
emi_mat <- generate_heatmatrix(emi_readspileupall)
emi_matG <- generate_heatmatrix(emi_readspileupG)
emi_matC <- generate_heatmatrix(emi_readspileupC)
emi_matA <- generate_heatmatrix(emi_readspileupA)
emi_matT <- generate_heatmatrix(emi_readspileupT)

heatmap(log10(ing_mat+1), Rowv= NA, Colv = NA)
ht1 <-  heatmap.2(log10(ing_mat+1), Rowv= NA, Colv = NA, dendrogram = "none", trace= "none")
ht1_breaks  <- ht1$breaks
ht2 <- heatmap.2(log10(emi_mat+1), Rowv= NA, Colv = NA, dendrogram = "none", trace= "none", breaks = ht1_breaks)
htG <- heatmap.2(log10(emi_matG+1), Rowv= NA, Colv = NA, dendrogram = "none", trace= "none", breaks = ht1_breaks)
htC <- heatmap.2(log10(emi_matC+1), Rowv= NA, Colv = NA, dendrogram = "none", trace= "none")
#htA <- heatmap.2(log10(emi_matA+1), Rowv= NA, Colv = NA, dendrogram = "none", trace= "none")
htA <- heatmap.2(emi_matA, Rowv= NA, Colv = NA, dendrogram = "none", trace= "none")
htT <- heatmap.2(log10(emi_matT+1), Rowv= NA, Colv = NA, dendrogram = "none", trace= "none")
htT <- heatmap.2(emi_matT, Rowv= NA, Colv = NA, dendrogram = "none", trace= "none")
htT_breaks <- htT$breaks[c(-8:-15)]
htT <- heatmap.2(emi_matT, Rowv= NA, Colv = NA, dendrogram = "none", trace= "none", breaks= htT_breaks)

ing_mods <- apply(ing_mat, 1, calculate_mod_sum )
emi_mods <- apply(emi_mat, 1, calculate_mod_sum)
p <- c(1/3, 1/3, 1/3)
apply (ing_mods, 2, chisq.test, p = p)
apply (emi_mods, 2, chisq.test, p = p)


#### BARPLOT OF ALIGNED READ LENGHTS AND FRACTION OF LAST NUCLEOTIDE
emi_A <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_A/appris_aln_Q2.sam_lengths_of_aligned')
emi_T <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_T/appris_aln_Q2.sam_lengths_of_aligned')
emi_C <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_C/appris_aln_Q2.sam_lengths_of_aligned')
emi_G <- read.table('~/Desktop/Emi_Ingolia/EMI_HEK_SAMCOUNTS/LAST_G/appris_aln_Q2.sam_lengths_of_aligned')

emi_last_nuc <- c()
for (i in 7:14) {
	emi_last_nuc <- cbind (emi_last_nuc, c(emi_A[i, 2], emi_T[i, 2], emi_C[i, 2] ,emi_G[i, 2]))
}
barplot (emi_last_nuc, names = c("24", "25", "26", "27", "28", "29", "30", "31"), beside = F, ylab = "Last_Nucleotide_Distribution")

ing_A <- read.table('~/Desktop/Emi_Ingolia/INGOLIA_HEK_SAM_COUNTS/LAST_A/appris_aln_Q2.sam_lengths_of_aligned')
ing_T <- read.table('~/Desktop/Emi_Ingolia/INGOLIA_HEK_SAM_COUNTS/LAST_T/appris_aln_Q2.sam_lengths_of_aligned')
ing_C <- read.table('~/Desktop/Emi_Ingolia/INGOLIA_HEK_SAM_COUNTS/LAST_C/appris_aln_Q2.sam_lengths_of_aligned')
ing_G <- read.table('~/Desktop/Emi_Ingolia/INGOLIA_HEK_SAM_COUNTS/LAST_G/appris_aln_Q2.sam_lengths_of_aligned')

ing_last_nuc <- c()
for (i in 7:14) {
	ing_last_nuc <- cbind (ing_last_nuc, c(ing_A[i, 2], ing_T[i, 2], ing_C[i, 2] ,ing_G[i, 2]))
}

barplot(ing_last_nuc, names = c("24", "25", "26", "27", "28", "29", "30", "31"),beside = F, ylab = "Last_Nucleotide_Distribution")
###################

####### RRNA CONTENT
ing_total <- 44583859
emi_total <- 14827585

emi_rRNA <- 9341991 + 219728
emi_trans <- 3355633 +722762 +28273 + 84323
emi_genome <- 408121 + 67930 - 275829 - 9616
emi_ncbi <- 275829 + 9616
emi_rest <- emi_total - (emi_rRNA + emi_trans + emi_genome + emi_ncbi)

ing_rRNA <- 33147308 + 47256
ing_trans <- 2414900 + 675553 + 27329 + 68844
ing_genome <- 6329732 + 492114 - 5866700 - 52216
ing_ncbi <- 5866700 + 52216
ing_rest <- ing_total - (ing_rRNA + ing_trans + ing_genome + ing_ncbi)

# CHANGE BESIDE = F for stacked
#barplot(cbind(c(ing_rRNA/ing_total, ing_ncbi/ing_total, ing_trans/ing_total, ing_genome/ing_total), c(emi_rRNA/emi_total, emi_ncbi/emi_total, emi_trans/emi_total, emi_genome/emi_total)),  names = c("RNase I", "RNase A + T1"), beside = T, ylab = "Fraction of Reads", ylim = c(0,1) )

barplot(cbind(c(ing_rRNA/ing_total, ing_ncbi/ing_total, ing_trans/ing_total, ing_genome/ing_total, ing_rest/ ing_total), c(emi_rRNA/emi_total, emi_ncbi/emi_total, emi_trans/emi_total, emi_genome/emi_total,  emi_rest/ emi_total)),  names = c("RNase I", "RNase A + T1"), beside = F, ylab = "Fraction of Reads", ylim = c(0,1) , col= c("dodgerblue4", "deepskyblue", "firebrick2", "deeppink", "black"))
#########################

##INGOLIA CLIPPED
clipped <- "
1	17.754365618452628	11.957469405907112	45.93916218942578	24.349002786214484
2	26.015237039631906	8.892913077745005	31.382103415941504	33.70974646668159
3	33.312760117961076	7.21420279029682	31.82431561162079	27.648721480121313
4	38.29815404718555	7.142039902826716	21.65150620990435	32.90829984008338
5	44.958035148998654	7.246165030263531	7.192127536559811	40.603672284178
6	46.92552014720603	7.716058576286418	15.968279447419324	29.39014182908823
7	47.61967458179147	10.997183079889625	18.88324276717648	22.49989957114242
8	54.22628624408667	12.612649793280568	13.208010997881543	19.953052964751212
9	43.51047763720947	20.300530736919832	10.264880839498439	25.92411078637226
10	26.34181352493511	19.723582025503894	16.067287939341455	37.86731651021954
11	13.00776812523115	24.950632918518785	23.21093604750544	38.830662908744614
12	14.25462295670727	29.147158392009093	26.448246662542154	30.149971988741488
13	16.299349502249232	37.95191887718827	17.813558938449003	27.935172682113496
14	28.272812813265	25.57067794423089	14.451945489958598	31.70456375254551
15	37.21863331749726	14.635307813978148	17.0328974887526	31.11316137977199
16	20.563513804401722	21.895581089111197	18.64362616075921	38.89727894572787
17	19.094414864357077	21.153810395820607	21.485973208375704	38.26580153144661
18	22.426293336339505	23.874806350881382	26.885656084638164	26.813244228140952
19	29.015802337125162	16.297955378723813	27.624031781597097	27.062210502553924
20	19.20369820425104	10.775187120156653	31.79796638011791	38.223148295474395
21	25.37307270544692	12.41334142894606	25.575019516045955	36.63856634956107
22	35.933441675985	17.920315944365413	26.748650951355664	19.397591428293925
23	42.32551647955466	27.14763884297567	9.479113396220132	21.04773128124954
24	38.737618102150286	30.98656262172496	12.589936342383416	17.685882933741336
25	42.56035215732822	18.507584027738993	15.483846054577008	23.44821776035578
26	35.57127581801902	18.584127335615662	15.844064395938595	30.000532450426725
27	31.09274173956133	28.871980732872405	14.622845741212167	25.412431786354105
28	33.83042049911432	22.698944898527778	22.54769666236978	20.922937939988124
29	43.338929081126246	15.938421336575948	12.039302670953248	28.68334691134456
30	46.33977059245424	16.193251161516827	13.050771040579265	24.416207205449673
31	23.42167700921601	13.67343118034091	18.147477386832662	44.75741442361042
"

ing_clipped <- read.table(textConnection(clipped))
plot (ing_clipped[,2], type= 'l', ylim = c(0, 60), xlab = "Position", ylab = "Percentage", main = "RNase I")
lines (ing_clipped[,3], col= 'green')
lines (ing_clipped[,4], col = 'red')
lines (ing_clipped[,5], col = 'blue')
grid (NULL, NA, lwd = 2)

emi_clipped <- "
1	14.221695576184523	38.9446427047965	20.10758326457073	26.72607845444825
2	25.297241593961523	14.89017260733963	22.031079235087844	37.781506563611
3	37.60329817701264	15.916030830374602	16.211662249786464	30.26900874282629
4	28.587676280392255	17.837921684481998	20.595660048483957	32.97874198664179
5	26.94166986734522	13.21482898260236	18.93542340172051	40.90807774833191
6	29.58057566353523	12.691594753967014	24.90232900367794	32.82550057881981
7	29.10584562489441	11.563224894681095	17.8306447071455	41.50028477327899
8	25.640062086981796	13.830869962977788	16.09325456573002	44.43581338431039
9	32.43611147735791	15.570424988290407	23.815462868700465	28.178000665651215
10	25.542460218572344	16.342539934857903	17.009863710105186	41.10513613646457
11	29.37315820479195	13.749137165627442	20.8841898394108	35.993514790169804
12	31.3460283653744	13.364711785499797	18.307917304132804	36.981342544993
13	30.959046938527074	17.886297734931212	19.571433918605084	31.583221407936623
14	33.02797003206996	19.524232673374165	20.557705152774723	26.890092141781157
15	36.10090921751586	13.36062480842295	15.669800577774465	34.868665396286715
16	40.738508664762335	15.88199292062733	16.687552288521697	26.691946126088638
17	31.78225584274175	13.871294617430957	14.27941906925504	40.06703047057225
18	30.74208645575122	18.55842337103446	12.899410119719429	37.80008005349489
19	26.985395019426633	16.386811822324347	16.630445235118763	39.99734792313026
20	31.985383851284233	16.038611499658753	15.01954225709015	36.95646239196686
21	25.059107280086206	20.627613695588625	13.20557239492747	41.1077066293977
22	31.816288121897536	20.420678999645215	12.866461792613165	34.89657108584409
23	26.89426889077487	18.77097456980393	12.92213814685845	41.412618392562756
24	28.404298261585083	20.994842763567384	11.087289119981257	39.513569854866276
25	26.55045151851664	17.390303188189506	13.832627226224353	42.2266180670695
26	24.83109771695338	15.502598004732365	17.35864456097717	42.30765971733709
27	28.863710569095748	16.639327371130317	15.674420053922866	38.822542005851076
28	28.175366977924355	23.74826385129268	17.727831612010803	30.348537558772165
29	29.88936834605804	8.971395944600802	15.4003607718699	45.738874937471266
30	32.33357803365357	6.698751175999906	11.555324469428264	49.41234632091825
31	23.561944506617184	4.873540645328982	7.0460567641134135	64.51845808394042
"
emi_clip <- read.table(textConnection(emi_clipped))

plot (emi_clip[,2], type= 'l', ylim = c(0, 65), xlab = "Position", ylab = "Percentage", main = "RNase A + T1")
lines (emi_clip[,3], col= 'green')
lines (emi_clip[,4], col = 'red')
lines (emi_clip[,5], col = 'blue')
grid (NULL, NA, lwd = 2)


########## UNUSED #################

# Ingolia Triplicity
ing_modulo_counts <- read.table('~/Desktop/Emi_Ingolia/INGOLIA_HEK_SAM_COUNTS/appris_aln_Q2.sam_modulo_Counts')
mod_counts_ing_24_31 <- ing_modulo_counts[289:726,]
mod_counts_ing_24_31 <- mod_counts_ing_24_31[-grep ("N", mod_counts_ing_24_31$V2), ]
mod_counts_ing_24_31 <- mod_counts_ing_24_31[!is.na(mod_counts_ing_24_31$V2),]
mod_max <- aggregate (V4 ~ V1 + V2, data= mod_counts_ing_24_31, FUN=max )
sum (mod_max$V4) / sum(mod_counts_ing_24_31$V4)
a2 <- aggregate(V4 ~ V3, data=mod_counts_ing_24_31, FUN=sum )
#  V3     V4
# 1  0 832319
# 2  1 328513
# 3  2 831244
a2 <- aggregate(V4 ~ V1+V3, data=mod_counts_ing_24_31, FUN=sum )
a2$V4 <- a2$V4/sum(mod_counts_emi_24_31$V4)
barplot (cbind(a1$V4/sum(mod_counts_emi_24_31$V4), a2$V4/sum(mod_counts_ing_24_31$V4)), names= c("RNase A + T1","RNase I"), beside= T, ylim = c(0,0.5)) 




## ALL _UNUSED PROTEOMICS STATS
prot_lfq_means <- apply(all[,7:9],1, function(x) {mean(x, na.rm=T)} )
prot_ibaq_medians <- apply(all[,4:6],1, function(x) {median(x, na.rm=T)} )
prot_lfq_medians <- apply(all[,7:9],1, function(x) {median(x, na.rm=T)} )
plot(log10(all$Ingolia),prot_ibaq_means, pch=20, col="blue", cex=0.2 )
plot(log10(all$EMI),prot_lfq_means, pch=20, col="blue", cex=0.2 )
plot(log10(all$Ingolia),prot_lfq_means, pch=20, col="blue", cex=0.2 )



#### QIAN -PNAS TIS Mapping HEK DATA
qian <- read.table ('~/Desktop/Emi_Ingolia/Qian_counts')
all_qian <- merge(emi_qian_merge, prot, by="ID")
set.seed(403)
s1 <- sort(sample(1:sum(qian$V2), sum(emi$V2)))
qian_sampled <- qian
qian_sampled[,2] <- 0
system.time(
for (i in 1:100) {
	qian_sampled$V2[i] <- sum(c(1+sum(qian$V2[1:i-1]):sum(qian$V2[1:i])) %in% s1)
}
)
sum (qian_sampled$V2)
sum(emi$V2)

qian10 <- qian[qian$V2 > 10, ]
emi_qian_merge <- merge(emi10, qian10, by="V1")
colnames(emi_qian_merge) <- c("ID", "EMI", "Qian")
cor.test(log10(emi_qian_merge$EMI), log10(emi_qian_merge$Qian))



# SubSampling
#m1 <- cbind(log10(all$EMI),log10(all$Ingolia),prot_ibaq_means )
#for (i in 1:1000) {
#	m2 <- m1[sample(nrow(m1), size = 2000, replace= TRUE), ]
	# SUBSAMPLED FINISH IMPLEMENTATION
#}

cor.test(log10(all$EMI),prot_lfq_means)
cor.test(log10(all$Ingolia),prot_lfq_means)
cor.test(log10(all$EMI),prot_ibaq_means, method="spearman")
cor.test(log10(all$Ingolia),prot_ibaq_means, method="spearman")
cor.test(log10(all$EMI),prot_lfq_means, method="spearman")
cor.test(log10(all$Ingolia),prot_lfq_means, method="spearman")

##COMPLETE PROTEOMICS DATA
complete <- all[complete.cases(all),]

prot_ibaq_medians <- apply(complete[,4:6],1, function(x) {median(x, na.rm=T)} )
prot_lfq_medians <- apply(complete[,7:9],1, function(x) {median(x, na.rm=T)} )
prot_ibaq_means <- apply(complete[,4:6],1, function(x) {mean(x, na.rm=T)} )
prot_lfq_means <- apply(complete[,7:9],1, function(x) {mean(x, na.rm=T)} )
plot(log10(complete$EMI),prot_ibaq_means, pch=20, col="blue", cex=0.2 )
plot(log10(complete$Ingolia),prot_ibaq_means, pch=20, col="blue", cex=0.2 )
plot(log10(complete$EMI),prot_lfq_means, pch=20, col="blue", cex=0.2 )
plot(log10(complete$Ingolia),prot_lfq_means, pch=20, col="blue", cex=0.2 )

cor.test(log10(complete$EMI),prot_ibaq_means)
cor.test(log10(complete$Ingolia),prot_ibaq_means)
cor.test(log10(complete$EMI),prot_lfq_means)
cor.test(log10(complete$Ingolia),prot_lfq_means)

cor.test(log10(complete$EMI),prot_ibaq_means, method="spearman")
cor.test(log10(complete$Ingolia),prot_ibaq_means, method="spearman")
cor.test(log10(complete$EMI),prot_lfq_means, method="spearman")
cor.test(log10(complete$Ingolia),prot_lfq_means, method="spearman")

## INGOLIA ALIGNED READS NUCLEOTIDE DISTR
aligned_reads <- "
1	12.033076507489703	28.667682958303363	42.96435982229808	16.334880711908852
2	24.92380618219094	19.34428597964954	26.47038588814816	29.261521950011353
3	31.828720425447017	24.53922519659037	21.099667705516914	22.5323866724457
4	30.140072788315344	23.569172978885884	23.41678534326777	22.873968889531003
5	30.95464145906861	19.906021891533026	16.551773957193866	32.5875626922045
6	29.35549119029109	20.97474389581083	24.45562810025387	25.21413681364421
7	30.82358121263407	25.3283385963826	21.473584996525698	22.37449519445763
8	26.375100617118324	29.62787145776145	21.17981741002938	22.817210515090846
9	27.30834588897374	25.166319236617063	24.3875255412685	23.137809333140698
10	27.867329879672248	28.076475889730517	20.996126670932146	23.06006755966509
11	25.096833226696386	26.903469484633963	24.396469285119675	23.60322800354998
12	27.930624066926725	24.25784125542644	22.761828101243182	25.04970657640365
13	27.41498283489161	26.003247266998276	23.164640564694228	23.41712933341589
14	28.308325249564852	22.429189628009052	27.995294214773693	21.267190907652406
15	31.530824957173227	22.760796130798813	25.726679187908058	19.9816997241199
16	30.980440720177775	26.520608449773995	25.718767414501247	16.78018341554698
17	27.164558007058677	30.613747222279553	23.82647760968126	18.395217160980508
18	25.19727834994806	33.7457775209318	20.653168493254352	20.403775635865788
19	28.926109809033797	30.53795054571678	19.27032049592966	21.26561914931976
20	26.553460354415463	27.755064831242414	19.976718255832708	25.714756558509404
21	28.67297822830002	27.00477086056124	20.15130747080637	24.170943440332366
22	29.550624490321475	26.62453324177003	19.42036138890081	24.404480879007686
23	27.320877335499596	29.23007491650871	20.44191714053615	23.007130607455547
24	29.73091339082704	26.429100623720352	20.87948584129255	22.960500144160058
25	31.262318264994416	27.797429263929153	19.873281938407747	21.066970532668687
26	28.97779234239727	30.266507133276434	20.29592427752713	20.459776246799173
27	31.2100150225338	27.082223335002503	21.97295943915874	19.73480220330496
28	33.0716731288121	26.701131442305513	19.208217128086524	21.018978300795865
29	29.086157488566457	28.411465240088233	22.722832801816427	19.77954446952888
30	27.533442167588507	24.39999025364879	29.695913842255305	18.370653736507396
31	22.14335009366065	25.732031943212068	32.65552597850734	19.469091984619936
"
ing_aligned <- read.table(textConnection(aligned_reads))
plot (ing_aligned[,2], type= 'l', ylim = c(0, 60), xlab = "Position", ylab = "Percentage", main = "RNase I")
lines (ing_aligned[,3], col= 'green')
lines (ing_aligned[,4], col = 'red')
lines (ing_aligned[,5], col = 'blue')
grid (NULL, NA, lwd = 2)

