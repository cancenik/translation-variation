number_NAs <- 1
linfeng_protein_na <- linfeng_protein_common[apply(is.na(linfeng_protein_common), 1, sum) < number_NAs, ]
linfeng_protein_ribo_rna <- merge (v3$E, linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_nonzero_variance <- 
  merge (v3$E[sig_ribo_rna_random,], linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_zero_variance <- 
  merge (v3$E[-sig_ribo_rna_random,], linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_strict_variance <- 
  merge (v3$E[sig_ribo_rna_random_strict,], linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_rnaonly_variance <- 
  merge (v3$E[sig_rna_only_random,], linfeng_protein_na, by.x="row.names", by.y="HGNC")
linfeng_protein_ribo_rna_rnaall_variance <- 
  merge (v3$E[rna_sig_random,], linfeng_protein_na, by.x="row.names", by.y="HGNC")

plot_density <- function (dataset) { 
  row.names(dataset) <- dataset[,1]
  dataset <- dataset[,-c(1,135)]
  type_prot <- c(type, rep("Prot", dim(dataset)[2] - length(type)))

  sample_labels_joint_prot <- c(sample_labels_joint, colnames(dataset)[134:161])
  rna_replicate_mean_prot <- apply (dataset[,type_prot=="RNA"], 1, function(x) {
    aggregate(x, by= list(as.factor(sample_labels_joint_prot[type_prot=="RNA"])), mean)  
  } )
  ribo_replicate_mean_prot <- apply (dataset[,type_prot=="Ribo"], 1, function(x) {
    aggregate(x, by= list(as.factor(sample_labels_joint_prot[type_prot=="Ribo"])), mean)  
  } )
  rna_samples <- match(as.character(rna_replicate_mean_prot[[1]]$Group.1), sample_labels_joint_prot[type_prot=="Prot"],)
  rna_samples <- rna_samples[!is.na(rna_samples)]
  rna_in_prot <- as.character(rna_replicate_mean_prot[[1]]$Group.1) %in% sample_labels_joint_prot[type_prot=="Prot"]
  across_ind_rna_correlation <- c()
  across_ind_rna_correlation_pval <- c()
  
  for (i in 1:length(rna_replicate_mean_prot)) { 
    #  cor1 <-  cor.test(rna_replicate_mean_prot[[i]]$x[rna_in_prot], as.numeric(as.matrix(linfeng_protein_ribo_rna[i,type_prot=="Prot"]))[rna_samples], use="pairwise.complete.obs")
    cor1 <-  cor.test(rna_replicate_mean_prot[[i]]$x[rna_in_prot], as.numeric(as.matrix(dataset[i,type_prot=="Prot"]))[rna_samples], use="pairwise.complete.obs", method="spearman")
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
    cor1 <-cor.test(ribo_replicate_mean_prot[[i]]$x[ribo_in_prot], as.numeric(as.matrix(dataset[i,type_prot=="Prot"]))[ribo_samples], use="pairwise.complete.obs", method="spearman")  
    across_ind_ribo_correlation <- c(across_ind_ribo_correlation, cor1$estimate)
    across_ind_ribo_correlation_pval <- c(across_ind_ribo_correlation_pval, cor1$p.value)
  } 
  # Remove one column which is not shared
  ribo_replicate_mean_rna <- lapply(ribo_replicate_mean_prot, function(x){x <- x[-24,]})
  c2 <- mapply (cbind, rna_replicate_mean_prot, ribo_replicate_mean_rna, SIMPLIFY=F)
  across_ind_rna_ribo <- as.numeric(lapply(c2, function(x){ cor(x[,2], x[,4],method="spearman") }))
  # Histograms of Across Ind Ribo-Prot, RNA-Prot and Ribo-RNA correlations
  d1 <- density(across_ind_rna_correlation, bw = 0.06)
  plot(d1, xlim=c(-1,1), ylim = c(0,2.5))
}
pdf ('~/Desktop/Presentation_CorrelationDistribution.pdf', height=12, width=4 )
par (mfrow = c(3,1))
plot_density (linfeng_protein_ribo_rna_zero_variance)
plot_density (linfeng_protein_ribo_rna_rnaonly_variance)
plot_density (linfeng_protein_ribo_rna_strict_variance)
dev.off()

#plot_density (linfeng_protein_ribo_rna_nonzero_variance)
# plot_density (linfeng_protein_ribo_rna_rnaall_variance)
# plot_density (linfeng_protein_ribo_rna_riboall_variance)

i = 4
## From rnaonly
plot(rna_replicate_mean_prot[[i]]$x[rna_in_prot]-4.2, as.numeric(as.matrix(dataset[i,type_prot=="Prot"]))[rna_samples])
barplot (rna_replicate_mean_prot[[i]]$x[rna_in_prot] - 4.2 , col = "Green2" )
barplot ( as.numeric(as.matrix(dataset[i,type_prot=="Prot"]))[rna_samples] , col = "Blue2" )
row.names(dataset)[i]
# AACS


i = 456

## From rnaonly
plot(rna_replicate_mean_prot[[i]]$x[rna_in_prot] - mean(rna_replicate_mean_prot[[i]]$x[rna_in_prot])
     , as.numeric(as.matrix(dataset[i,type_prot=="Prot"]))[rna_samples])
barplot (rna_replicate_mean_prot[[i]]$x[rna_in_prot] - mean(rna_replicate_mean_prot[[i]]$x[rna_in_prot]) , col = "Green2" )
barplot ( as.numeric(as.matrix(dataset[i,type_prot=="Prot"]))[rna_samples] , col = "Blue2" )




# Venn Diagram Colors 
# fill = c("Yellow2", "Green2")
#   18951 and 19240
reps18951 = grep ("18951", colnames(v3$E) ) 
reps19240 = grep ("19240", colnames(v3$E))
barplot( v3$E[2495,reps18951] - min(v3$E[2495,reps18951]-1),
         col= c("Green2", "Green2" ,"Yellow2" ,"Yellow2", "Yellow2" ) ) 
barplot( v3$E[2495,reps19240] - min(v3$E[2495,reps19240] - 1), 
         col= c("Green2", "Green2" , "Green2", "Green2" , "Green2", "Yellow2" ,"Yellow2", "Yellow2" ) ) 


