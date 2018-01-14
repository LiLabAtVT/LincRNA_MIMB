########################################################################################################################
# By Alex Qi Song
# 01/12/18
# usage: Rscript makeCoexpNet.R <lincRNA_FPKM_file> <protein_coding_FPKM_file> <number_of_cluster> <output_edgelist_file>
#######################################################################################################################

# Process commandline arguments
args = commandArgs(TRUE)
lincAveFPKM <- read.csv(args[1],row.names = 1)
proAveFPKM <- read.csv(args[2],row.names = 1)
k <- args[3]
outFile <- args[4]

# Remove lincRNA with all expression < 0.05
cat("Filtering the FPKM matrices...\n")
lincIDs <- apply(lincAveFPKM > 0.05,1,any)

# Find the 5000 most variable gebes
mostVarIDs <- names(sort(apply(proAveFPKM,1,var),decreasing = TRUE)[1:5000])
cat("Done\n")

# Perform kmeans clustering using the most variable genes and filtered lincRNAs
cat("Performing kmeans clustering...\n")
combFPKM <- rbind(lincAveFPKM[lincIDs,],proAveFPKM[mostVarIDs,])
clust <- kmeans(combFPKM,k,iter.max = 30)
cat("Done\n")

# Find the cluster with most balanced number of lincRNA and protein coding genes
cat("Performing the correlation test for most balanced cluster...\n")
minDiff = 0.5;mostBalancedClust = NA
for(clustID in unique(clust$cluster)){
  clustGenes <- names(clust$cluster[clust$cluster == clustID])
  Diff = abs(0.5 - length(grep("^AT.*",clustGenes))/length(clustGenes))
  if(Diff<minDiff){
    minDiff = Diff
    mostBalancedClust = clustID
  }
}

# Perform correlation test in this cluster. Find significant gene pairs to generate edge list
coGenes <- names(clust$cluster[clust$cluster == mostBalancedClust])
n <- length(coGenes)
r <- cor(t(combFPKM[coGenes,]),use = "pairwise.complete.obs") # compute correlation
t <- (r*sqrt(n-2))/sqrt(1-r^2) # compute t statistic
p <- 1 - pt(t,(n-2)) # compute probability (p-value). Single tail test (coexp > mean(coexp))
pIndex <- which(upper.tri(p),arr.ind = TRUE)
edgeInd <- pIndex[p.adjust(p[pIndex],"bonferroni") < 0.05,] # Adjust pvalues and select pairs with p-value < 0.05
edgeList <- data.frame("Gene1"= rownames(p)[edgeInd[,1]],"Gene2" = colnames(p)[edgeInd[,2]])  # Generate edgelist
cat("Done\n")
write.csv(edgeList,outFile,row.names = F, quote = F)
