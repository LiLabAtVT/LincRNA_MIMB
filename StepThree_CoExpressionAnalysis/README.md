# Step 3. Co-expression analysis
This is the third step in this analysis, we will use expression data from RNA-seq experiments for both coding and non-coding RNA to perform co-expression analysis.
## Download R script
The R script can be downloaded from [github repository](https://github.com/LiLabAtVT/LincRNA_MIMB). No required packages need to be installed to run the script.
## Run R script
makeCoexpNet.R needs two input files: 
1. A FPKM file containing the expression level of lncRNAs
2. a FPKM file containing the expression level of protein-coding genes.
Additionally, one parameter k which specifies the number of desired clusters, and one argument giving the output file name should be provided.

makeCoexpNet.R will process the input files and construct co-expression network by the following four steps:
1. Remove those having expression < 0.05 across all columns (conditions) in lncRNA FPKM file.
2. For protein-coding gene FPKM file, calculate variance of each gene, keep the top 5000 genes with greatest variance.
3. Perform K-means clustering using filtered data set of lncRNAs and protein-coding genes, number of cluster is specified by command line argument. 
4. Find the cluster with most balanced number of lncRNA and protein-coding genes. Calculate Pearson Correlation Coefficient (PCC) for each gene pair, and p values for the PCC using Student’s t-distribution. Select gene pairs with p values < 0.05. An edge list of these gene pairs will then be generated and saved into the output file.
Below is an example of running the Rscript:  
```shell
$ Rscript makeCoexpNet.R AllLincAveFPKM.v4.1.csv \
aveFPKM.pro.v4.csv \
60 \
edgeList.csv
```
This will generate the co-expresssion network and save the edge list into edgeList.csv
