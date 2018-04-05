# these are the main steps for generating FPKM file.

library(DESeq2)
library(edgeR)

dds  = DESeqDataSetFromMatrix(countData = InputDF2,
                              colData = sampleInfo, 
                              design = ~Treatment)
dds<-DESeq(dds)

dds.EstSizFac <- estimateSizeFactors(dds)

norm.count=counts(dds.EstSizFac, normalized=TRUE)

AnnoData <- read.table(GeneidLengthFile, 
                       sep='\t',
                       as.is=T, 
                       header=T)

norm.count.DGEList <- DGEList(counts=norm.count,
                              genes = AnnoData[,c("Geneid","Length")])
    
norm.rpkm <- rpkm(norm.count.DGEList,
                  norm.count.DGEList$genes$Length)

