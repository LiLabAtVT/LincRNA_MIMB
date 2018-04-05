# This step allows counting number of reads for each gene including protein coding genes and lincRNAs

featureCounts -T 8 \ 
              -t exon \ 
              -g gene_id \ 
              -p \ 
              -a $GTF \ 
              -o out.readcount.txt \ 
              $mappingdir/$file1  


