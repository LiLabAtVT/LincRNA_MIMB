# This script use stringtie to construct gene structures
# 

stringtie SRR3664408.bam \
           -o SRR3664408.gtf \
           -p 4 \
           -G Araport11_GFF3_genes_transposons.201606.gtf


# once stringtie was used to analyze each file, use -merge to merge all the GTF files to generate final prediction


stringtie --merge \
          -G Araport11_GFF3_genes_transposons.201606.gtf \
          -o merged.Arath.gtf ../gtf_list.txt


