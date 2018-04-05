# Step 1 in read mapping is to contruct a mapping index.
# replay "Gnmdir" with a specific directory name
# 
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir Gnmdir \
     --genomeFastaFiles TAIR10_Chr.all.fasta \
     --sjdbGTFfile Araport11_GFF3_genes_transposons.201606.gtf


# Step 2 in read mapping is to actually map reads to the index
# This dataset is paired-end reads
# Mapping requires providing names for two files, representing both ends of reads 

STAR   --runThreadN 16 \
       --genomeDir Gnmdir \
       --readFilesIn SRR3664408.1.fastq SRR3664408.2.fastq \
       --outSAMstrandField intronMotif \
       --outFileNamePrefix output/bam/SRR3664408 \
       --outSAMtype BAM SortedByCoordinate;

