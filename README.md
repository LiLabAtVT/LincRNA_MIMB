# LincRNA_MIMB
findNonCoding.py is a python script for identifying the non-coding genes and makeCoexpNet.R is an R script for identifying the co-expression network for lncRNAs and protein-coding genes
## Identification of non-coding genes
### Installation of required Python packages
Before running the script, two Python packages, Pandas and Biopython need to be installed. The two packages can be conveniently installed using pip or conda. Please refer to https://pip.pypa.io/en/stable/installing/ for installing or upgrading pip, and https://conda.io/docs/user-guide/install/index.html for installing conda.
To install Pandas and Biopython through pip, use the two commands below:
```shell
$ pip install pandas
$ pip install biopython
```
Alternatively, if you prefer using conda to install and manage the packages, use the two commands below:
```shell
$ conda install pandas
$ conda install -c anaconda biopython
```
All package dependencies will be automatically handled by pip or conda.
### Run python script
Download the two python scripts (myGTF.py and findNonCoding.py) from the [github repository](https://github.com/LiLabAtVT/LincRNA_MIMB) myGTF.py defines functions for parsing the GTF file. findNonCoding.py uses these functions to parse input file and find non-coding genes. findNonCoding.py requires four input files: 
1. a GTF format file containing annotations for genes of interest
2. BED format file that stores annotations for protein-coding genes
3. BED format file that stores annotations for transposable elements
4. FASTA format file for assembled chromosome sequences.

findNonCoding.py will ignore any transcript ID with the format of “ATXGXXXXX” and take other transcripts to check if they are the non-coding transcripts. The output file of findNonCoding.py is a list of identified non-coding gene IDs. Gene is considered as non-coding gene when all isoforms of that gene are non-coding transcripts. The non-coding transcripts are selected from the input GTF file using the four criteria:
1. Mature transcript length > 200 bps (including UTRs and exons, excluding introns).
2. Longest open reading frame does not encode more than 100 amino acids.
3. Transcribed regions do not overlap with any transposable elements
4. Transcribed regions are 500 bp away from any protein coding gene.
Before running the script, make sure myGTF.py and findNonCoding.py are put into the same directory, otherwise findNonCoding.py will fail to load the required functions.

`-h` / `--help` option can be used to display usage and help information for findNonCoding.py:
```shell
$ python findNonCoding.py --help
usage: ../src/findNonCoding.py [-h] -o OUTFILE GTF_file Pro_BED Transpos_BED FASTA

Identify non-coding genes

positional arguments:
  GTF_file              GTF file for transcripts of interest
  Pro_BED               BED file for protein coding genes
  Transpos_BED          BED fle for transposable elements
  FASTA                 FASTA file for chromosome sequences

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --out OUTFILE
                        Output file name
```
Below is an example of running the script:
```shell
$ python findNonCoding.py -o out \
merged.Arath.gtf \ 
Araport11_protein_coding.201606.bed \
Araport11_transposable_element_gene.201606.bed \
TAIR10_Chr.all.fasta
```
`-o` \ `--out` specifies the output file name. In this example, the output gene list was saved into a file named ‘out’.
## Co-expression analysis for lncRNAs and protein-coding genes
### Download R script
The R script can be downloaded from [github repository](https://github.com/LiLabAtVT/LincRNA_MIMB). No required packages need to be installed to run the script.
### Run R script
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
