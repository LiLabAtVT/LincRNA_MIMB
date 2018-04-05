# Change working directory to the LincRNA_ATH directory
cd LincRNA_ATH/software

# Download and install STAR alignment software
wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz
tar -xzf 2.5.2b.tar.gz
STAR-2.5.2b/bin/Linux_x86_64_static/STAR --version

# Download and install featureCounts software
wget https://sourceforge.net/projects/subread/files/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz/download 
tar -zxvf download
subread-1.5.1-Linux-x86_64/bin/featureCounts -v

# Download and install StringTie Software
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz 
tar -zxvf stringtie-1.3.3b.Linux_x86_64.tar.gz
stringtie-1.3.3b.Linux_x86_64/stringtie --version

# Download and install SRA toolkit
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xzf sratoolkit.current-centos_linux64.tar.gz
./sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --version

# The following commands have to be executed under R interactive session to install DESeq2 and edgeR
#
# source('https://bioconductor.org/biocLite.R')
# biocLite('DESeq2')
# biocLite('edgeR')
#
#

# Install Python packages
pip install pandas
pip install biopython

# Alternative way to install Python packages
#
# conda install pandas
# conda install -c anaconda biopython


