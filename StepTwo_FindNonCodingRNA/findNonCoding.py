"""
By Alex Qi Song
01/10/18
usage:  python findNonCoding.py [-o OUTFILE] \
        <GTF_for_transcript> \
        <BED_for_pro_coding_genes> \
        <BED_for_transpos_elements> \
        <genome_assembly_fasta>
example: python findNonCoding.py merged.Arath.gtf Araport11_protein_coding.201606.bed Araport11_transposable_element_gene.201606.bed TAIR10_Chr.all.fasta > nonCodingTranscripts
"""
__author__ = 'Alex Qi Song'
from collections import defaultdict,namedtuple
from Bio import SeqIO,BiopythonWarning
from Bio.Seq import Seq
from pandas import read_table
from getpass import getuser
from time import strftime
from socket import gethostname
import sys
import os
import argparse
import warnings
import myGTF

# Suppress warnings
warnings.simplefilter('ignore',BiopythonWarning)
# Define region class
reg = namedtuple('region',['start','end'])

# Parse input arguments and print help message
def parseArgs(args):
    parser = argparse.ArgumentParser(prog = args[0],description='Identify non-coding genes',epilog = "Author: Alex Qi Song alexsong@vt.edu")
    parser.add_argument("-o","--out",dest = "out",metavar = "OUTFILE", help="Output file name. If provided, this will suppress the output of gene list on screen",type=str,required = True)
    parser.add_argument("gtf",metavar = "GTF_file", help="GTF file for transcripts of interest",type=str)
    parser.add_argument("bedPro",metavar = "Pro_BED", help="BED file for protein coding genes",type=str)
    parser.add_argument("bedTrans",metavar = "Transpos_BED", help="BED fle for transposable elements",type=str)
    parser.add_argument("fasta",metavar = "FASTA", help="FASTA file for chromosome sequences",type=str)
    parsedArgs = parser.parse_args(args[1:])

    # Check if all input files exist
    fileNames = [val for opt,val in parsedArgs.__dict__.items() if opt != "out"]
    for fileName in fileNames:
        try:
            if not os.path.exists(fileName):
                raise Exception('Error: file does not exist: {}'.format(fileName))
        except Exception as excOut:
            parser.print_usage()
            print(excOut)
            sys.exit(1)

    return(parsedArgs)

'''
Annotation info for transcirpts
Class attributes:
    ts2len: transcript ID => transcript length
    ts2pos: transcript ID => [[exon1Start,exon1End],[exon2Start,exon2End]...]
    ts2chr: transcript ID => chromosome name e.g. 'Chr1'
    ts2seq: transcript ID => transcript sequence (concatenated exon sequences)
    gene2ts: gene ID => transcript ID
Class methods:
    parseGTF: parsing gtf file, sending the data to self.ts2len, self.ts2pos
    self.ts2chr
    parseFasta: parsing chromosome fasta, sending the data to ts2chr
'''
class annoInfo:
    def __init__(self):
        self.ts2len = defaultdict(int)
        self.ts2pos = defaultdict(list)
        self.gene2ts = defaultdict(set)
        self.ts2chr = dict()
        self.ts2seq = dict()

    def parseGTF(self,GTFfile):
        for line in myGTF.lines(GTFfile):
            geneID = line['gene_id']
            tsID = line['transcript_id']
            chrname = line['seqname']
            startPos = int(line['start'])
            endPos = int(line['end'])
            if not geneID.startswith("AT") and (line['feature'] == 'exon' or line['feature'] == '5UTR' or line['feature'] == '3UTR'):
                self.ts2len[tsID] += (endPos-startPos+1)
                self.ts2pos[tsID].append(reg._make([startPos,endPos]))
                self.gene2ts[geneID].add(tsID)
            self.ts2chr[tsID]=chrname

    # Convert fasta file into a dictionary key=> seqID, value=> seq
    def parseFasta(self,fastaFile):
        chr2seq = {seq_record.id:seq_record.seq for seq_record in SeqIO.parse(fastaFile, "fasta")}
        for tsID,tsChr in self.ts2chr.items():
            self.ts2seq[tsID] = Seq("".join([str(chr2seq[tsChr][start:end+1]) for start,end in self.ts2pos[tsID]]))
'''
Gene regions to be searched against
Class attributes:
    proCodingRegion: chromosome name => [[gene1Start,gene1End],[gene2Start,gene2End]...]
    transposRegion: chromosome name => [[gene1Start,gene1End],[gene2Start,gene2End]...]
Class methods:
    parseAllBED: parsing all bed files, sending all data into self.proCodingRegion and
    self.transposRegion
'''
class regions:
    def __init__(self): # To remove the duplicated regions, the initial regions were put in a set
        self.proCodingRegion = defaultdict(set)
        self.transposRegion = defaultdict(set)

    def parseAllBED(self,proBedFile,transposBedFile):
        self._parseBED(proBedFile,self.proCodingRegion)
        self._parseBED(transposBedFile,self.transposRegion)

    def _parseBED(self,bedFile,attr):
        for row in read_table(bedFile).iterrows():
            chrName,startPos,endPos,geneID = row[1][0],row[1][1],row[1][2],row[1][3]
            attr[chrName].add(reg._make([startPos,endPos]))

        # sort gene postitions for each chromosome and convert to list
        for chrName,allPos in attr.items():
            attr[chrName] = sorted(allPos,key = lambda x: x.start)
            attr[chrName].sort(key=lambda x: x.end) # Break the tie by end pos

# Find orf length from six frameshifts
def maxOrfLen(seq):
    seqRev = seq.reverse_complement()
    allSeqs = [seq,seq[1:],seq[2:],seqRev,seqRev[1:],seqRev[2:]]
    return(max([orfLen(seq) for seq in allSeqs]))

def orfLen(seq):
    startCodon = Seq("ATG")
    stopCodon = set([Seq("TAG"),Seq("TAA"),Seq("TGA")])
    startPos = endPos = -1
    for i in range(len(seq)//3):
        if seq[3*i:3*i+3] == startCodon:
            startPos = 3*i
        if seq[3*i:3*i+3] in stopCodon:
            endPos = 3*i
            break
    if startPos != -1 and endPos != -1:
        return(endPos-startPos)
    else:
        return(0)

'''
See if the given position is close to a gene.
Return 0 for not close to any gene, 1 for overlapping with gene(s),
2 for close to a non-overlapping gene.
'''
def closeToGene(tsStart,tsEnd,chrName,dist,genePos):
    if chrName not in genePos:
        return(3) # No chromosome info available
    lIndex1,rIndex1 = biSearch(tsStart,chrName,genePos)
    lIndex2,rIndex2 = biSearch(tsEnd,chrName,genePos)
    if lIndex1 == rIndex1 or lIndex2 == rIndex2 or rIndex2 - lIndex1 > 1:
        return(1) # Ovelap with existing gene

    # Get end position of left closest gene and start position of right closest gene
    lPos = 0 if lIndex1 == -1 else genePos[chrName][lIndex1].end
    rPos = float('Inf') if rIndex2 == len(genePos[chrName]) else genePos[chrName][rIndex1].start
    if lPos > tsStart - dist or rPos < tsEnd + dist:
        return(2) # Close to an existing gene
    return(0) # None of the above apply

# Binary search for the closest gene
def biSearch(pos,chrName,genePos):
    left = 0
    right = len(genePos[chrName])-1
    while left <= right:
        mid = (left+right)/2
        if genePos[chrName][mid].start > pos:
            right = mid-1
        elif genePos[chrName][mid].end < pos:
            left = mid+1
        elif genePos[chrName][mid].start <= pos and genePos[chrName][mid].end >= pos:
            return(mid,mid)
    if genePos[chrName][mid].start > pos:
        return(mid-1,mid)
    else:
        return(mid,mid+1)

#Main function. Non-coding transcripts are selected based on four criteria
def main(args):
    sys.stdout.write(strftime("%d/%m/%Y %H:%M:%S {0}@{1}: Parsing input files ...\n").format(getuser(),gethostname()))
    parsedArgs = parseArgs(args)
    anno = annoInfo()
    anno.parseGTF(parsedArgs.gtf)
    anno.parseFasta(parsedArgs.fasta)
    geneRegions = regions()
    geneRegions.parseAllBED(parsedArgs.bedPro,parsedArgs.bedTrans)
    sys.stdout.write(strftime("%d/%m/%Y %H:%M:%S {0}@{1}: Done\n").format(getuser(),gethostname()))
    outFile = open(parsedArgs.out,"w")
    '''
    Process transcripts for each gene. If all transcripts of a gene are non-coding
    transcripts. The gene will be considerred as a non-coding gene
    '''
    totalCount = len(anno.gene2ts)
    curCount = 0
    sys.stdout.write(strftime("%d/%m/%Y %H:%M:%S {0}@{1}: Processing genes ...\n").format(getuser(),gethostname()))
    for geneID,tsIDs in anno.gene2ts.items():
        curCount += 1
        sys.stdout.write("\r{0}/{1} ({2}%) finished ".format(curCount,totalCount,round(100*curCount/float(totalCount),2)))
        sys.stdout.flush()
        isNonCoding = True
        for tsID in tsIDs:
            tsLen = anno.ts2len[tsID]
            tsStart = min([pos.start for pos in anno.ts2pos[tsID]])
            tsEnd = max([pos.end for pos in anno.ts2pos[tsID]])
            chrName = anno.ts2chr[tsID]
            '''
            Four criteria for non-coding transcript:
                1) exon region length < 200
                and 2) maximum orf length/3 < 100
                and 3) not close to any protein-coding gene within +/- 500 bp
                and 4) not overlapping with any transposable elements
            '''
            if tsLen < 200 or maxOrfLen(anno.ts2seq[tsID])/3 > 100:
                isNonCoding = False
                break
            else:
                closeToPro = closeToGene(tsStart,tsEnd,chrName,500,geneRegions.proCodingRegion)
                closeToTranspos = closeToGene(tsStart,tsEnd,chrName,0,geneRegions.transposRegion)
                if closeToPro or closeToTranspos in [0,2]:
                    isNonCoding = False
                    break
        if isNonCoding:
            outFile.write(geneID+"\n")
    sys.stdout.write("\n")
    sys.stdout.write(strftime("%d/%m/%Y %H:%M:%S {0}@{1}: All done!\n").format(getuser(),gethostname()))
    outFile.close()

if __name__ == "__main__":
    main(sys.argv)
