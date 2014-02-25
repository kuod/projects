#Input
#   epigenetic bigwig file 
#   
#Output
#   text file of signal over gene positions

#Systems Imports
import os
import sys
from sys import argv
import time
from bx.bbi.bigwig_file import BigWigFile

#Numpy
import numpy as np

timeStart = time.time()

#Example
#python geneExpSums /cbio/shared/encode/UCSC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqHelas3CellPapMinusRawSigRep1.bigWig
script, bigWigFp  = argv

fpStart = 'histone'
#Check Starting path
if fpStart.upper() not in bigWigFp.upper():
    print "Full filepath not found. Exiting."
    sys.exit()
#set strand

bigWigNm = bigWigFp.split('/')[7].split('.')[0]

#print 'Working on ' + bigWigNm
#output file
date = time.strftime("%d%m%y/")
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/epiSums/'
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)

#start printing output 
outputFilename = resultsFolder + bigWigNm + '_DistOutput.txt'
#if os.path.exists(outputFilename):
#    print 'file completed, exiting'
#    sys.exit()
    
output = open(outputFilename, "w+")

#Gene information
geneListFile = "/cbio/grlab/share/databases/genomes/H_sapiens/GENCODE/release18/gencode.v18.annotation.genesPcLncrna.gtf"
#Load gene information in GTF format
genesRaw = np.loadtxt(geneListFile, delimiter='\t',
    dtype={'names':('chromosome', 'annotation','element', 'start', 'end',
    'score', 'strand', 'frame', 'attribute'),'formats': ('S5', 'S10','S10',
    'u4','u4','S1','S1','S1','S100')})

#all genes
starts = []
stops = []
genes = []
chromosomes = []

#Get rnaSeq values 2kb from promotor
for strandSym in ['+','-']:
    for g in genesRaw:
        if not g[0] == 'chrY' and not g[0] == 'chrM':
            if g[6] == strandSym:
                chromosomes.append(g[0])
                if g[3] > 2000:
                    starts.append(g[3]-2000)
                else:
                    starts.append(g[3])
                stops.append(g[4])
                curGene = g[8].split()[1].translate(None, '";')
                genes.append(curGene)

### Functions
#return sum
def retSum(arr):
    arr[np.isnan(arr)] = 0
    return np.sum(arr)

#load bigwig
def loadBw(file):
     bw = BigWigFile(open(file))
     return bw

#load file
bigWigFile = loadBw(bigWigFp)

timeStart = time.time()
strLong = 'Gene' + '\t' + 'Sum' + '\n'
for i in range(len(genes)):
    snglChrom = chromosomes[i]
    snglStart = starts[i]
    snglStop = stops[i]

    #get array
    bwArr = bigWigFile.get_as_array(snglChrom,snglStart,snglStop)
    bwSum = retSum(bwArr)

    #
    geneStr = genes[i] + '\t' + str(bwSum) + '\n'
    strLong = strLong + geneStr

    #if (i + 1) % 1000 == 0:
        #print 'Completed 1000 Genes'
output.write(strLong)
output.close()
timeEnd = time.time()
print 'File written to: ' + outputFilename

