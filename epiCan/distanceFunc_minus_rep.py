#Input
#	gene positions
#	list of bigwig file
#Output
#	text file of nearest neighbors based on euclidean distance

import h5py
import os
import numpy as np
from scipy import spatial as spsp
import matplotlib.pyplot as plt
from matplotlib import cm
import time
from bx.bbi.bigwig_file import BigWigFile

#geneFolder contains things with 2kb in promoter
#geneFolder = '/cbio/grlab/share/databases/encode/hdf5_2kb'
#geneFiles = []
#for geneFile in os.listdir(geneFolder):
#	if geneFile.endswith(".hdf5"):
#		geneFiles.append(geneFolder + '/' + geneFile)

#output file
date = time.strftime("%d%m%y/")
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/' + date
if not os.path.exists(resultsFolder):
	os.makedirs(resultsFolder)

#Gene information
geneListFile = "/cbio/grlab/share/databases/genomes/H_sapiens/GENCODE/release18/gencode.v18.annotation.genesPcLncrna.gtf"
#Load gene information in GTF format
genesRaw = np.loadtxt(geneListFile, delimiter='\t',
	dtype={'names':('chromosome', 'annotation','element', 'start', 'end',
	'score', 'strand', 'frame', 'attribute'),'formats': ('S5', 'S10','S10',
	'u4','u4','S1','S1','S1','S100')})
#minus strand
starts_minus = []
stops_minus = []
genes_minus = []
chromosomes_minus = []

#Get rnaSeq values 2kb from promotor
for g in genesRaw:
	if not g[0] == 'chrY' and not g[0] == 'chrM':
		if g[6] == '-':
			chromosomes_minus.append(g[0])
			if g[3] > 2000:
				starts_minus.append(g[3]-2000)
			else:
				starts_minus.append(g[3])
			stops_minus.append(g[4])
			curGene = g[8].split()[1].translate(None, '";')
			genes_minus.append(curGene)

#Distance Functions
def matDist(inputMat, mode):
	if mode == 'euc':
		Y = spsp.distance.pdist(inputMat, 'euclidean')
		return spsp.distance.squareform(Y)
	if mode == '':
		pass
quitInd = False

def vecDist(vec1, vec2, mode):
	if mode == 'euc':
		return np.sqrt(np.square(vec1) + np.square(vec2))
	if mode == '':
		pass

#return sum
def rmNanRetSum(arr):
	arr[np.isnan(arr)] = 0
	return np.sum(arr)

#load bigwig
def loadBw(file):
	 bw = BigWigFile(open(file))
	 return bw

testCellType = 'k562'
compareCellType = ['gm12878','hela','hepg2', 'huvec']

#test cell line
k562FileMinusRep1 = '/cbio/grlab/share/databases/encode/k562/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqK562CytosolPapMinusRawSigRep1.bigWig'
k562FileMinusRep2 = '/cbio/grlab/share/databases/encode/k562/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqK562CytosolPapMinusRawSigRep2.bigWig'
gm12878FileMinusRep1 = '/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapMinusRawSigRep1.bigWig'
gm12878FileMinusRep2 = '/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapMinusRawSigRep2.bigWig'
helaFileMinusRep1 = '/cbio/grlab/share/databases/encode/hela/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHelas3CytosolPapMinusRawSigRep1.bigWig'
helaFileMinusRep2 = '/cbio/grlab/share/databases/encode/hela/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHelas3CytosolPapMinusRawSigRep2.bigWig'
hepg2FileMinusRep1 ='/cbio/grlab/share/databases/encode/hepg2/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHepg2CytosolPapMinusRawSigRep1.bigWig'
hepg2FileMinusRep2 ='/cbio/grlab/share/databases/encode/hepg2/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHepg2CytosolPapMinusRawSigRep2.bigWig'
huvecFileMinusRep1 = '/cbio/grlab/share/databases/encode/huvec/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHuvecCytosolPapMinusRawSigRep3.bigWig'
huvecFileMinusRep2 = '/cbio/grlab/share/databases/encode/huvec/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHuvecCytosolPapMinusRawSigRep4.bigWig'

k562FileMinusRep1F = loadBw(k562FileMinusRep1)
k562FileMinusRep2F = loadBw(k562FileMinusRep2)

gm12878FileMinusRep1F = loadBw(gm12878FileMinusRep1)
gm12878FileMinusRep2F = loadBw(gm12878FileMinusRep2)

helaFileMinusRep1F = loadBw(helaFileMinusRep1)
helaFileMinusRep2F = loadBw(helaFileMinusRep2)

hepg2FileMinusRep1F = loadBw(hepg2FileMinusRep1)
hepg2FileMinusRep2F = loadBw(hepg2FileMinusRep2)

huvecFileMinusRep1F = loadBw(huvecFileMinusRep1)
huvecFileMinusRep2F = loadBw(huvecFileMinusRep2)


#minus strand
minusOutputName = resultsFolder + 'minusOutput_withRep.txt'
minusOutput = open(minusOutputName, "w+")
for i in range(len(starts_minus)):
	chrom = chromosomes_minus[i]
	start = starts_minus[i]
	stop = stops_minus[i]
	k562Rep1Arr = k562FileMinusRep1F.get_as_array(chrom,start,stop)
	k562Rep1Sum = rmNanRetSum(k562Rep1Arr)

	k562Rep2Arr = k562FileMinusRep2F.get_as_array(chrom,start,stop)
	k562Rep2Sum = rmNanRetSum(k562Rep2Arr)
	
	gm12878Rep1Arr = gm12878FileMinusRep1F.get_as_array(chrom,start,stop)
	gm12878Rep1Sum = rmNanRetSum(gm12878Rep1Arr)

	gm12878Rep2Arr = gm12878FileMinusRep2F.get_as_array(chrom,start,stop)
	gm12878Rep2Sum = rmNanRetSum(gm12878Rep2Arr)
	
	helaRep1Arr = helaFileMinusRep1F.get_as_array(chrom,start,stop)
	helaRep1Sum = rmNanRetSum(helaRep1Arr)

	helaRep2Arr = helaFileMinusRep2F.get_as_array(chrom,start,stop)
	helaRep2Sum = rmNanRetSum(helaRep2Arr)
	
	hepg2Rep1Arr = hepg2FileMinusRep1F.get_as_array(chrom,start,stop)
	hepg2Rep1Sum = rmNanRetSum(hepg2Rep1Arr)

	hepg2Rep2Arr = hepg2FileMinusRep2F.get_as_array(chrom,start,stop)
	hepg2Rep2Sum = rmNanRetSum(hepg2Rep2Arr)

	huvecRep1Arr = huvecFileMinusRep1F.get_as_array(chrom,start,stop)
	huvecRep1Sum = rmNanRetSum(huvecRep1Arr)

	huvecRep2Arr = huvecFileMinusRep2F.get_as_array(chrom,start,stop)
	huvecRep2Sum = rmNanRetSum(huvecRep2Arr)
	
	disList = [k562Rep2Sum,gm12878Rep1Sum,gm12878Rep2Sum,helaRep1Sum,helaRep2Sum,hepg2Rep1Sum,hepg2Rep2Sum,huvecRep1Sum,huvecRep2Sum]
	minValIdx = disList.index(min(disList))
	if len(set(disList)) == 1:
		minValCt = "all"
	elif minValIdx == 0:
		minValCt = "k562Rep2"
	elif minValIdx == 1:
		minValCt = "gm12878Rep1"
	elif minValIdx == 2:
		minValCt = "gm12878Rep2"
	elif minValIdx == 3:
		minValCt = "helaRep1"
	elif minValIdx == 4:
		minValCt = "helaRep2"
	elif minValIdx == 5:
		minValCt = "hepg2Rep1"
	elif minValIdx == 6:
		minValCt = "hepg2Rep2"
	elif minValIdx == 7:
		minValCt = "huvecRep1"
	elif minValIdx == 8:
		minValCt = "huvecRep2"

	geneLine = genes_minus[i] + '\t' + str(vecDist(k562Rep1Sum,k562Rep2Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,gm12878Rep1Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,gm12878Rep2Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,helaRep1Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,helaRep2Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,hepg2Rep1Sum,'euc')) + '\t'+ \
		str(vecDist(k562Rep1Sum,hepg2Rep2Sum,'euc')) + '\t'+ \
		str(vecDist(k562Rep1Sum,huvecRep1Sum,'euc')) + '\t'+ \
		str(vecDist(k562Rep1Sum,huvecRep2Sum,'euc')) + '\t'+ \
		str(minValCt) + '\n'
	minusOutput.write(geneLine)
minusOutput.close()