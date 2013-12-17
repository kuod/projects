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
#Plus strand
starts_Plus = []
stops_Plus = []
genes_Plus = []
chromosomes_Plus = []

#Get rnaSeq values 2kb from promotor
for g in genesRaw:
	if not g[0] == 'chrY' and not g[0] == 'chrM':
		if g[6] == '-':
			chromosomes_Plus.append(g[0])
			if g[3] > 2000:
				starts_Plus.append(g[3]-2000)
			else:
				starts_Plus.append(g[3])
			stops_Plus.append(g[4])
			curGene = g[8].split()[1].translate(None, '";')
			genes_Plus.append(curGene)

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
k562FilePlusRep1 = '/cbio/grlab/share/databases/encode/k562/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqK562CytosolPapPlusRawSigRep1.bigWig'
k562FilePlusRep2 = '/cbio/grlab/share/databases/encode/k562/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqK562CytosolPapPlusRawSigRep2.bigWig'
gm12878FilePlusRep1 = '/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapPlusRawSigRep1.bigWig'
gm12878FilePlusRep2 = '/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapPlusRawSigRep2.bigWig'
helaFilePlusRep1 = '/cbio/grlab/share/databases/encode/hela/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHelas3CytosolPapPlusRawSigRep1.bigWig'
helaFilePlusRep2 = '/cbio/grlab/share/databases/encode/hela/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHelas3CytosolPapPlusRawSigRep2.bigWig'
hepg2FilePlusRep1 ='/cbio/grlab/share/databases/encode/hepg2/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHepg2CytosolPapPlusRawSigRep1.bigWig'
hepg2FilePlusRep2 ='/cbio/grlab/share/databases/encode/hepg2/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHepg2CytosolPapPlusRawSigRep2.bigWig'
huvecFilePlusRep1 = '/cbio/grlab/share/databases/encode/huvec/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHuvecCytosolPapPlusRawSigRep3.bigWig'
huvecFilePlusRep2 = '/cbio/grlab/share/databases/encode/huvec/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqHuvecCytosolPapPlusRawSigRep4.bigWig'

k562FilePlusRep1F = loadBw(k562FilePlusRep1)
k562FilePlusRep2F = loadBw(k562FilePlusRep2)

gm12878FilePlusRep1F = loadBw(gm12878FilePlusRep1)
gm12878FilePlusRep2F = loadBw(gm12878FilePlusRep2)

helaFilePlusRep1F = loadBw(helaFilePlusRep1)
helaFilePlusRep2F = loadBw(helaFilePlusRep2)

hepg2FilePlusRep1F = loadBw(hepg2FilePlusRep1)
hepg2FilePlusRep2F = loadBw(hepg2FilePlusRep2)

huvecFilePlusRep1F = loadBw(huvecFilePlusRep1)
huvecFilePlusRep2F = loadBw(huvecFilePlusRep2)


#Plus strand
PlusOutputName = resultsFolder + 'PlusOutput_withRep.txt'
PlusOutput = open(PlusOutputName, "w+")
for i in range(len(starts_Plus)):
	chrom = chromosomes_Plus[i]
	start = starts_Plus[i]
	stop = stops_Plus[i]
	k562Rep1Arr = k562FilePlusRep1F.get_as_array(chrom,start,stop)
	k562Rep1Sum = rmNanRetSum(k562Rep1Arr)

	k562Rep2Arr = k562FilePlusRep2F.get_as_array(chrom,start,stop)
	k562Rep2Sum = rmNanRetSum(k562Rep2Arr)
	
	gm12878Rep1Arr = gm12878FilePlusRep1F.get_as_array(chrom,start,stop)
	gm12878Rep1Sum = rmNanRetSum(gm12878Rep1Arr)

	gm12878Rep2Arr = gm12878FilePlusRep2F.get_as_array(chrom,start,stop)
	gm12878Rep2Sum = rmNanRetSum(gm12878Rep2Arr)
	
	helaRep1Arr = helaFilePlusRep1F.get_as_array(chrom,start,stop)
	helaRep1Sum = rmNanRetSum(helaRep1Arr)

	helaRep2Arr = helaFilePlusRep2F.get_as_array(chrom,start,stop)
	helaRep2Sum = rmNanRetSum(helaRep2Arr)
	
	hepg2Rep1Arr = hepg2FilePlusRep1F.get_as_array(chrom,start,stop)
	hepg2Rep1Sum = rmNanRetSum(hepg2Rep1Arr)

	hepg2Rep2Arr = hepg2FilePlusRep2F.get_as_array(chrom,start,stop)
	hepg2Rep2Sum = rmNanRetSum(hepg2Rep2Arr)

	huvecRep1Arr = huvecFilePlusRep1F.get_as_array(chrom,start,stop)
	huvecRep1Sum = rmNanRetSum(huvecRep1Arr)

	huvecRep2Arr = huvecFilePlusRep2F.get_as_array(chrom,start,stop)
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

	geneLine = genes_Plus[i] + '\t' + str(vecDist(k562Rep1Sum,k562Rep2Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,gm12878Rep1Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,gm12878Rep2Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,helaRep1Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,helaRep2Sum,'euc')) + '\t' + \
		str(vecDist(k562Rep1Sum,hepg2Rep1Sum,'euc')) + '\t'+ \
		str(vecDist(k562Rep1Sum,hepg2Rep2Sum,'euc')) + '\t'+ \
		str(vecDist(k562Rep1Sum,huvecRep1Sum,'euc')) + '\t'+ \
		str(vecDist(k562Rep1Sum,huvecRep2Sum,'euc')) + '\t'+ \
		str(minValCt) + '\n'
	PlusOutput.write(geneLine)
PlusOutput.close()