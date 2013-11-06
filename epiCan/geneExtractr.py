import numpy as np
import sys, os, csv
import time

from bx.bbi.bigwig_file import BigWigFile

#set test file paths
bwTestFile = "/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapMinusRawSigRep1.bigWig"
geneListFile = "/cbio/grlab/share/databases/genomes/H_sapiens/GENCODE/release18/gencode.v18.annotation.chr21.gtf"
methylFile = "/cbio/grlab/share/databases/encode/k562/methyl/wgEncodeHaibMethylRrbsK562HaibSitesRep1_chr21.bed"

#TODO Take in path of HM folder  
hm_1 = "/cbio/grlab/share/databases/encode/k562/histonePeaks/wgEncodeBroadHistoneK562H3k4me3StdAln_chr21.bed"


#Load gene information in GTF format
genesRaw = np.loadtxt(geneListFile, delimiter='\t', dtype={'names':('chromosome', 'annotation', 
	'element', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
	'formats': ('S5', 'S10','S10','u4','u4','S1','S1','S1','S100')})

#Get methylation patterns
#methylation in bed
#methylRaw = np.loadtxt(methylFile, delimiter='\t', usecols=(0,1,2,4), dtype={'names':('chrom', 'chromStart', 'chromEnd',
#	'score'), 'formats': ('S5', 'u4','u4','u4')}) 
#methylRaw = np.loadtxt(methylFile, delimiter='\t', usecols=[1,4])
#methylMed = np.median()
#methylSites = []
#subset based on score
#for m in methylRaw:
#	if m[3] >  



#Parse start stop and gene names
starts = []
stops = []
genes = []
for g in genesRaw:
	starts.append(g[3]-2000)
	stops.append(g[4]+2000)
	curGene = g[8].split()[1].translate(None, '";')
	genes.append(curGene)

#turn histone track into coordinate pairs
hm1Arr = loadHm(hm_1)
hm1Coord = []
for h in hm1Arr:
	hm1Coord.append([h[1], h[2]])


#Summarize gene expression for length of gene?
bw = load(bwTestFile)

#generate gene by gene matrix
bigList = []
start = time.time()
for i in range(len(genes)):
	#pos index
	posMat = []
	posMat = np.mat(np.arange(starts[i], stops[i], 1))
	#geneExp summary

	#ge from bigwig
	geMat = []
	geMat = np.mat(bw.get_as_array(chromosome, starts[i], stops[i]).T)
	geMat[np.isnan(geMat)] = 0

	hmTrack = htGen(starts[i], stops[i], hm1Coord)
	hmTrack = np.mat(hmTrack) 

	temp = np.hstack((posMat.T,geMat.T, hmTrack.T))
	fileName = '/cbio/grlab/projects/epiCan/results/chr21/' + genes[i] 
	numpy.savez(fileName, temp)
	#methylation bed
end = time.time()
print end - start


#Needs htCoord
def htGen(start, stop, htCoord):
	hTrack = np.zeros(stop-start)
	
	htCoord = np.mat(htCoord)

	for posCnt in range(len(hTrack)):
		pos = posCnt + start
		htSourceMin = htCoord[pos > htCoord[:,0]][0,0]
		htSourceMax = htCoord[pos > htCoord[:,0]][0,1]
		if pos > htSourceMin and pos < htSourceMax:
			hTrack[pos] = 1
			print 'aok'
	return hTrack



def loadHm(file):
	hmArray = np.loadtxt(file, delimiter='\t', usecols=(0,1,2,4), dtype={'names':('chrom', 'chromStart', 'chromEnd',
	'score'), 'formats': ('S5', 'u4','u4','u4')}) 
	return hmArray

curGene = bw.get_as_array(chrom, start, end)

def load(file):
	bw = BigWigFile(open(file))
	return bw


