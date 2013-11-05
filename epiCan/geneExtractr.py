import numpy as np
import sys, os, csv
import time

from bx.bbi.bigwig_file import BigWigFile

bwTestFile = "/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapMinusRawSigRep1.bigWig"
geneListFile = "/cbio/grlab/share/databases/genomes/H_sapiens/GENCODE/release18/gencode.v18.annotation.chr21.gtf"
methylFile = "/cbio/grlab/share/databases/encode/k562/methyl/wgEncodeHaibMethylRrbsK562HaibSitesRep1_chr21.bed"


#Load gene information in GTF format
genesRaw = np.loadtxt(geneListFile, delimiter='\t', dtype={'names':('chromosome', 'annotation', 
	'element', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
	'formats': ('S5', 'S10','S10','u4','u4','S1','S1','S1','S100')})

#methylation
#methylRaw = np.loadtxt(methylFile, delimiter='\t', usecols=(0,1,2,4), dtype={'names':('chrom', 'chromStart', 'chromEnd',
#	'score'), 'formats': ('S5', 'u4','u4','u4')}) 

#methylRaw = np.loadtxt(methylFile, delimiter='\t', usecols=[1,4])



#methylMed = np.median()

#methylSites = []



#subset based on score
for m in methylRaw:
	if m[3] >  


#Parse out start stop
starts = []
stops = []
genes = []
for g in genesRaw:
	starts.append(g[3]-2000)
	stops.append(g[4]+2000)
	curGene = g[8].split()[1].translate(None, '";')
	genes.append(curGene)

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




	temp = np.hstack((posMat.T,geMat.T))
	fileName = '/cbio/grlab/projects/epiCan/results/chr21/' + genes[i] 
	numpy.savez(fileName, temp)
end = time.time()
print end - start
	#methylation bed













#Get methylation patterns



curGene = bw.get_as_array(chrom, start, end)


def load(file):
	bw = BigWigFile(open(file))
	return bw


def geneFeatExtractr (chrom, start, histoneFile, geneExpFile):

