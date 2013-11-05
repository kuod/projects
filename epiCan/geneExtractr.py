import numpy as np
import sys, os, csv

from bx.bbi.bigwig_file import BigWigFile

bwTestFile = "/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapMinusRawSigRep1.bigWig"
geneListFile = "/cbio/grlab/share/databases/genomes/H_sapiens/GENCODE/release18/gencode.v18.annotation.chr21.gtf"

#Load gene information in GTF format
genesRaw = np.loadtxt(geneListFile, delimiter='\t', dtype={'names':('chromosome', 'annotation', 
	'element', 'start', 'end', 'score', 'strand', 'frame', 'attribute'),
	'formats': ('S5', 'S10','S10','u4','u4','S1','S1','S1','S100')})

#Parse out start stop
starts = []
stops = []
genes = []
for g in genesRaw:
	starts.append(g[3]-2000)
	stops.append(g[4]+2000)
	curGene = g[8].split()[1].translate(None, '";')
	genes.append(curGene)

#genes now indexed
pos = []
for i in range(len(genes)):
	pos = np.mat(np.arange(starts[i], stops[i], 1)).T




for s in range(len(starts)):
	print starts[s]
	print stops[s]


bw = load(bwTestFile)




curGene = bw.get_as_array(chrom, start, end)


def load(file):
	bw = BigWigFile(open(file))
	return bw


def geneFeatExtractr (chrom, start, histoneFile, geneExpFile):

