import numpy as np
import sys, os, csv

from bx.bbi.bigwig_file import BigWigFile

bwTestFile = "/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapMinusRawSigRep1.bigWig"
geneListFile = "/cbio/grlab/share/databases/encode/genes/gencodeV12/wgEncodeGencodeCompV12Chr21_5col.gp"

#Load gene information
genes = np.loadtxt(geneListFile, dtype={'names':('transcript', 'chromosome', 'strand', 'start', 'end'),
	'formats': ('S15', 'S10','S1','u4','u4')})

#Parse out start stop
starts = []
stops = []
both = []
for g in genes:
	starts.append(g[3])
	stops.append(g[4])
	both.append([g[3],g[4]])

for s in range(len(starts)):
	print starts[s]
	print stops[s]


bw = load(bwTestFile)




curGene = bw.get_as_array(chrom, start, end)


def load(file):
	bw = BigWigFile(open(file))
	return bw


def geneFeatExtractr (chrom, start, histoneFile, geneExpFile):

