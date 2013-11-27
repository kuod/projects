#Input
#	training matrix of values including ge
#	test vector of comparable gene expression
#Output
#	predicted histone track based on euclidian distance
import h5py
import os
import numpy as np
from numpy import linalg as la
from scipy import spatial as spsp
import matplotlib.pyplot as plt
import time

#All gene files in geneFiles
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/131127/'
geneFolder = '/cbio/grlab/home/dkuo/temp'
geneFiles = []
os.chdir( geneFolder )
geDir = os.getcwd() 
for files in os.listdir("."):
	if files.endswith(".hdf5"):
		geneFiles.append(geDir + '/' + files)


def dist(inputMat, mode):
	if mode == 'euc':
		Y = spsp.distance.pdist(inputMat, 'euclidean')
		return spsp.distance.squareform(Y)
		
#getsGeneNames
uniqGenes = []
for geneFile in geneFiles:
	uniqGenes.append(geneFile.split('/')[-1].split('_')[0])
uniqGenes = np.unique(uniqGenes)

start = time.time()
for uGene in uniqGenes:
	#find the files that correspond to that gene
	genesByCellTypes = filter(lambda x: uGene in x, geneFiles)
	gexMat = []
	geMatName = []
	for geneOneCellType in genesByCellTypes:
		#print 'now on ' + geneCellTypeFile
		f = h5py.File(str(geneOneCellType), 'r')
		#print f.values()
		g = f.__getitem__(uGene)
		#get gene expression indices
		indices = [i for i, s in enumerate(f.attrs.values()) if 'Gene expression' in s]
		
		if len(gexMat) == 0:
			gexMat = np.mat(g[:,idx]).T
			geMatName = []
		else:
			for idx in indices:
				gexMat = np.hstack((gexMat, np.mat(g[:,idx]).T))
				#print gexMat.shape
				geMatName.append(f.attrs.values()[idx])
		f.close()
		#print geMatName
		gexMat2 = gexMat[:,1:].T
		#print gexMat2.shape
		distMat = dist(gexMat2,'euc')
		#print distMat.shape
		#print 'oneCellTypeDone'
	fig = plt.figure()
	fig.suptitle(uGene,fontsize=14, fontweight='bold')

	ax = fig.add_subplot(111)
	cax = ax.matshow(distMat, interpolation='nearest')
	#ax.set_xticklabels(['']+geMatName, rotation=45)
	ax.set_yticklabels(['']+geMatName)
	figName = resultsFolder + uGene + '_euc.png'
	fig.savefig(figName, bbox_inches='tight')
	print 'onto the next one'
	pyplot.clf()
	#plt.show()
end = time.time()
print 'Total time: ' + str(end - start)








	#grab the relevant gene expression tracks
	#calculate distances
	#save output of distances and show assignment based on filename

	if uGene in geneFiles:










#output: for each gene, show the distances from the 





def eucDist()


#hardcoded
kVal = 2

#def kNn(kVal, test, training, distFunc):

#def dotProd(test, training):
#	return np.dot(test, training)

def eucDist(input, training):
	d = np.sum((input-trainin)^2)
	return d

	#return la.norm(test-training)

def sqeucDist(input,training):
	return spsp.distance.sqeuclidian(test, training)








#for s in 


#import matplotlib.pyplot as plt
#plt.plot([1,2,3,4])
#plt.ylabel('some numbers')
#plt.show()