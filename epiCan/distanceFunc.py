#Input
#	training matrix of values including ge
#	test vector of comparable gene expression
#Output
#	predicted histone track based on euclidian distance
import h5py 
import numpy as np
from numpy import linalg as la
from scipy import spatial as spsp

#All gene files in geneFiles
geneFolder = '/cbio/grlab/home/dkuo/temp'
geneFiles = []
os.chdir( geneFolder )
geDir = os.getcwd() 
for files in os.listdir("."):
	if files.endswith(".hdf5"):
		geneFiles.append(geDir + '/' + files)

#getsGeneNames
uniqGenes = []
for f in geneFiles:
	uniqGenes.append(f.split('/')[-1].split('_')[0])

uniqGenes = np.unique(uniqGenes)

for uGene in uniqGenes:
	#find the files that correspond to that gene
	geneByCellTypeList = filter(lambda x: uGene in x,geneFiles)

	gexMat = []
	geMatName = []
	for geneCellTypeFile in geneByCellTypeList:
		print 'now on' + geneCellTypeFile
		
		f = h5py.File(geneCellTypeFile,'r')
		g = f.__getitem__(uGene)
		#get gene expression indices
		indices = []
		idx = []
		indices = [i for i, s in enumerate(attrs) if 'Gene expression' in s]
		if len(gexMat) == 0:
			gexMat = np.mat(g[:,0]).T
			geMatName = []
		else:
			for idx in indices:
				gexMat = np.hstack((gexMat, np.mat(g[:,idx]).T))
				print gexMat.shape
				geMatName.append(f.attrs.values()[idx])
		f.close()
		print geMatName
		print 'oneCellTypeDone'
		







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