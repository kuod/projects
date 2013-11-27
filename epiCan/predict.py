#Input
#	training matrix of values including ge
#	test vector of comparable gene expression
#Output
#	predicted histone track based on euclidian distance

import numpy as np
from numpy import linalg as la
from scipy import spatial as spsp
import h5py 

#	Mode 1: directory with gene files 
#	Mode 2: gene file name
#	Mode 3: gene expression track
geneFolder = '/cbio/grlab/home/dkuo/temp'
geneFiles = []
os.chdir( geneFolder )
geDir = os.getcwd() 
for files in os.listdir("."):
	if files.endswith(".hdf5"):
		geneFiles.append(geDir + '/' + files)

#getUniqGenes
uniqGenes = []
for f in geneFiles:
	uniqGenes.append(f.split('/')[-1].split('.')[0])

uniqGenes = np.unique(uniqGenes)

for uGene in uniqGenes:
	if uGene in geneFiles:
		print geneFiles

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