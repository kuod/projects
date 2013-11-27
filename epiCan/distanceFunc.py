#Input
#	training matrix of values including ge
#	test vector of comparable gene expression
#Output
#	predicted histone track based on euclidian distance
import h5py
import os
import numpy as np
from scipy import spatial as spsp
import matplotlib.pyplot as plt
import time

#All gene files in geneFiles
geneFolder = '/cbio/grlab/share/databases/encode/hdf5/'
geneFiles = []
os.chdir( geneFolder )
geDir = os.getcwd() 
for geneFile in os.listdir("."):
	if geneFile.endswith(".hdf5"):
		geneFiles.append(geDir + '/' + geneFile)
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/131128/'
if not os.path.exists(resultsFolder):
	os.makedirs(resultsFolder)

def dist(inputMat, mode):
	if mode == 'euc':
		Y = spsp.distance.pdist(inputMat, 'euclidean')
		return spsp.distance.squareform(Y)
		
#getsGeneNames
uniqGenes = []
for gF in geneFiles:
	uniqGenes.append(gF.split('/')[-1].split('_')[0])
uniqGenes = np.unique(uniqGenes)

subsetUniqGenes = uniqGenes[2:10]

start = time.time()
for uGene in subsetUniqGenes:
	#find the files that correspond to that gene
	geneCellTypesList = filter(lambda x: uGene in x, geneFiles)
	gexMat = []
	gexMat = np.mat(gexMat)
	geMatName = []
	for geneCellTypeOne in geneCellTypesList:
		print 'now on ' + geneCellTypeOne
		f = h5py.File(geneCellTypeOne, 'r')
		g = f.__getitem__(uGene)
		print 'so far so good on ' + geneCellTypeOne
		#get gene expression indices
		indices = [i for i, s in enumerate(f.attrs.values()) if 'Gene expression' in s]
		print indices
		if len(gexMat) == 0:
			gexMat = np.mat(g[:,0]).T
			geMatName = []
		else:
			for idx in indices:
				geMatName.append(f.attrs.values()[idx])
				print 'gexmat shape is ' + gexMat.shape
				print 'index is ' + str(idx)
				gexMat = np.hstack((gexMat, np.mat(g[:,idx]).T))
				
		f.close()
	print gexMat.shape
	gexMat2 = gexMat[:,1:]
	distMat = dist(gexMat2,'euc')
	print "dismat shape is " + str(distMat.shape)
	print geMatName
	fig = plt.figure()
	fig.suptitle(uGene,fontsize=14, fontweight='bold')
	ax = fig.add_subplot(111)
	cax = ax.matshow(distMat, interpolation='nearest')
	#ax.set_xticklabels(['']+geMatName, rotation=45)
	ax.set_yticklabels(['']+geMatName)
	figName = resultsFolder + uGene + '_euc.png'
	fig.savefig(figName)
	print 'onto the next one'
	plt.close()
	plt.clf()
	
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