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
from matplotlib import cm
import time

#All gene files in geneFiles
geneFolder = '/cbio/grlab/share/databases/encode/hdf5/'
geneFiles = []
os.chdir( geneFolder )
geDir = os.getcwd() 
for geneFile in os.listdir("."):
	if geneFile.endswith(".hdf5"):
		geneFiles.append(geDir + '/' + geneFile)
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/131203/'
if not os.path.exists(resultsFolder):
	os.makedirs(resultsFolder)

def dist(inputMat, mode):
	if mode == 'euc':
		Y = spsp.distance.pdist(inputMat, 'euclidean')
		return spsp.distance.squareform(Y)
	if mode == '':
		pass

quitInd = False

#getsGeneNames
uniqGenes = []
for gF in geneFiles:
	uniqGenes.append(gF.split('/')[-1].split('_')[0])
uniqGenes = np.unique(uniqGenes)

disMatList = []
geMatNameList = []

for uGene in uniqGenes:
	start = time.time()
	#find the files that correspond to that gene
	#gets only the unique genes
	uGeneCellTyps = filter(lambda x: uGene in x, geneFiles)
	#process one gene at a time
	gexMat = []
	gexMat = np.mat(gexMat)
	geMatName = []
	for uGeneCellTyp in uGeneCellTyps:
		#print 'now on ' + uGeneCellTyp
		#try to open the file
		try: 
			f = h5py.File(uGeneCellTyp, 'r')
		#if file corrupted, go to next 
		except (IOError):
			continue
		g = f.__getitem__(uGene)
		#print 'so far so good on ' + uGeneCellTyp
		#get gene expression indices
		indices = [i for i, s in enumerate(f.attrs.values()) if 'Gene expression' in s]
		if len(indices) == 0:
			#f it, no gene expression, no prediction
			continue
		#print indices
		#print 'indices is' + str(indices)
		#check columns
		if gexMat.shape[1] == 0:
			#print 'shape of first col is' + str(g[:,0].shape)
			gexMat = np.mat(g[:,0]).T
			#print 'gexmat shape is ' + str(gexMat.shape)
			geMatName = []
		#iterate through indices
		for idx in indices:
			geMatName.append(f.attrs.values()[idx])
		#	print 'gexmat shape is ' + str(gexMat.shape)
		#	print 'index is ' + str(idx)
		#	print 'shape of idx is' + str(g[:,idx].shape)
			gexMat = np.hstack((gexMat, np.mat(g[:,idx]).T))
		#	print "here's what gexMat is" + str(gexMat.shape)
		f.close()
		#print "closed file"
	gexMat2 = gexMat[:,1:]
	distMat = dist(gexMat2.T,'euc')
	#print "dismat shape is " + str(distMat.shape)
	#print geMatName
	disMatList.append(distMat)
	geMatNameList.append(geMatName)


	fig = plt.gcf()
	fig.suptitle(uGene,fontsize=14, fontweight='bold')
	ax = fig.add_subplot(111)
	cax = ax.matshow(distMat, interpolation='nearest', cmap = cm.Greys_r,rasterized=True)
	cbar = fig.colorbar(cax)
	fig.set_size_inches(15,5)
	ax.set_frame_on(False)
	ax.set_yticks(np.arange(distMat.shape[0]+0.5),minor=False)
	ax.set_yticklabels(geMatName, minor=False)
	for tick in ax.yaxis.get_major_ticks():
		tick.tick1On = False
		tick.tick2On = False
	plt.tick_params(axis='x', which='both',bottom='off',top='off',labelbottom='off')
	figName = resultsFolder + uGene + '_euc.png'
	fig.savefig(figName, bbox_inches='tight',dpi = 200)
	plt.close()
	plt.clf()
	end = time.time()
	print 'Total time: ' + str(end - start)
	if quitInd:
		break
	