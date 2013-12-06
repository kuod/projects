#Input
#   hdf5 file per gene
#
#Output
#   plot showing gene expression against histone modification
import h5py
import os
import numpy as np
from scipy import spatial as spsp
import matplotlib
import time

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
geneFolder = '/cbio/grlab/share/databases/encode/hdf5/'
geneFiles = []
os.chdir( geneFolder )
geDir = os.getcwd()
for geneFile in os.listdir("."):
    if geneFile.endswith(".hdf5"):
        geneFiles.append(geDir + '/' + geneFile)
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/131206/'
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)

cellTypeList = ['k562','gm12878','h1hesc','hela','hepg2', 'huvec']


start = time.time()
#one cell type at a time
for cellType in cellTypeList:
    cellTypeGeVec = []
    cellTypeHmVec = []
    #one gene at a time for celltype
    for geneFile in geneFiles:
        if geneFile.endswith(cellType + ".hdf5"):
            try:
                f = h5py.File(geneFile,'r')
                geneName = geneFile.split('/')[-1].split('_')[0]
                g = f.__getitem__(geneName)
                #NEED TO FIX!
                geIdx = [i for i, s in enumerate(f.attrs.values()) if 'Gene expression' in s]
                hmIdx = [i for i, s in enumerate(f.attrs.values()) if 'H3k4me3' in s]
                if len(geIdx) == 0 or len(hmIdx) == 0:
                    continue
                cellTypeGeVec.append(np.sum(g[0:1000,geIdx[0]]))
                cellTypeHmVec.append(np.sum(g[0:1000,hmIdx[0]]))
            except(IOError):
                print str(geneFile) + " doesn't work"
                continue
    fig, ax = plt.subplots()
    ax.scatter(cellTypeGeVec, cellTypeHmVec)
    ax.set_xlabel('Gene Expression')
    ax.set_ylabel('h3k4 Methylation')
    ax.set_ylim(0,2E6)
    ax.set_title(cellType)
    figName = resultsFolder + cellType + '_chr1promoter-zoom.png'
    fig.savefig(figName, bbox_inches='tight',dpi = 200)
    plt.close()
    plt.clf()
end = time.time()
print 'Total time: ' + str(end - start)