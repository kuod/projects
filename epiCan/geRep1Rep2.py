#Input
#   hdf5 file per gene
#
#Output
#   plot showing gene expression against histone modification
import h5py
import os
import numpy as np
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
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/131208/'
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)

#cellTypeList = ['k562','gm12878','h1hesc','hela','hepg2', 'huvec']
cellTypeList = ['k562','gm12878','hela','hepg2', 'huvec']

start = time.time()
#one cell type at a time
for cellType in cellTypeList:
    cellTypeRep1Vec = []
    cellTypeRep2Vec = []
    #one gene at a time for celltype
    for geneFile in geneFiles:
        if geneFile.endswith(cellType + ".hdf5"):
            try:
                f = h5py.File(geneFile,'r')
                geneName = geneFile.split('/')[-1].split('_')[0]
                g = f.__getitem__(geneName)
                #NEED TO FIX!
                geIdx = [i for i, s in enumerate(f.attrs.values()) if 'GE' in s]
                if len(geIdx) == 0:
                    continue
                cellTypeRep1Vec.append(np.sum(g[:,geIdx[0]]))
                cellTypeRep2Vec.append(np.sum(g[:,geIdx[1]]))
                f.close()
            except(IOError):
                print str(geneFile) + " doesn't work"
                continue
    cellTypeCorr = np.corrcoef(cellTypeRep1Vec,cellTypeRep2Vec)[0][1]
    cellTypeCorr = '%s' % float('%.4g' % cellTypeCorr)
    fig, ax = plt.subplots()
    ax.scatter(cellTypeRep1Vec, cellTypeRep2Vec)
    ax.set_xlabel('Gene Expression Rep 1 (fpkm)')
    ax.set_ylabel('Gene Expression Rep 2 (fpkm)')
    ax.set_title(cellType + ': Correlation of ' + str(cellTypeCorr))
    figName = resultsFolder + cellType + 'geRep1Rep2.png'
    fig.savefig(figName, bbox_inches='tight',dpi = 200)
    plt.close()
    plt.clf()
    end = time.time()
    print 'one cell: ' + str(end-start)
end = time.time()
print 'Total time: ' + str(end - start)

