import h5py
import os
import numpy as np
from scipy import spatial as spsp
import scipy as sp
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
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/131209/'
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)

#cellTypeList = ['k562','gm12878','h1hesc','hela','hepg2', 'huvec']
cellTypeList = ['k562','gm12878','hela','hepg2', 'huvec']

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
                geIdx = [i for i, s in enumerate(f.attrs.values()) if 'GE-' in s]
                hmIdx = [i for i, s in enumerate(f.attrs.values()) if 'H3k4me3' in s]
                if len(geIdx) == 0 or len(hmIdx) == 0:
                    continue
                cellTypeGeVec.append(np.sum(g[:,geIdx[0]]))
                cellTypeHmVec.append(np.sum(g[:,hmIdx[0]]))
                f.close()
            except(IOError):
                print str(geneFile) + " doesn't work"
                continue
    cellTypeGeVec = sp.log(cellTypeGeVec)
    cellTypeGeVec[np.isneginf(cellTypeGeVec)] = 0
    cellTypeHmVec = sp.log(cellTypeHmVec)
    cellTypeHmVec[np.isneginf(cellTypeHmVec)] = 0
    print 'geVec shape is: ' + str(len(cellTypeGeVec))
    cellTypeCorr = np.corrcoef(cellTypeGeVec,cellTypeHmVec)[0][1]
    cellTypeCorr = '%s' % float('%.4g' % cellTypeCorr)
    fig, ax = plt.subplots()
    ax.scatter(sp.log(cellTypeGeVec), sp.log(cellTypeHmVec))
    ax.set_xlabel('Gene Exp across gene (log fpkm)')
    ax.set_ylabel('H3k4me3 across gene (log fpkm)')
    #ax.set_ylim(0,2E6)
    ax.set_title(cellType + ': Correlation of ' + str(cellTypeCorr))
    figName = resultsFolder + cellType + '_expGeGene_HmGene_corr.png'
    fig.savefig(figName, bbox_inches='tight',dpi = 200)
    plt.close()
    plt.clf()
    print 'saved file ' + figName
    end = time.time()
    print 'one celltype takes: ' + str(end-start)
end = time.time()
print 'Total time: ' + str(end - start)
