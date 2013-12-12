import h5py
import os
import numpy as np
from scipy import spatial as spsp
import scipy as sp
import matplotlib
import time
import sys

if len(sys.argv) != 2:
    print 'Histone Mark required.'
    print 'Please specify one of the following: H3k9ac, H4k20me1, H3k27me3, H3k4me3, Control, Ctcf, H3k36me3, H3k79me2, H3k4me2, Pol2b, InputStd, H3k27ac, or H3k4me3.'
    sys.exit()

hm_list=['H3k9ac','H4k20me1','H3k27me3','H3k4me3','Control','Ctcf','H3k36me3','H3k36me3','H3k79me2','H3k4me2','Pol2b','InputStd','H3k27me3','H3k27ac','H3k4me3']

hMark = str(sys.argv[1])

if hMark in hm_list:
    print 'Plotting ' + hMark 
else:
    print 'Histone Mark not found in list. Exiting'
    sys.exit()

date = time.strftime("%d%m%y/")

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
geneFolder = '/cbio/grlab/share/databases/encode/hdf5_2kb/'
geneFiles = []
os.chdir( geneFolder )
geDir = os.getcwd()
for geneFile in os.listdir("."):
    if geneFile.endswith(".hdf5"):
        geneFiles.append(geDir + '/' + geneFile)
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/' + date
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
                hmIdx = [i for i, s in enumerate(f.attrs.values()) if hMark in s]
                if len(geIdx) == 0 or len(hmIdx) == 0:
                    continue
                cellTypeGeVec.append(np.sum(g[2000:,geIdx[0]])+1)
                cellTypeHmVec.append(np.sum(g[:,hmIdx[0]])+1)
                f.close()
            except(IOError):
                print str(geneFile) + " doesn't work"
                continue
    cellTypeGeVec = sp.log(cellTypeGeVec)
    cellTypeGeVec[np.isneginf(cellTypeGeVec)] = 0
    cellTypeHmVec = sp.log(cellTypeHmVec)
    cellTypeHmVec[np.isneginf(cellTypeHmVec)] = 0
    if len(cellTypeGeVec) == 0:
        print 'GeVec length is zero, skipping ' + str(cellType)
        continue
    if len(cellTypeHmVec) == 0:
        print 'HmVec length is zero, skipping ' + str(cellType)
        continue
    print 'geVec shape is: ' + str(len(cellTypeGeVec))
    cellTypeCorr = np.corrcoef(cellTypeGeVec,cellTypeHmVec)[0][1]
    cellTypeCorr = '%s' % float('%.4g' % cellTypeCorr)
    fig, ax = plt.subplots()
    ax.scatter(cellTypeGeVec, cellTypeHmVec)
    ax.set_xlabel('RNA Seq across gene (log fpkm)')
    ax.set_ylabel(hMark +' across gene (log fpkm)')
    #ax.set_ylim(0,2E6)
    ax.set_title(cellType + ': Correlation of ' + str(cellTypeCorr))
    figName = resultsFolder + cellType + '_expGeGene' + hMark + '.png'
    fig.savefig(figName, bbox_inches='tight',dpi = 200)
    plt.close()
    plt.clf()
    print 'saved file ' + figName
    end = time.time()
    print 'one celltype takes: ' + str(end-start)
end = time.time()
print 'Total time: ' + str(end - start)
