'''
Script that plots the matching cell type's epigenetic sums
'''
import prettyplotlib as ppl
import numpy as np
import scipy.spatial.distance as spsd
import matplotlib.pyplot as plt
import matplotlib as mpl
from prettyplotlib import brewer2mpl
from sys import argv
import time
import os
import sys

script, baseCellType, matchCellType, histoneMark = argv

#testing code
#baseCellType = 'hela'
#matchCellType = 'gm12878'
#histoneMark = 'h3k4me3'

#Some sanity checks to make sure the inputs are right
if baseCellType not in ['gm12878','hela','hepg2','huvec', 'k562']:
    print 'Matching cell type not in k562, gm12878, hela, hepg2 or huvec'
    sys.exit()

if matchCellType not in ['gm12878','hela','hepg2','huvec', 'k562']:
    print 'Matching cell type not in k562, gm12878, hela, hepg2 or huvec'
    sys.exit()

#### Setting paths
dataFolder = '/cbio/grlab/home/dkuo/projects/epican/results/epiSums/'
filepaths = [
'/cbio/grlab/home/dkuo/projects/epican/results/240214/k562Rep1_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/k562Rep2_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/gm12878Rep1_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/gm12878Rep2_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/helaRep1_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/helaRep2_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/hepg2Rep1_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/hepg2Rep2_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/huvecRep1_neither_EucDistAll.txt',
'/cbio/grlab/home/dkuo/projects/epican/results/240214/huvecRep2_neither_EucDistAll.txt'
]

#set result folder
date = time.strftime("%d%m%y/")
resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/' + date
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)

#Return a matrix of values from the histoneMarkList
def concatStrandExpSum(histoneMarkList):
    x = []
    if len(histoneMarkList) < 2:
        fp = dataFolder + str(histoneMarkList[0])
        tempCol = []
        for line in open(fp):
            if 'Sum' not in line:
                fileCol = line.split('\t')
                if len(fileCol) >= 2:
                    tempCol.append(float(fileCol[1]))
        return tempCol
    else:
        for hmp in histoneMarkList:
            fp = dataFolder + hmp
            tempCol = []
            for line in open(fp):
                if 'Sum' not in line:
                    fileCol = line.split('\t')
                    if len(fileCol) >= 2:
                        tempCol.append(float(fileCol[1]))
            if len(x) == 0:
                x = tempCol
            else:
                x = np.vstack((x,tempCol))
        print 
        return x.T

#plotting the two vectors
def plotXyHexbins(x, y, rep):
    #Base is X, Match is Y
    x = np.array(x)
    y = np.array(y)
    x = np.log(x+1)
    y = np.log(y+1)
    corr = np.corrcoef(x,y)[0][1]
    corr = '%s' % float('%.4g' % corr)
    plt.subplot(111)
    plt.hexbin(x, y, gridsize=75, bins='log', cmap=plt.cm.Spectral_r, mincnt = 1)
    plt.title("%s%s and %s; %s, %s Correlation for %s Genes" % \
        (baseCellType,rep, matchCellType, str(histoneMark.upper()), str(corr), \
        str(len(matchArr))))
    plt.xlabel("%s mean sum, gene length adj (log)" % (baseCellType) )
    plt.ylabel("%s mean sum, gene length adj (log)" % (matchCellType) )
    cb = plt.colorbar()
    cb.set_label('log10(N)')
    fig = plt.gcf()
    fig.set_size_inches(11,8.5)
    nmParam = '%s%s_%s' % (baseCellType, 
        rep, matchCellType)
    nmParam = nmParam + '_' + histoneMark.upper()
    pltNm = resultsFolder + nmParam 
    plt.savefig(pltNm, dpi=200)
    plt.clf()
    plt.close()
    metaDataNm = pltNm + '.txt'
    with open(metaDataNm, 'wb') as f:
        f.write('The base cell type was %s%s \n' % (baseCellType,rep))
        f.write('The match cell type was %s \n' % (matchCellType))
        f.write('The HM considered was %s \n' % (histoneMark))
        f.write('Total Number of Genes was %s \n' % (str(len(matchArr))))
        f.write("--------------------------------------\n")
        f.write('Base cell type files considered were: \n')
        f.write(str(baseHml))
        f.write('\n')
        f.write("--------------------------------------\n")
        f.write('Match cell type files considered were: \n')
        f.write(str(matchHml))

#GeneWidth
def retGenWidths():
    widthFiles = ['/cbio/grlab/home/dkuo/projects/epican/results/geSums/geneWidth_minus.txt', \
    '/cbio/grlab/home/dkuo/projects/epican/results/geSums/geneWidth_plus.txt']
    #Get gene widths
    minWidth = np.loadtxt(widthFiles[0], usecols = (1,))
    pluWidth = np.loadtxt(widthFiles[1], usecols = (1,))
    #Add 2kb because of promoter
    minWidth += 2000
    pluWidth += 2000
    return(pluWidth,minWidth)

geneWidthPlu,geneWidthMin = retGenWidths()
geneWidths = np.concatenate((geneWidthPlu,geneWidthMin))

#Get file names
listing = os.listdir(dataFolder)

#Getting filepaths
fpList = []
for fp in filepaths:
    if baseCellType.upper() in fp.upper():
        fpList.append(fp)

#Main function
for oneFp in fpList:
    if 'Rep1' in oneFp:
        rep = 'Rep1'
    else:
        rep = 'Rep2'
    inputRaw = np.loadtxt(oneFp, delimiter='\t',
    dtype={'names':('geneName', 'cellType'),'formats':('S20','S20')})
    #get indices
    matchArr = []
    for i in range(len(inputRaw)):
        if matchCellType.upper() in inputRaw[i][1].upper():
            matchArr.append(i)
    #parse file names for histone mark name
    histoneMarkList = []
    for l in listing:
        if histoneMark.upper() in l.upper():
            histoneMarkList.append(l)
    baseHml = []
    matchHml = []
    for hml in histoneMarkList:
        if baseCellType.upper() in hml.upper():
            baseHml.append(hml)
        elif matchCellType.upper() in hml.upper():
            matchHml.append(hml)   
    if len(baseHml) == 0 or len(matchHml) == 0:
        print 'Insufficient samples for %s' % (histoneMark)
        sys.exit()
    baseMat = np.array(concatStrandExpSum(baseHml))
    matchMat = np.array(concatStrandExpSum(matchHml))
    basePts = []
    matchPts = []
    if len(baseHml) == 1:
        for i in matchArr:
            basePts.append(baseMat[i]/geneWidths[i])
    elif len(baseHml) != 1:
        for i in matchArr:
            basePts.append(np.mean(baseMat[i,:])/geneWidths[i])
    if len(matchHml) == 1:
        for i in matchArr:
            matchPts.append(matchMat[i]/geneWidths[i])
    elif len(matchHml) != 1:
        for i in matchArr:
            matchPts.append(np.mean(matchMat[i,:])/geneWidths[i])
    plotXyHexbins(basePts,matchPts,rep)