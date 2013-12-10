#generates h1hesc hdf5 files
import numpy as np
import sys, os, csv
import os.path
import time
import h5py
import re
from bx.bbi.bigwig_file import BigWigFile

#set test file paths
bwTestFile = "/cbio/grlab/share/databases/encode/gm12878/rnaSeqPolACyto/wgEncodeCshlLongRnaSeqGm12878CytosolPapMinusRawSigRep1.bigWig"
geneListFile = "/cbio/grlab/share/databases/genomes/H_sapiens/GENCODE/release18/gencode.v18.annotation.genesPcLncrna.gtf"
methylFile = "/cbio/grlab/share/databases/encode/k562/methyl/wgEncodeHaibMethylRrbsK562HaibSitesRep1_chr21.bed"

#TODO Take in path of HM folder  
hm_1 = "/cbio/grlab/share/databases/encode/k562/histonePeaks/wgEncodeBroadHistoneK562H3k4me3StdAln_chr21.bed"

#Load gene information in GTF format
genesRaw = np.loadtxt(geneListFile, delimiter='\t', 
    dtype={'names':('chromosome', 'annotation','element', 'start', 'end', 
    'score', 'strand', 'frame', 'attribute'),'formats': ('S5', 'S10','S10',
    'u4','u4','S1','S1','S1','S100')})

chromosome = 'chr21'
outputDir = '/cbio/grlab/share/databases/encode/hdf5/'
#Get methylation patterns
#methylation in bed
#methylRaw = np.loadtxt(methylFile, delimiter='\t', usecols=(0,1,2,4), dtype={'names':('chrom', 'chromStart', 'chromEnd',
#   'score'), 'formats': ('S5', 'u4','u4','u4')}) 
#methylRaw = np.loadtxt(methylFile, delimiter='\t', usecols=[1,4])
#methylMed = np.median()
#methylSites = []
#subset based on score
#for m in methylRaw:
#   if m[3] >  

#Parse start stop and gene names
starts = []
stops = []
genes = []
for g in genesRaw:
    starts.append(g[3]-500)
    stops.append(g[4]+500)
    curGene = g[8].split()[1].translate(None, '";')
    genes.append(curGene)
#cellTypeList = ['k562','gm12878','h1hesc','hela','hepg2', 'huvec']
#cellTypeList = ['k562','gm12878','h1hesc','hela']
cellTypeList = ['h1hesc']
#loadBW func
def loadBw(file):
    bw = BigWigFile(open(file))
    return bw

#generate gene by gene matrix
def saveCellTypeDat(genes, starts, stops, cellType, geFileList, hmFileList):
    start = time.time()
    for i in range(len(genes)):
        middle1 = time.time()
        fileName = outputDir + genes[i] + '_' + cellType +'.hdf5'
        
        if os.path.isfile(fileName):
            try:
                h = h5py.File(fileName)
                print "The file exists: " + str(fileName)
                os.remove(fileName)
                h.close()
                continue
            except(IOError):
                print str(fileName) + 'doesnot open so deleting. Please rerun'
                os.remove(fileName)
        else:
            #pos index
            posMat = []
            posMat = np.mat(np.arange(starts[i], stops[i], 1))
            #geneExp summary
            #ge from bigwiglist
            geMat = []
            hmMat = []
            geTrackNames = []
            hmTrackNames = []
            for ge_count in range(len(geFileList)):
                if ge_count == 0:
                    geFileName = geFileList[ge_count].split("/")[-1]
                    geFileName = re.sub('wgEncodeCshlLongRnaSeq','',geFileName)
                    geFileName = re.sub('.bigWig','',geFileName)
                    geTrackNames.append(geFileName)
                    geFile = loadBw(geFileList[ge_count])
                    oneGeMat = np.mat(geFile.get_as_array(chromosome, starts[i], stops[i]).T)
                    oneGeMat[np.isnan(oneGeMat)] = 0
                    geMat = np.mat(geMat)
                    geMat = np.mat(oneGeMat.T)
                else:
                    geFileName = geFileList[ge_count].split("/")[-1]
                    geFileName = re.sub('wgEncodeCshlLongRnaSeq','',geFileName)
                    geFileName = re.sub('.bigWig','',geFileName)
                    geTrackNames.append(geFileName)
                    geFile = loadBw(geFileList[ge_count])
                    oneGeMat = np.mat(geFile.get_as_array(chromosome, starts[i], stops[i]).T)
                    oneGeMat[np.isnan(oneGeMat)] = 0
                    geMat = np.hstack((geMat, oneGeMat.T))
            for hm_count in range(len(hmFileList)):
                if hm_count == 0:
                    hmTrackNames.append(hmFileList[hm_count].split("/")[-1])
                    hmFile = loadBw(hmFileList[hm_count])
                    oneHmMat = np.mat(hmFile.get_as_array(chromosome, starts[i], stops[i]).T)
                    oneHmMat[np.isnan(oneHmMat)] = 0
                    hmMat = np.mat(hmMat)
                    hmMat = np.mat(oneHmMat.T)
                else:
                    hmTrackNames.append(hmFileList[hm_count].split("/")[-1])
                    hmFile = loadBw(hmFileList[hm_count])
                    oneHmMat = np.mat(hmFile.get_as_array(chromosome, starts[i], stops[i]).T)
                    oneHmMat[np.isnan(oneHmMat)] = 0
                    hmMat = np.hstack((hmMat, oneHmMat.T))
            
            temp = np.hstack((posMat.T,geMat,hmMat))
            #print temp.shape
            f = h5py.File(fileName, 'a')
            f.create_dataset(name=genes[i], data=temp)
            for colNum in range(temp.shape[1]):
                if colNum == 0:
                    f.attrs.create(name=str(colNum), data='position')
                elif 0 < colNum < (len(geFileList) - 1):
                #replace len(geFileList) with number of GE Tracks
                    f.attrs.create(name=str(colNum), data='GE-' + geTrackNames[colNum-1])
                    #TODO: include methylation tracks as well
                else:
                    f.attrs.create(name=str(colNum), data='HM-' + hmTrackNames[colNum - len(geFileList) -1])
            f.close()
            middle2 = time.time()
#            print 'One gene took: ' + str(middle2 - middle1)
        #np.savez(fileName, temp)
    end = time.time()
    print 'Total time for genes in a cell line: ' + str(end - start)

#Subset parameters
startSubset = starts[0:1000]
stopSubset = stops[0:1000]
geneSubset = genes[0:1000]

for ct in cellTypeList:
    geFileList = []
    gePath = "/cbio/grlab/share/databases/encode/" + ct + "/rnaSeqPolACyto/"
    os.chdir( gePath )
    geDir = os.getcwd() 
    for files in os.listdir("."):
        if files.endswith(".bigWig"):
            geFileList.append(geDir + '/' + files)
    hmFileList = []
    os.chdir("/cbio/grlab/share/databases/encode/" + ct + "/histone/bw/")
    hmDir = os.getcwd()
    for files in os.listdir("."):
        if files.endswith(".bw"):
            hmFileList.append(hmDir + '/' + files)
    saveCellTypeDat(geneSubset, startSubset, stopSubset, ct, geFileList, hmFileList)

#for i in range(10):
#   #pos index
#   posMat = []
#   posMat = np.mat(np.arange(startSubset[i], stopSubset[i], 1))
#   #geneExp summary
#
#   #ge from bigwig
#   geMat = []
#   geMat = np.mat(bwf.get_as_array(chromosome, startSubset[i], stopSubset[i]).T)
#   geMat[np.isnan(geMat)] = 0
#
#   ts1 = time.time()
#   hmTrack = htGen(startSubset[i], stopSubset[i], hm1Coord)
#   hmTrack = np.mat(hmTrack) 
#   ts2 = time.time()
#   print ts2 - ts1
#
#   temp = np.hstack((posMat.T,geMat.T, hmTrack.T))
#   #fileName = '/cbio/grlab/projects/epiCan/results/chr21/' + genes[i] 
#   #numpy.savez(fileName, temp)
#   print np.max(temp[:,2])
#   print '\n'
#   print np.max(temp[:,1])
#########################################

#def loadHm(file):
#   hmArray = np.loadtxt(file, delimiter='\t', usecols=(0,1,2,4), 
#   dtype={'names':('chrom', 'chromStart', 'chromEnd','score'), 
#   'formats': ('S5', 'u4','u4','u4')}) 
#   return hmArray

###older test code
#load ge bw file 
#geFileList = []
#os.chdir("/cbio/grlab/share/databases/encode/k562/rnaSeqPolACyto/")
#geDir = os.getcwd() 
#for files in os.listdir("."):
#    if files.endswith(".bigWig"):
#        geFileList.append(geDir + '/' + files)
#
#hmFileList = []
#os.chdir("/cbio/grlab/share/databases/encode/k562/histone/bw/")
#hmDir = os.getcwd()
#for files in os.listdir("."):
#   if files.endswith(".bw"):
#       hmFileList.append(hmDir + '/' + files)#
