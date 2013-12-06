import h5py
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
import time

eucData = '/cbio/grlab/home/dkuo/projects/epican/results/eucData.hdf5'
eucDataH5 = h5py.File(eucData, 'r')

disMatList = eucDataH5.__getitem__('disMats')
geneList = eucDataH5.__getitem__('genelist')
nameList = eucDataH5.__getitem__('namelist')

#which sample to care about
testSamp = 0

#slice up the distancematrix
disMatSlice = disMatList[:,testSamp,:]

resultsFolder = '/cbio/grlab/home/dkuo/projects/epican/results/131204/'
if not os.path.exists(resultsFolder):
	os.makedirs(resultsFolder)

for i in range(disMatSlice.shape[0]):
	sort_index =[]
	curSlice = []
	curSlice = disMatSlice[i]
	sort_index = numpy.argsort(curSlice)
	print namelist[i][sort_index][0]
	print 'Nearest neighbor: ' + str(namelist[i][sort_index[1]])
	if np.sum(curSlice) == 0:
		print 'skipped'
		pass
	else:
		print sorted(range(len(curSlice)), key=curSlice.__getitem__)

fig = plt.gcf()
ax = fig.add_subplot(111)
fig.set_size_inches(15,5)
ax.set_frame_on(False)
plt.tick_params(axis='x', which='both',bottom='off',top='off',labelbottom='off')
for i in range(disMatList.shape[0]):
	start = time.time()
	fig.suptitle(geneList[i],fontsize=14, fontweight='bold')
	cax = ax.matshow(disMatList[i,:,:], interpolation='nearest', cmap = cm.Greys_r,rasterized=True)
	cbar = fig.colorbar(cax)
	ax.set_yticks(np.arange(disMatList[i,:,:].shape[0]+0.5),minor=False)
	ax.set_yticklabels(nameList[i], minor=False)
	for tick in ax.yaxis.get_major_ticks():
		tick.tick1On = False
		tick.tick2On = False
	for tick in ax.xaxis.get_major_ticks():
		tick.tick1On = False
		tick.tick2On = False
	figName = resultsFolder + geneList[i] + '_euc.png'
	fig.savefig(figName, bbox_inches='tight',dpi = 200)
	end = time.time()
	print 'Total time: ' + str(end - start)
	plt.clf()
	raw_input()
plt.close()
plt.clf()