import matplotlib.pyplot as plt
from matplotlib import cm



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