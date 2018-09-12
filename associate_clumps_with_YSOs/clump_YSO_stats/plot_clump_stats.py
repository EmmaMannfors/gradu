#######################################################################
#		Plot histograms of clump statistics
#######################################################################
#			Jul 31, 2018
#######################################################################

import os
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import norm
import matplotlib.mlab as mlab

##################################
# Paths
##################################
file_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/clumps_vs_all/clump_association/total_clumps_vs_YSO_clumps.txt'

thesis = '/home/emma/gradu/Thesis file/'


##################################
# Code
##################################
ratioList = []
ysoList = []
starLessList = []
totList = []

with open(file_in) as f: 
	next(f)
	for line in f: 
		if line[0] == '#':
			break
		p = line.split()
		field = p[0]
		totClumps = int(p[1])
		nClumps = int(p[2])
		ratio = 100*float(p[3])
		starless = totClumps - nClumps

		if starless > 60: 
			print field, starless,nClumps
		if nClumps > 60: 
			
			print field, starless,nClumps
		

		starLessList.append(starless)
		ysoList.append(nClumps)
		totList.append(totClumps)
		ratioList.append(ratio)

##############
# YSO/all
##############
sigma = np.std(ratioList)
mu = np.mean(ratioList)
nbins = 20
n, bins, patches = plt.hist(ratioList,nbins,color='teal',alpha=0.7,normed=1)
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'k--', linewidth=2)
plt.title('protostellar / all clumps')
plt.xlabel('Cores which contain YSOs [%]')
plt.ylabel('Frequency')

s = str(round(sigma,3))
m = str(round(mu,3))
plt.text(70,0.025,'$\sigma$ = '+s+'%\n'+'$\mu$ = '+m+'%')

plt.show()
#plt.savefig(thesis+'stellar_all_cores_ratio.png')


###############
# YSO/starless
###############
sigma = np.std(ysoList)
mu = np.mean(ysoList)
nbins = 25

sigmaS = np.std(starLessList)
muS =  np.mean(starLessList)

n, bins, patches = plt.hist([ysoList,starLessList],nbins,label=['Protostellar clumps','Starless clumps'],alpha=0.6,range=(0,60),normed=1)

y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'b-', linewidth=2)

x = mlab.normpdf( bins, muS, sigmaS)
k = plt.plot(bins, x, 'g-', linewidth=2)

plt.title('Clumps per field (n$_{clumps}$ < 60)')
plt.xlabel('Number of clumps in a field')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
s = str(round(sigma,3))
m = str(round(mu,3))
plt.text(45,0.09,'$\sigma_{YSO}$ = '+s+'\n'+'$\mu_{YSO}$ = '+m)

sS = str(round(sigmaS,3))
mS = str(round(muS,3))
plt.text(45,0.075,'$\sigma_{starless}$ = '+sS+'\n'+'$\mu_{starless}$ = '+mS)
plt.show()
#plt.savefig(thesis+'stellar_starless_cores.png')




	
#	y = mlab.normpdf( bins, mu, sigma)
#
#	(mu, sigma) = norm.fit(clumpStarList)
#	f = mlab.normpdf( bins, mu, sigma)
#	
#	l = plt.plot(bins, y, 'b--', linewidth=2,label='best fit of all YSOs')
#	plt.plot(bins, f, color='darkgreen', linewidth=2,label='best fit of clump YSOs'),#'g-', linewidth=2,label='distribution of clump stars')











#(mu, sigma) = norm.fit(ratioList)
#	n, bins, patches = plt.hist(ratioList, 40, facecolor='midnightblue', alpha=0.6) #normed=1, 
#	y = mlab.normpdf( bins, mu, sigma)
#	l = plt.plot(bins, y, 'k--', linewidth=2)#,label='Distribution function ($\sigma$)')# = %s)', %s(sigma))


