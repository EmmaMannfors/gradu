###############################################################################
#		 How many stars each individual clump is associated with
###############################################################################
#				jul 30, 2018
###############################################################################

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax

from scipy.stats import norm
import matplotlib.mlab as mlab

######################################
# Paths
######################################
files_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/clumps_vs_all/'

hist_out = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/clumps_vs_all/clump_association/histograms/'

text_out = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/clumps_vs_all/clump_association/'
######################################
# Functions
######################################
def str_to_int(val):
	if val == 'nan':
		num = 0
	else: 
		v = float(val)
		num = int(v)
	return num
######################################
# Init
######################################
clumps = dict()

######################################
# Code
######################################
with open(files_in+'ratio_of_stars_in_clumps.txt') as f: 
	for line in f: 
		p = line.split()


		if p[0] == 'Field:':
			field = p[1]
			#print p[1]
		elif p[0] == '#' or line[0] == '#':
			continue
		elif p[1] == 'Total stars in field':
			next(f)
			next(f)
		else: 
		#	if field != 'G093.51-04.31':
		#		continue


			field = field
			l = line.split('|')
			clump_str = l[0].strip()
			#if clump_str != 'nan':
				#cl = float(clump)
				#c = int(cl)
			try:
				# Gets rid of the total stars in field-variable, since that is an integer
				# and so can be converted straight from string to int
				clump = int(clump_str)
				continue
			except:
				clump = str_to_int(clump_str) # nan-values are set to zero
		
			occurrences = str(l[1].strip())
			
			#print clump,occurrences


			if field in clumps: # Not the first instance of this field 
				c = clumps[field]
				#print type(c)
				c.append(clump)
				c.append(occurrences)
				clumps[field] = c
			else: 
				clumps[field] = [clump,occurrences]

		

	#	if field != 'G093.51-04.31':
	#		continue
			
		#print line





#G026.54+00.72

allClumpsDic = dict()
meanList = []
maxList = []
minList = []
stdList = []

ysoDic = dict()

for key in clumps: 
	#print clumps[key]
	
	valList = []
	occList = []
	
	c = clumps[key]

	text = open(text_out+key+'.txt','w')
	text.write('# '+key+'\n# clump	occurrences\n')

	nClumps = 0
	
	for i in range(len(c)):
		if i%2 == 0: # even
			clump = int(c[i])
			occ = ''
		else: 
			occ = int(c[i])
			if clump != 0:
				occList.append(occ)
				nClumps += 1		# 1 more YSO-bearing clump (nan is not a clump)
		if occ != '':
			#if occ == 28:
			#	print clump,occ,key
			for i in range(occ):
				if key == 'G026.54+00.72':
					if clump != 0:
						valList.append(clump)
				else: 
					valList.append(clump)
			if clump == 0:
				text.write('  nan	'+str(occ)+'\n')
			else: 
				text.write('  '+str(clump)+'	'+str(occ)+'\n')

			if clump in allClumpsDic:
				o = allClumpsDic[clump]
				o += occ
				allClumpsDic[clump] = o
			else: 
				allClumpsDic[clump] = occ
			
			

			continue
		else: 
			continue

	#print key, nClumps
	ysoDic[key] = nClumps
	
	if len(occList) != 0:
		ma = np.max(occList)
		mi = np.min(occList)
		mean = np.mean(occList)
		std = np.std(occList)
	
		text.write(5*'#'+'\n')
		text.write('* max: '+str(ma)+'\n'+'* min: '+str(mi)+'\n')
		text.write('* mean: '+str(mean)+'\n'+'* std: '+str(std)+'\n')

		meanList.append(mean)
		maxList.append(ma)
		minList.append(mi)
		stdList.append(std)
	

#	nbins = np.max(valList)	
#	if nbins == 0:
#		continue
#	plt.hist(valList,nbins)
#	plt.title(key)

#	#plt.show()
#	plt.savefig(hist_out+key+'.png')		
#	plt.close()





allFile = open(text_out+'general_stats.txt','w')
allFile.write('# Overall stats of all star-associated clumps\n')
allFile.write('# Includes only clumps which have associated YSOs, but ignores "nan"-clumps.\n')

allFile.write('Max: \n')
ma = np.max(maxList)
mi = np.min(maxList)
mean = np.mean(maxList)
std = np.std(maxList)
allFile.write('max: '+str(ma)+'\n'+'min: '+str(mi)+'\n')
allFile.write('mean: '+str(mean)+'\n'+'std: '+str(std)+'\n')


allFile.write('Min: \n')
ma = np.max(minList)
mi = np.min(minList)
mean = np.mean(minList)
std = np.std(minList)
allFile.write('max: '+str(ma)+'\n'+'min: '+str(mi)+'\n')
allFile.write('mean: '+str(mean)+'\n'+'std: '+str(std)+'\n')


allFile.write('Mean: \n')
ma = np.max(meanList)
mi = np.min(meanList)
mean = np.mean(meanList)
std = np.std(meanList)
allFile.write('max: '+str(ma)+'\n'+'min: '+str(mi)+'\n')
allFile.write('mean: '+str(mean)+'\n'+'std: '+str(std)+'\n')


allFile.write('Stadard Deviation: \n')
ma = np.max(stdList)
mi = np.min(stdList)
mean = np.mean(stdList)
std = np.std(stdList)
allFile.write('max: '+str(ma)+'\n'+'min: '+str(mi)+'\n')
allFile.write('mean: '+str(mean)+'\n'+'std: '+str(std)+'\n')


A = open(text_out+'total_clumps_vs_YSO_clumps.txt','w')
A.write('field		total clumps	YSO clumps	ratio (yso/tot)\n')

ratioList = []
for key in ysoDic: 
	field = key
	n = ysoDic[key]
	#print key, n
#	ysoDic[key] = nClumps
	with open(files_in+'total_clumps.txt') as tot: 
		next(tot)
		for line in tot: 
			p = line.split()
			P = p[0].split('_')
			current = P[0]+'.'+P[1]+'.'+P[2]
			#print current
			if current == field:
				allC = float(p[1])
				ratio = str(n/allC)
				ratioList.append(float(ratio))
				#print field, allC,n,ratio
				A.write(field+'	'+p[1]+'		'+str(n)+'		'+ratio+'\n')




A.write(50*'#'+'\n')
A.write('# 50 hashtags\n')
A.write('# Stats of the ratios over all the fields\n')
A.write('Max: '+str(np.max(ratioList))+'\nMin: '+str(np.min(ratioList))+'\nMean: '+str(np.mean(ratioList))+'\nStd: '+str(np.std(ratioList))+'\n')
















