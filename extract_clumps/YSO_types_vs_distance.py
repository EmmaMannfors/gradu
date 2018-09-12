################################################################################
#		Distance to field vs. type of YSO
################################################################################
#				Aug 1, 2018
################################################################################

import numpy as np
import os
import matplotlib.pyplot as plt

#################################
# Paths
#################################
distFile = '/home/emma/gradu/data/All_data/SCUBA_850/distances_to_sources/distances_to_SCUBA.txt'

dir_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/no_duplicates/'
#################################
# Code
#################################
distList = [] #x-axis
derrList = [] # xerr

FISList = []	
PSCList = []
YList = []
YcandList = []
typI_II_List = []
typIIIList = []

with open(distFile) as f: 
	for i in range(3):
		next(f)
	for line in f: 
		p = line.split()
		field = p[1]
		fld = field.split('_')
		FIELD = fld[0]+'.'+fld[1]+'.'+fld[2]
		dist = float(p[2])

		if p[3] != 'pm':
			pl = float(p[4])
			m = float(p[6])
			derr = (m+pl)/2.0

		else:
			pl = float(p[4])
			m = float(p[4])
			derr = pl
		
		#print pl,m
		fName = FIELD+'_no_duplicates.txt'
		#print fName
		
		distList.append(dist)
		derrList.append(derr)

		#if fName != 'G010.21+02.40_no_duplicates.txt':
		#	continue

		d = open(dir_in+fName,'r')

		for line in d: 
			if line[0] == '#' and line[1] == '#':
				continue

			e = line.split()

			if line[0] == '#' and e[1] == 'Name':
				continue

			if e[0] == 'match':
				continue

			if line[0] == '#':
				clump = e[4]
			else: 
				clump = e[3]

			if clump == 'nan': # Taking only clump-associated YSOs
				continue

			#print clump

		
			if line[0] == '#':
				typ = e[5]
			else: 
				typ = e[4] 

			#print typ
			if typ == 'FIS':
				FISList.append(dist)
			elif typ == 'PSC':
				PSCList.append(dist)
			elif typ == 'sim_Y*?':
				YcandList.append(dist)
			elif typ == 'sim_Y*O':
				YList.append(dist)
			elif typ == 'viz_III':
				typIIIList.append(dist)
			#elif typ == 'viz_':
			#	YcandList.append(dist)
			else: 
				print typ

#for item in FISList: 
#	print item

FISavg = np.mean(FISList)
PSCavg = np.mean(PSCList)
Ycandavg = np.mean(YcandList)
Yavg = np.mean(YList)
typIIIavg = np.mean(typIIIList)


FISstd = np.std(FISList)
PSCstd = np.std(PSCList)
Ycandstd = np.std(YcandList)
Ystd = np.std(YList)
typIIIstd = np.std(typIIIList)

x = [FISList,PSCList,YcandList,YList,typIIIList]
lab = ['Akari FIS','Akari PSC','Sim YSO cand','Sim YSO','Vizier type III']
col = ['darksalmon','peru','olivedrab','royalblue','midnightblue']

nbins = 10
plt.hist(x,stacked=True,bins=nbins,label=lab,color=col)#, histtype='bar'
plt.legend(loc='upper right')
plt.xlabel('Distance [kpc]')
#plt.show()

#FISList = []	
#PSCList = []
#YList = []
#YcandList = []
#typI_II_List = []
#typIIIList = []

print 'Type of YSO & mean distance[kpc] & std[kpc]'
print 'Akari FIS & '+ str(FISavg)+' & '+str(FISstd)
print 'Akari PSC & '+ str(PSCavg)+' & '+str(PSCstd)
print 'Simbad YSO candidate & '+ str(Ycandavg)+' & '+str(Ycandstd)
print 'Simbad YSO & '+ str(Yavg)+' & '+str(Ystd)
print 'Vizier type III & '+ str(typIIIavg)+' & '+str(typIIIstd)







#ax1.hist(x, n_bins, normed=1, histtype='bar', stacked=True)
#ax1.set_title('stacked bar')













