###############################################################################
#		 Calculate number of stars associated with clumps
#			vs. overall stars in field
###############################################################################
#				jul 27, 2018
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
all_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/no_duplicates/'
clumps_in = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/no_duplicates/no_nan_or_duplicates/'

dir_out = '/home/emma/gradu/data/All_data/YSO/clump_and_YSO_list/clumps_vs_all/'
######################################
# Init
######################################
doc = open(dir_out+'ratio_of_stars_in_clumps.txt','w')
tex = open(dir_out+'ratios_for_latex.txt','w')
tex.write('SCUBA-2 field & Total YSOs & YSOs associated with clumps & Percent of YSOs in clumps \\\ \n\hline \n')

fieldList = []
ratioList = []
allStarList = []
clumpStarList = []
######################################
# Code
######################################


for filename in os.listdir(all_in):
	if filename.endswith('.txt'):

		##############################
		# Name of only clump-ass. file
		##############################
		clumpFile = filename.strip('.txt')+'_no_nan.txt'

		field = filename.strip('_no_duplicates.txt')
	#	print field
		fieldList.append(field)

		##############################
	#	if filename != 'G001.36+20.96_no_duplicates.txt' and filename != 'G005.93-01.00_no_duplicates.txt':
	#		continue
		##############################



		##############################
		# Write to file
		##############################
		doc.write('Field: '+field+'\n')
		
		doc.write('# All stars file: '+filename+'\n')
		doc.write('# Clump-associated stars file: '+clumpFile+'\n')
		##############################
		# Init
		##############################
		allStars = 0
		clumpStars = 0
		clDic = dict()






		################################
		# All YSOs
		################################
		with open(all_in+filename) as ALL: 
			for line in ALL:
				##########################
				# Skips irrelevant lines
				##########################
				if line[0] == '#':
					continue
				##########################
				# Another star
				##########################
				allStars += 1

				

				################################
				# List of all clumps with YSOs
				################################
				p = line.split()
				clump = p[3]
			#	print clump
				################################
				# How many stars each clump has
				# If value exists, adds 1 to it
				# if nonexistent, creates entry
				# for the clump.
				################################
				try:
					val = clDic[clump]
					val += 1
					clDic[clump] = val
				except: 
					clDic[clump] = 1


		################################
		# clump-associated YSOs
		################################
		with open(clumps_in+clumpFile) as CL:
			for line in CL:

				##########################
				# Skips irrelevant lines
				##########################
				if line[0] == '#':
					continue


				##########################
				# New star
				##########################
				clumpStars += 1



		#for clump in clDic:
		#	print clump
		#	print clDic[clump]
		#	print 5*'*'

		#print allStars
		#print clumpStars
		ratio = float(clumpStars) / float(allStars)
		ratioList.append(ratio)
		allStarList.append(allStars)
		clumpStarList.append(clumpStars)
		
		if allStars > 350: 
			print field, allStars, clumpStars, ratio
		

		spc = '		|		'
		doc.write('# Total stars in field	|	Stars in clumps		|	portion of stars in clumps\n')
		doc.write('	'+str(allStars)+spc+str(clumpStars)+spc+str(ratio)+'\n')
		doc.write('# Clump number	|	Occurrences\n')
		for clump in sorted(clDic):
			doc.write(clump+spc+str(clDic[clump])+'\n') #
			

		doc.write(60*'#'+'\n')

		################################
		# Document for tex
		################################
		per = ratio*100
		percent = round(per,3)
		tex.write(field+' & '+str(allStars)+' & '+str(clumpStars)+' & '+str(percent)+' \\\ \n')
		tex.write('\hline \n')

####################################################
# Ratios and plot ratios
####################################################
# Removes zero-values
ratio_noZero = [i for i in ratioList if i > 0]		
	
stat = open(dir_out+'statistics.txt','w')



ma = np.max(ratioList)	
mi = np.min(ratio_noZero)#[ratioList > 0])
mean = np.mean(ratioList)
mean_noZeros = np.mean(ratio_noZero)
std = np.std(ratioList)

#print ratioList
#print 10*'*'
#print ratio_noZero

print ma,mi,mean,mean_noZeros,std
# 0.75 0.116279069767 0.428785163112 0.437031031634

stat.write('	RATIOS\n')
stat.write('Max: '+str(ma)+'\nMin: '+str(mi)+'\nMean: '+str(mean)+'\nMean, ignoring zeros: '+str(mean_noZeros)+'\n')
stat.write('Standard deviation: '+str(std)+'\n')
def plot_stats():
	n_bins = 40

	# From internet
	(mu, sigma) = norm.fit(ratioList)
	n, bins, patches = plt.hist(ratioList, 40, facecolor='midnightblue', alpha=0.6) #normed=1, 
	y = mlab.normpdf( bins, mu, sigma)
	l = plt.plot(bins, y, 'k--', linewidth=2)#,label='Distribution function ($\sigma$)')# = %s)', %s(sigma))
	######

	#plt.hist(ratioList,bins=n_bins,color='midnightblue')#,orientation='horizontal')#,cumulative=True)

	plt.title('Fraction of YSOs in clumps')
	# Draw line at mean?
	# Line color must change at height of bin

	#plt.vlines(mean,0,5,color='w')
	#plt.vlines(mean,5,10)
	plt.vlines(mean,0,9)

	plt.xlabel('Fraction (clump-associated/all YSOs)')
	plt.ylabel('Occurrences')
	plt.text(mean+0.005,8,'mean ='+str(mean))#,transform=plt.transAxes)
	s =  round(sigma,3)
	m = round(mu,3)
	plt.text(0.05,8,'$\sigma$ = '+str(s)+'\n'+'$\mu$ = '+str(m))
	#sigma = std
	#y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
	#     np.exp(-0.5 * (1 / sigma * (n_bins - mean))**2))
	##ax.plot(bins, y, '--')
	#plt.plot(y,'--')

	#plt.legend(loc='upper left')

	stat.write('sigma = std: '+str(sigma)+'\n'+'mu = mean: '+str(mu)+'\n')

	plt.show()


#plot_stats()









####################################################
# MOre stats, number of stars
####################################################
stat.write('\n	ALL STARS\n')
maA = np.max(allStarList)
miA = np.min(allStarList)
meanA = np.mean(allStarList)
stdA = np.std(allStarList)
stat.write('Max: '+str(maA)+'\nMin: '+str(miA)+'\nMean: '+str(meanA)+'\n'+'Standard deviation: '+str(stdA)+'\n')








stat.write('\n	CLUMP-ASSOCIATED STARS\n')
maC = np.max(clumpStarList)
miC = np.min(clumpStarList)
meanC = np.mean(clumpStarList)
stdC = np.std(clumpStarList)
stat.write('Max: '+str(maC)+'\nMin: '+str(miC)+'\nMean: '+str(meanC)+'\n'+'Standard deviation: '+str(stdC)+'\n')


def plot_stars():
	n_bins = 40
	
	#plt.xscale('log')

	# From internet
	(mu, sigma) = norm.fit(allStarList)
	bins = n_bins
	#n, bins, patches = plt.hist(allStarList, bins, normed=1, alpha=0.5, label='Total YSOs',range=(0,350))
	#plt.hist(clumpStarList, n_bins, normed = 1,alpha=0.5, label='Clump-associated YSOs',range=(0,350))
	
	n, bins, patches = plt.hist([allStarList,clumpStarList],bins,label=['Total YSOs','Clump-associated YSOs'],range=(0,300), alpha=0.5, normed = 1)
	
	y = mlab.normpdf( bins, mu, sigma)

	(mu, sigma) = norm.fit(clumpStarList)
	f = mlab.normpdf( bins, mu, sigma)
	
	l = plt.plot(bins, y, 'b--', linewidth=2,label='best fit of all YSOs')
	plt.plot(bins, f, color='darkgreen', linewidth=2,label='best fit of clump YSOs'),#'g-', linewidth=2,label='distribution of clump stars')






	#plt.plot(bins, f, 'k--', linewidth=1)
	######
	
	#plt.hist(allStarList, n_bins, alpha=0.5, label='Total stars',range=(0,350))
	#plt.hist(clumpStarList, n_bins, alpha=0.5, label='Clump stars',range=(0,350))
	plt.legend(loc='upper right')

	plt.xlabel('Number of YSOs per field')

	# If normed = 1 is not used
	#plt.title('Total and clump-associated stars (fields with Nstar < 350)')

	# If normed = 1 is used
	plt.title('Distribution of YSOs (fields with $N_{YSO}$ < 300)')

	#plt.vlines(mean,0,15,color='blue')#,5,10,color='blue')
	# Total = blue, clump= green
	plt.show()

plot_stars()



#x = [allStarList,clumpStarList]
#fig, axes = plt.subplots(nrows=1, ncols=2)
#plt.hist(ratioList,bins=n_bins,color='midnightblue')
#ax0, ax1 = axes.flatten()
#colors = ['red', 'tan', 'lime']
#ax0.hist(allStarList, bins=n_bins,range=(0,350))#,log=True)#,orientation='horizontal')#, density=True, histtype='bar')
#ax0.legend(prop={'size': 10})
#ax0.set_title('Number of stars')

#plt.xscale('log')
#ax1.hist(clumpStarList, bins=n_bins,range=(0,350))#,log=True)#,orientation='horizontal')#, density=True, histtype='bar')
#ax1.set_title('Number of stars in clumps')
#fig.tight_layout()
#plt.show()




