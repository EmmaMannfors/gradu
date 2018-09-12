################################################################################
#		Plot distance to field vs. real size of clumps
#			use average clump size 
################################################################################
#				Aug 1, 2018
################################################################################

import os
import matplotlib.pyplot as plt

#################################
# Paths
#################################
dir_in = '/home/emma/gradu/data/All_data/SCUBA_850/distances_to_sources/size_of_clumps/'

GCC_in = '/home/emma/gradu/data/All_data/SCUBA_850/distances_to_sources/size_of_clumps/GCC_cat/'
#################################
# Code
#################################
distList = []
sizeList = []
yerrList = []	# field size error
xerrList = []	# distance error

GCCdistList = []
GCCsizeList = []
GCCyerrList = []	# field size error
GCCxerrList = []

for filename in os.listdir(dir_in):
	if filename.endswith('.txt'):
	#if filename.endswith('G006_01+36_74_size_of_clumps.txt'):
		f = open(dir_in+filename,'r')

		#print filename

		next(f)
		for line in f: 
			q = line.split()
			if q[1] == 'Distance' :
				dist = float(q[3])
				xerr = float(q[5])
			if q[0] == '#':
				continue

			# ERROR BARS
			
			#print dist
			distList.append(dist)
			xerrList.append(xerr)
			

			if q[2] != 'pm':
				x = float(q[1].split('(')[1])
				y = float(q[5].split(',')[1])
				raerr  = (float(q[3])+float(q[5].split(',')[0]))/2
				decerr = (float(q[7])+float(q[9].strip(')')))/2
				yerr = (raerr+decerr)/2
				#print xerr


#[0]		[1] [2] q[3] - [5] + q[7] - 0.608)

			else: 
				x = float(q[1].split('(')[1])
				#print x
				
				y = float(q[3].split(',')[1])
				#print y

				yerr = (float(q[3].split(',')[0])+float(q[5].strip(')')))/2
				#if filename == 'G006_01+36_74_size_of_clumps.txt':
				#	print float(q[3].split(',')[0])

			# ERROR BARS

			

			avgSize = (x+y)/2

			#print avgSize
			

			sizeList.append(avgSize)
			yerrList.append(xerr)
			


for filename in os.listdir(GCC_in):
	if filename.endswith('.txt'):
		f = open(GCC_in+filename,'r')
		#print filename

		next(f)
		for line in f: 
			q = line.split()
			if q[1] == 'Distance' :
				dist = float(q[3])
			if q[0] == '#':
				continue

			
			#print dist
			GCCdistList.append(dist)
			GCCxerrList.append(0.0)
			
			xy = q[1]
			x = float(xy.split(',')[0].split('(')[1])
			y = float(xy.split(',')[1].strip(')'))
			#y = float(q[1].split(',')[1])

			#print x,y

			avg = (x+y)/2.0
			GCCsizeList.append(avg)
			GCCyerrList.append(0.0)

			print dist
			


#for item in yerrList: 
#	print item
#	print type(item)
#print 5*'*'
#for item in distList: #
#	print item
#	print type(item)

#for i in range(len(xerrList)):
#	print sizeList[i]
#	print xerrList[i]
#	print 5*'*'
#
#print len(xerrList)
#print len(sizeList)

#plt.plot(distList,sizeList,'ro')#,range=(0,2.7))
plt.errorbar(distList,sizeList,yerr=yerrList,fmt='.',color='r',label='GCCIV distances')#,'ro')#,range=#(0,2.7))#,xerr=xerrList
plt.errorbar(GCCdistList,GCCsizeList,yerr=GCCyerrList,fmt='.',color='b',label='GCC catalog distances')
plt.xlabel('Distance to field [kpc]')
plt.legend(loc='upper left')
plt.ylabel('Size of clump [pc]')
#plt.title('Clump sizes ')
plt.ylim(0,25)
plt.xlim(0,4.7)
plt.show()
# clumpSize_vs_distance.png


#if q[1] == 'Distance':
#					dist = q[3]
#				else: 	
#					continue




