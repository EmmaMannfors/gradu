####################################################################
#		Plot d, M, N(H2), density, T, aspect ratio (a/b)
####################################################################
#				Dec 4, 2018
####################################################################

#import aplpy
import seaborn as sns#;sns.set
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

plt.rcParams.update({'font.size':14})


################################
# Paths
################################
d_in = '/home/emma/gradu/data/All_data/SCUBA_850/removed_edges/clumps/size_of_clumps/'	# d, aspect ratio
# G006_01+36_74_size_of_clumps.txt

M_in = '/home/emma/gradu/data/All_data/SCUBA_850/removed_edges/clumps/masses/' # M, n
# G006_01+36_74_clump_masses.txt

NH2_in = '/home/emma/gradu/data/All_data/SCUBA_850/removed_edges/clumps/clumpInfo/' # NH2, T
# G006_01+36_74_clump_info.txt



################################
# Init
################################

distList = []
aspList = []
mList = []
nList = []
nh2List = []
tList = []


#G162_46-08_69_size_of_clumps.txt : NH(2), T estimate do not exist
#G215_87-17_50_size_of_clumps.txt : NH(2), T estimate do not exist
#G107_18+05_44_size_of_clumps.txt : NH(2), T estimate do not exist
#G107_26+05_71_size_of_clumps.txt : NH(2), T estimate do not exist
#G173_15+02_40_size_of_clumps.txt : NH(2), T estimate do not exist

noSPIRE = ['G162_46-08_69','G215_87-17_50','G107_18+05_44','G107_26+05_71','G173_15+02_40']


nFields = 0

################################
# Functions
################################


def getName(feld):
	MN = field +'_clump_masses.txt'
	TN = field +'_clump_info.txt'
	return MN,TN

def arrayify(l):
	A = np.array(l).reshape(len(l),1)
	return A


################################
# Code
################################

for filename in os.listdir(d_in):
	if filename.endswith('.txt'):
	#if filename == 'G006_01+36_74_size_of_clumps.txt':
		field = filename.strip('_size_of_clumps.txt')
		MN,TN = getName(field)
	
		if field in noSPIRE: 	
			continue
		nFields += 1

		with open(d_in+filename) as dF: 
			for line in dF: 
				p = line.split()
				if p[1] == 'Distance':
					d = float(p[3])
				elif p[0] != '#':
					if filename == 'G082_40-01_84_size_of_clumps.txt' or filename == 'G082_42-01_84_size_of_clumps.txt': 
						RA,DEC = p[1],p[5]
					else: 
						RA,DEC = p[1],p[3]

					dec = float(DEC.split(',')[1])
					ra = float(RA.split('(')[1])
					if ra > dec: 
						asp = dec/ra	#ra/dec
					else: 
						asp = ra/dec	#dec/ra
					distList.append(d)
					aspList.append(asp)
				else: 
					continue
				



		with open(M_in+MN) as mF: 
			next(mF)
			for line in mF: 
				p = line.split()
				M,n = float(p[1]),float(p[2])
				mList.append(M)
				nList.append(n)





		try: 
			with open(NH2_in+TN) as tF: 
				next(tF)
				for line in tF: 
					p = line.split()
					NH2,T = float(p[1]),float(p[2])
					nh2List.append(NH2)
					tList.append(T)


		except: 
			print filename, ': NH(2), T estimate do not exist'
			continue



#G162_46-08_69_size_of_clumps.txt : NH(2), T estimate do not exist
#G215_87-17_50_size_of_clumps.txt : NH(2), T estimate do not exist
#G107_18+05_44_size_of_clumps.txt : NH(2), T estimate do not exist
#G107_26+05_71_size_of_clumps.txt : NH(2), T estimate do not exist
#G173_15+02_40_size_of_clumps.txt : NH(2), T estimate do not exist


#print np.mean(aspList)
#print np.max(aspList)
#print np.min(aspList)


print nFields

new_nh2 = []
for nh2 in nh2List: 
	new = nh2/(10e21)
	new_nh2.append(new)

new_n = []
for n in nList: 
	new = n/(10e5)
	new_n.append(new)

D = arrayify(distList)
ASP = arrayify(aspList)
MA = arrayify(mList)
T = arrayify(tList)
N = arrayify(new_nh2)
RHO = arrayify(new_n)
#ASP = np.array(aspList).reshape(len(distList),1)

ARR = np.concatenate((D,ASP,MA,T,N,RHO),axis=1)


arr = pd.DataFrame(ARR,columns = (r'$d \/ \/ \rm (kpc)$','Aspect ratio',r'$M \/ \/ \rm (M_{\odot})$','$T \\rm \/ (K) $',r'$N \rm (H_2) \/ \/ (10^{21} \,cm^{-2})$',r'$n \rm (H_2) \/ \/ (cm^{-3})$')) 
#print arr

#sns.set(style="darkgrid")

#sns.scatterplot(distList,tList)
g = sns.pairplot(arr)

#rc={'axes.labelsize': 10, 'font.size': 8, 'legend.fontsize': 10.0, 'axes.titlesize': 10}
rc = {'axes.axisbelow': True, 'axes.grid': True, 'axes.spines.right': True, 'axes.spines.top': True}
sns.set(rc=rc)
#sns.set_style(rc) 

print sns.axes_style()
for ax in g.axes.flat: 	
    plt.setp(ax.get_xticklabels(), rotation=65)		# rotates x-axis labels by 45 degrees
    #plt.setp(ax.set_xlim(-0.1))
    #plt.setp(ax.set_ylim(-0.1))

g.map_upper(sns.scatterplot) # sns.regplot, sns.residplot, sns.kdeplot = density plot
g.map_lower(sns.scatterplot) 
g.map_diag(plt.hist) 

#plt.savefig('/home/emma/gradu/Thesis_file/picture_folder/multiPlot.pdf')
plt.savefig('/home/emma/gradu/A_A_paper/picture_folder/multiPlot.pdf')
plt.show()




#{'axes.axisbelow': False, 'font.sans-serif': [u'Bitstream Vera Sans', u'DejaVu Sans', u'Lucida Grande', u'Verdana', u'Geneva', u'Lucid', u'Arial', u'Helvetica', u'Avant Garde', u'sans-serif'], 'axes.labelcolor': u'k', 'font.family': [u'sans-serif'], 'axes.grid': False, 'axes.edgecolor': u'k', 'grid.color': u'k', 'patch.edgecolor': u'k', 'ytick.color': u'k', 'figure.facecolor': u'0.75', 'xtick.color': u'k', 'xtick.direction': u'in', 'lines.solid_capstyle': u'projecting', 'grid.linestyle': u':', 'image.cmap': u'jet', 'axes.facecolor': u'w', 'text.color': u'k', 'ytick.direction': u'in'}












