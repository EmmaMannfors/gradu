import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm







def PLOT1(): 
	# Histogram of temperature and spectral index beta
	# FIGURE 3 in my paper

	tempS,bS = [],[]			# Lists of temperatures and beta for my dataset 1
	tempP,bP = [],[]			# Lists of temperatures and beta for my dataset 2
	# Creating random data, original lists were from my actual FITS images
	tempS,bS = np.random.randint(5,35,size=1000), np.random.uniform(1.0,3.0,size=1000)
	tempP,bP = np.random.randint(15,45,size=1000), np.random.uniform(1.0,3.8,size=1000)


	labs = ['Temperature (K)','$\\beta$)']	# Labels
	cols = ['violet','blue']		# Set colors of histograms


	fig, axs = plt.subplots(1,2,sharey=True)# sharey = both plots share the same y-axis
	b = np.arange(5,45,3)			# Set bin range 
	axs[0].hist(tempP,bins=b,color=cols[0],histtype=u'stepfilled',label='160-500 $\mu$m')
	axs[0].hist(tempS,bins=b,color=cols[1],histtype='step', hatch='/',linewidth=1, edgecolor='k',fill=False,label='250-500 $\mu$m')


	b = np.arange(0.5,4,0.25)		# Set bin range again
	axs[1].hist(bP,bins=b,color=cols[0],histtype=u'stepfilled',label='160-500 $\mu$m')
	axs[1].hist(bS,bins=b,color=cols[1],histtype=u'step',hatch='/',linewidth=1, edgecolor='k',fill=False,label='250-500 $\mu$m')

	
	axs[0].set_ylabel('Pixels')
	axs[1].set_xlabel('$\\beta$')
	axs[0].set_xlabel('T (K)')


	axs[0].text(0.8, 0.5, 'a', horizontalalignment='center', verticalalignment='center', transform=axs[0].transAxes,size=16)
	axs[1].text(0.8, 0.5, 'b', horizontalalignment='center', verticalalignment='center', transform=axs[1].transAxes,size=16)


	axs[0].legend(loc='upper left')
	#plt.savefig('sdfghgj',bbox_inches='tight')


	plt.show()



def beta30(T):
	A,alpha = 1.71, 1.33
	return A*(T/20.0)**(-1.0*alpha)

def betaarch(T): 
	return 11.5*(T**-0.66)

def betapro(T): 
	return 1.0/(0.4+(0.008*T))

def result(T,A,alpha):
	return A*(T**-alpha)



def PLOT2(x,y,X,Y,avgT0,avgB0,avgT1,avgB1):

	xxx,yyy = X+x,Y+y

	s = 18
	lab = ['S','S+P','avg,S','avg,S+P']
	fig, axs = plt.subplots(2,2)#,figsize=(1,15))
	col = ['darkturquoise','#110068','yellow','deeppink']
	c = col


	things = [x,X,avgT0,avgT1]

	axs[0,0].hist(X,color='k',histtype='step',linewidth=3)
	axs[0,0].hist(X,color='#110068',label=lab[1],histtype='step',linewidth=2)
	axs[0,0].hist(x,color='k',histtype='step',linewidth=4)
	axs[0,0].hist(x,color='darkturquoise',label=lab[0],histtype='step',linewidth=2)

	axs[0,0].set_ylabel('Pixels')
	axs[0,0].set_xlabel('T (K)')
	axs[0,0].xaxis.set_label_position("top")
	axs[0,0].xaxis.tick_top()
	axs[0,0].grid(alpha=0.5)
	axs[0,0].text(0.08, 0.9,'a',horizontalalignment='center',verticalalignment='center',transform = axs[0,0].transAxes,fontsize=s)
	axs[0,0].set_yscale('log')
	axs2=axs[0,0].twinx()
	axs2.set_ylabel('Fields')
	axs2.hist([avgT1],color='k',bins=5,histtype='step',linewidth=4)
	axs2.hist([avgT0],color='k',bins=5,histtype='step',linewidth=4)
	axs2.hist([avgT1],color=['deeppink'],bins=5,label=lab[3],histtype='step',linewidth=2)
	axs2.hist([avgT0],color=['yellow'],bins=5,label=lab[2],histtype='step',linewidth=2)
	axs2.set_ylim(0,16)
	axs2.tick_params(axis='y', labelcolor='k')
	axs[0,0].legend(numpoints=1,frameon=False,prop={'size': 10},loc='upper right',bbox_to_anchor=(0.67, 1.02))
	axs2.legend(numpoints=1,frameon=False,prop={'size': 10},loc='center right',bbox_to_anchor=(1.03, 0.88))



	axs[1,1].hist(y,color='k',label=lab,histtype='step',linewidth=4)
	axs[1,1].hist(y,color='darkturquoise',label=lab,histtype='step',linewidth=2)
	axs[1,1].hist(Y,color='k',label=lab,histtype='step',linewidth=3)
	axs[1,1].hist(Y,color='#110068',label=lab,histtype='step',linewidth=2)
	axs[1,1].set_ylabel('Pixels')
	axs[1,1].set_xlabel('$\\beta$')
	axs[1,1].grid(alpha=0.5)
	axs[1,1].text(0.08, 0.9,'d',horizontalalignment='center',verticalalignment='center',transform = axs[1,1].transAxes,fontsize=s)
	axs[1,1].set_yscale('log')
	axs1=axs[1,1].twinx()
	axs1.set_ylabel('Fields')
	axs1.hist([avgB0],color=['k'],bins=3,label=lab,histtype='step',linewidth=4)
	axs1.hist([avgB1],color=['k'],bins=3,label=lab,histtype='step',linewidth=4)
	axs1.hist([avgB0],color=['yellow'],bins=3,label=lab,histtype='step',linewidth=2)
	axs1.hist([avgB1],color=['deeppink'],bins=3,label=lab,histtype='step',linewidth=2)
	axs1.set_ylim(0,20)

	Y10,X10,y10,x10 = [],[],[],[]
	for i in range(len(Y)):
		if i%10 == 0:
			Y10.append(Y[i])
			X10.append(X[i])
	for i in range(len(y)):
		if i%10 == 0:
			y10.append(y[i])
			x10.append(x[i])

	axs[0,1].plot(Y10,X10,'b.',markersize=2,color=col[1],label=lab[1])
	axs[0,1].plot(y10,x10,'r.',markersize=2,color=col[0],label=lab[0])

	axs[0,1].plot(avgB1,avgT1,'rs',markersize=6,color=col[3],label=lab[3])
	axs[0,1].plot(avgB0,avgT0,'rs',markersize=6,color=col[2],label=lab[2])
	
	axs[0,1].set_ylabel('T (K)')
	axs[0,1].set_xlabel('$\\beta$')
	axs[0,1].yaxis.set_label_position("right")
	axs[0,1].yaxis.tick_right()
	axs[0,1].xaxis.set_label_position("top")
	axs[0,1].xaxis.tick_top()
	axs[0,1].grid(alpha=0.5)
	axs[0,1].legend(numpoints=1,frameon=True,prop={'size': 10},ncol=2,markerscale=2)
	axs[0,1].text(0.08, 0.9,'b',horizontalalignment='center',verticalalignment='center',transform = axs[0,1].transAxes,fontsize=s,color='white')


	CC = '0.7'

	h = axs[1,0].hist2d(xxx,yyy,bins=50,norm=LogNorm(),cmap='viridis')
	plt.colorbar(h[3], ax=axs[1,0])

	t = np.arange(10,41)
	axs[1,0].plot(t,beta30(t),'k--',linewidth=2,label='(i)')
	axs[1,0].plot(t,betaarch(t),'k-.',linewidth=3,label='(ii)')
	axs[1,0].plot(t,betapro(t),'k:',linewidth=3,label='(iii)')

	AL = [7.01]
	alphaL = [0.47]
	titles=['Our fit']
	colors=['crimson']
	colors=['yellow']
	widths= [2]
	for i in range(len(AL)):
		axs[1,0].plot(t,result(t,AL[i],alphaL[i]),'y-',linewidth=widths[i]+2,color='k',linestyle = '-')
		axs[1,0].plot(t,result(t,AL[i],alphaL[i]),'y-',linewidth=widths[i],color=colors[i],linestyle = '-',label=titles[i])

	axs[1,0].set_ylabel('$\\beta$')
	axs[1,0].set_xlabel('T (K)')
	axs[1,0].grid(alpha=0.5)
	axs[1,0].set_ylim(0.5,3.5)
	axs[1,0].legend(numpoints=1,frameon=True,prop={'size': 10},ncol=2)
	axs[1,0].text(0.08, 0.9,'c',horizontalalignment='center',verticalalignment='center',transform = axs[1,0].transAxes,fontsize=s)


	"""
	####################################
	# Clump locations and YSOs
	####################################
	####################################
	# Fig 6 in paper
	####################################

	data = np.loadtxt('TB_YSO')
	arr = np.array(data)
	TLy,BLy = arr[:,0],arr[:,1]
	Ay,alphay = 13.679819676898292,0.7194112587476388
	ny = 113

	data = np.loadtxt('TB_clump')
	arr = np.array(data)
	TLc,BLc = arr[:,0],arr[:,1]
	rr = np.nonzero((TLc>8)&(BLc<3.4))
	TLc=TLc[rr]
	BLc=BLc[rr]
	Ac,alphac = 5.760967344138566,0.46250944944515876
	nc = 15920

	fff = np.arange(5,41,1)

			
	dotcol = '#ad7ec4'
	lcol = '#2b0040'
	axs[2,0].plot(TLy,BLy,'g.',color='k',markersize=12)
	axs[2,0].plot(TLy,BLy,'g.',color=dotcol,markersize=8,label='$n$ = 113')
	axs[2,0].plot(fff,func(fff,Ay,alphay),'k-',linewidth=2)
	axs[2,0].set_xlabel('Temperature at YSO (K)') 
	axs[2,0].set_ylabel('$\\beta$')
	axs[2,0].set_ylim(0.5,3.5)
	axs[2,0].text(0.95, 0.9,'e',horizontalalignment='center',verticalalignment='center',transform = axs[2,0].transAxes,fontsize=s)

	axs[2,1].plot(TLc,BLc,'.',markersize=10,color='k',label='$n$ = 15920')
	axs[2,1].plot(TLc,BLc,'.',color=dotcol,label='$n$ = 15920')
	axs[2,1].plot(fff,func(fff,Ac,alphac),'k-',linewidth=2)
	axs[2,1].set_xlabel('Temperature at clumps (K)')
	axs[2,1].text(0.95, 0.88,'f',horizontalalignment='center',verticalalignment='center',transform = axs[2,1].transAxes,fontsize=s)




	# Legend
	custom_lines = [Line2D([0], [0], color='k', linestyle = '--',lw=4),
	                Line2D([0], [0], color='k', linestyle = '-.', lw=4),
	                Line2D([0], [0], color='k', linestyle = ':', lw=4),
	                #Line2D([0], [0], color='darkorchid', linestyle = '-.', lw=4),
	                #Line2D([0], [0], color='midnightblue', linestyle = ':', lw=4),
	                #Line2D([0], [0], color='crimson', lw=4),
	                #Line2D([0], [0], color='darkorchid', lw=4),
	                #Line2D([0], [0], color='midnightblue', lw=4),
	                Line2D([0], [0], color='yellow', lw=4),
	                Line2D([0], [0], marker = '.', color = 'w',markerfacecolor='darkturquoise',markersize = 15, lw=1),
	                Line2D([0], [0], marker = '.', color = 'w', markerfacecolor='#110068',lw=1,markersize = 15),
	                Line2D([0], [0], marker = 's', color = 'w', markerfacecolor='yellow', lw=4),
	                Line2D([0], [0], marker = 's', color = 'w', markerfacecolor='deeppink', lw=4),
	                Line2D([0], [0], marker = 's', color = 'w', markerfacecolor=dotcol, lw=4),
	                Line2D([0], [0], linestyle = '-', color = 'k', lw=2),
	                Line2D([0], [0], linestyle = '-', color = 'k', lw=2)]




	#axs[1,0].legend(custom_lines, ['1.71(T/20 K)$^{-1.33}$ (i)','11.5(T)$^{-0.66}$ (ii)', '1/(0.4+0.008T) (iii)','250-500 $\mu$m fit','160-500 $\mu$m fit', 'Both datasets', '250-500 $\mu$m','160-500 $\mu$m','average, 250-500 $\mu$m','average, 160-500 $\mu$m'],numpoints=1,ncol=2,loc='upper right',bbox_to_anchor=(2.1, -0.2))

	axs[1,0].legend(custom_lines, ['1.71(T/20 K)$^{-1.33}$ (i)','11.5(T)$^{-0.66}$ (ii)', '1/(0.4+0.008T) (iii)', 'Fit to our data', '250-500 $\mu$m','160-500 $\mu$m','average, 250-500 $\mu$m','average, 160-500 $\mu$m','160-500 $\mu$m','13.7 $T^{-0.7}$ (frame e)','5.8 $T^{-0.5}$ (frame f)'],numpoints=1,ncol=2,loc='upper right',bbox_to_anchor=(2.1, -2.))

	"""
	fig.tight_layout()

	#axs[0,0].legend(numpoints=1,ncol=2,loc='upper right',bbox_to_anchor=(1.8, -1.4))
	#axs[1,0].legend(numpoints=1,ncol=2,loc='upper right',bbox_to_anchor=(2.1, -0.2))

	#fig.set_size_inches(8, 10)

	plt.show()


def makePlot():

	####################################
	# Fig 5 in paper
	####################################
	DIR = '/home/emma/gradu/A_A_paper/codes/'

	TL0,BL0 = np.loadtxt('spire.dat',unpack=True)
	avgT0,avgB0 = np.loadtxt(DIR+'avg_spire_TB.dat',unpack=True)

	TL1,BL1  = np.loadtxt('spire_pacs.dat',unpack=True)
	avgT1,avgB1 = np.loadtxt(DIR+'avg_spire_pacs_TB.dat',unpack=True)

	PLOT2(TL0,BL0,TL1,BL1,avgT0,avgB0,avgT1,avgB1)
	





def random_histogram_function_I_found(): 
	# Legend
	custom_lines = [Line2D([0], [0], color='k', linestyle = '--',lw=4),
	                Line2D([0], [0], color='k', linestyle = '-.', lw=4),
	                Line2D([0], [0], color='k', linestyle = ':', lw=4),
	                #Line2D([0], [0], color='darkorchid', linestyle = '-.', lw=4),
	                #Line2D([0], [0], color='midnightblue', linestyle = ':', lw=4),
	                #Line2D([0], [0], color='crimson', lw=4),
	                #Line2D([0], [0], color='darkorchid', lw=4),
	                #Line2D([0], [0], color='midnightblue', lw=4),
	                Line2D([0], [0], color='yellow', lw=4),
	                Line2D([0], [0], marker = '.', color = 'w',markerfacecolor='darkturquoise',markersize = 15, lw=1),
	                Line2D([0], [0], marker = '.', color = 'w', markerfacecolor='#110068',lw=1,markersize = 15),
	                Line2D([0], [0], marker = 's', color = 'w', markerfacecolor='yellow', lw=4),
	                Line2D([0], [0], marker = 's', color = 'w', markerfacecolor='deeppink', lw=4),
	                Line2D([0], [0], marker = 's', color = 'w', markerfacecolor=dotcol, lw=4),
	                Line2D([0], [0], linestyle = '-', color = 'k', lw=2),
	                Line2D([0], [0], linestyle = '-', color = 'k', lw=2)]




	#axs[1,0].legend(custom_lines, ['1.71(T/20 K)$^{-1.33}$ (i)','11.5(T)$^{-0.66}$ (ii)', '1/(0.4+0.008T) (iii)','250-500 $\mu$m fit','160-500 $\mu$m fit', 'Both datasets', '250-500 $\mu$m','160-500 $\mu$m','average, 250-500 $\mu$m','average, 160-500 $\mu$m'],numpoints=1,ncol=2,loc='upper right',bbox_to_anchor=(2.1, -0.2))

	axs[1,0].legend(custom_lines, ['1.71(T/20 K)$^{-1.33}$ (i)','11.5(T)$^{-0.66}$ (ii)', '1/(0.4+0.008T) (iii)', 'Fit to our data', '250-500 $\mu$m','160-500 $\mu$m','average, 250-500 $\mu$m','average, 160-500 $\mu$m','160-500 $\mu$m','13.7 $T^{-0.7}$ (frame e)','5.8 $T^{-0.5}$ (frame f)'],numpoints=1,ncol=2,loc='upper right',bbox_to_anchor=(2.1, -2.))







makePlot()
#PLOT1()







