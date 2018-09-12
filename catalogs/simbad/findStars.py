##########################################################################
#	Go through files created from Simbad, and find all 
#				YSOs
#		Also include YSO candidates, pre_main sequence stars
#		and candidates, and TT stars and candidates
##########################################################################
#		Also create list of all stars (except those above)
##########################################################################
#			Jul 9, 2018
##########################################################################
import os

files_in = '/home/emma/gradu/data/All_data/catalogs/simbad/'

dir_out = './'

YSO = open(dir_out+'YSO.txt','w') 
YSO.write('name			type		ra		dec\n')

stars = open(dir_out+'stars.txt','w')
stars.write('name			type		ra		dec\n')
#YSO_cand = open(dir_out+'YSO_candidates.txt','w')
#YSOList = []
#yso_cand_list = []
#pmsList = []
#TTList = []
scubaList = []
for filename in os.listdir(files_in):
	if filename.endswith('.txt'):
		scubaList.append(filename)
scuba_sort = sorted(scubaList)

for i in range(0,len(scuba_sort)):
		filename = scuba_sort[i]
		with open(files_in+filename) as f:
			# 9 lines of into
			YSO.write('\n'+5*'-'+'\n'+filename+'\n'+5*'-'+'\n')
			stars.write('\n'+5*'-'+'\n'+filename+'\n'+5*'-'+'\n')
			print 10*'*'
			for i in range(0,9):
				next(f)
			for line in f:
				if line[0] == '=':
					break
				p = line.split()
				# p[4] = type
				name = p[2]+' '+p[3]
				ra,dec = p[5]+':'+p[6]+':'+p[7],p[8]+':'+p[9]+':'+p[10]
				co = ra+'	'+dec
				if p[4] == 'Y*O':
					if p[2] == '2MASS' or p[2] == 'SSTGLMA':
						YSO.write(name+'	'+'confirmed YSO'+'	'+co+'\n')
					else:
						YSO.write(name+'		'+'confirmed YSO'+'	'+co+'\n')
				elif p[4] == 'Y*?':
					if p[2] == '2MASS' or p[2] == 'SSTGLMA':
						YSO.write(name+'	'+'candidate YSO'+'	'+co+'\n')
					else:
						YSO.write(name+'		'+'candidate YSO'+'	'+co+'\n')
				elif p[4] == 'TT*':
					if p[2] == '2MASS' or p[2] == 'SSTGLMA':
						YSO.write(name+'	'+'TT star'+'	'+co+'\n')
					else:
						YSO.write(name+'		'+'TT star'+'	'+co+'\n')
				elif p[4] == 'TT?':
					if p[2] == '2MASS' or p[2] == 'SSTGLMA':
						YSO.write(name+'	'+'candidate TT'+'	'+co+'\n')
					else:
						YSO.write(name+'		'+'candidate TT'+'	'+co+'\n')
				elif p[4] == 'pr*':
					if p[2] == '2MASS' or p[2] == 'SSTGLMA':
						YSO.write(name+'	'+'PMS star'+'	'+co+'\n')
					else:
						YSO.write(name+'		'+'PMS star'+'	'+co+'\n')
				elif p[4] == 'pr?':
					if p[2] == '2MASS' or p[2] == 'SSTGLMA':
						YSO.write(name+'	'+'candidate PMS'+'	'+co+'\n')
					else:
						YSO.write(name+'		'+'candidate PMS'+'	'+co+'\n')

				elif p[4].endswith('*'):
					if p[2] == '2MASS' or p[2] == 'SSTGLMA':
						stars.write(name+'	'+p[4]+'	'+co+'\n')
					else:
						stars.write(name+'		'+p[4]+'	'+co+'\n')

				#elif p[4] == 'X' or p[4] == 'ULX' or p[4] == 'UX?':
				#	print filename+'\n'+name+'	'+p[4]+'	'+co



#Candidate_YSO     	Y*?
# YSO     	Y*O
# pMS*     	pr*  	Pre-main sequence Star 
#TTau*     	TT*     	T Tau-type Star

# Candidate_pMS*     	pr?     	Pre-main sequence Star Candidate    
# Candidate_TTau*     	TT?     	T Tau star Candidate 



#starSymbolList = ['*','*iC','*iN','*iA','*i*','V*?','Pe*','HB*','Ae*','Em*','Be*','BS*','RG*','AB*','C*','S*','sg*','s*r','s*y','s*b','HS*','pA*','WD*']
#starSymbolList_2 = ['ZZ*','LM*','BD*','N*','OH*','CH*','WR*','PM*','HV*','V*','Ir*','Or*','RI*','Er*','Fl*','FU*','RC*','RC?','Ro*','a2*','Psr','BY*','RS*','Pu*']

#starSymbolList_3 = ['RR*','Ce*','dS*','RV*','WV*','bC*','cC*','gD*','SX*','LP*','Mi*','sr*','SN*','su*','Pl?','Pl']


 # if name.endswith('*') = star  
      
  



















