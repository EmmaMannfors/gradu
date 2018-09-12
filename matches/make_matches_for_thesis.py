OUT = open('gradu_matches.txt','w')
#OUT.write('')
with open('cropMatches.txt') as crop: 
	next(crop)
	for line in crop: 
		p = line.split()
		
		# SCUBA = p[0]
		scu = p[0]
		scuba = scu.split('_')[0]
		#print scuba


		# PACS 70 = p[1]
		p70 = p[1]
		if p[1] == 'NA':
			pacs70 = p[1]
		else: 
			if not p70.endswith('_crop.fits'):
				pacs70 =  p70.split ('_')[0]+'\_ '+p70.split ('_')[5].strip('.fits')
			else: 
				pacs70 = p70.split ('_')[0]
		#print pacs70


		# PACS 100 = p[2]
		p100 = p[2]
		if p100 == 'NA':
			pacs100 = p100
		else: 
			if not p100.endswith('_crop.fits'):
				pacs100 = p100.split ('_')[0]+'\_ '+p100.split ('_')[5].strip('.fits')
			else: 
				pacs100 = p100.split ('_')[0]
		#print pacs100
		

		# PACS 160 = p[3]
		p160 = p[3]
		if p160 == 'NA':
			pacs160 = p160
		else: 
			if not p160.endswith('_crop.fits'):
				pacs160 = p160.split ('_')[0]+'\_ '+p160.split ('_')[5].strip('.fits')
			else: 
				pacs160 = p160.split ('_')[0]
		#print pacs160
		
		
		# SPIRE = p[4]
		psw = p[4]
		if psw == 'NA':
			spire = psw
		else: 
			if not psw.endswith('_crop.fits'):
				spire = psw.split ('_')[0]+'\_ '+psw.split ('_')[4].strip('.fits')
			else: 
				spire = psw.split ('_')[0]
		
		print spire



		OUT.write(scuba+' & '+spire+' & '+pacs160+' & '+pacs100+' & '+pacs70+'\\\ '+'\n')









