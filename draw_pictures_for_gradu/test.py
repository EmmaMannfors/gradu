import os

f = open('Herschel_centers.txt','w')
hList = []
for filename in os.listdir('/home/emma/gradu/data/All_data/SPIRE_PSW_250_Herschel/cropped_250'):
	if filename.endswith('crop.fits'):
		hList.append(filename)
	elif filename.endswith('crop_2.fits'):
		hList.append(filename)

	elif filename.endswith('crop_3.fits'):
		hList.append(filename)
		
	else: 
		continue

hList_s = sorted(hList)
for line in hList_s: 
	f.write(line+'\n')
