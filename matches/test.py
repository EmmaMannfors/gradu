with open('500_160_matches.txt') as f:
	next(f)
	for line in f:
		p = line.split()
		if p[0] != p[2]:
			print p[0]
		else:
			print 'i'
