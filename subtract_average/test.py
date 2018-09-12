import numpy as np

a = np.eye(4,4)
m = np.nonzero(a>0)

print np.mean(a)

print 5*'*'

print np.mean(a[m])

catList = ['meow','purr','hiss']
print catList
catList = []
print catList
