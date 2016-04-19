import numpy as np
import matplotlib.pyplot as plt
import scipy.misc

gridnum = 127
placecells = 100

f = open('winners.txt','r')
winners = np.zeros((gridnum,gridnum))
i=1
for l in f:
	#print i, int(l)
	
	if(i%gridnum!=0): 
		winners[gridnum - i%gridnum, (i - i%gridnum)/gridnum] = int(l)	
		#print i, gridnum - i%gridnum, (i - i%gridnum)/gridnum
	elif (i%gridnum==0): 
		winners[0, i/gridnum-1] = int(l)
		#print i, 0, i/gridnum-1
	i+=1
f.close()

#print winners


plt.figure(1)

plt.subplot(121)
plan = plt.imread('/home/kirill/Thesis/Stage/worlds/bitmaps/cave_filled.png')
plan127 = scipy.misc.imresize(plan, (127,127))

plt.imshow(plan, cmap='gray')

winners[np.where(plan127==0)]=-10

plt.subplot(122)
img = plt.imshow(winners, interpolation='none', cmap='hot')
plt.colorbar(img)

plt.show()
