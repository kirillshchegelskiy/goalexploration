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
u,c = np.unique(winners, return_counts=True)

plt.figure(1)

#plt.subplot(121)
plan = plt.imread('/home/kirill/Thesis/Stage/worlds/bitmaps/cave_filled.png')
plan_resized = scipy.misc.imresize(plan, (gridnum,gridnum))

#plt.imshow(plan, cmap='gray')

winners[np.where(plan_resized==0)]=-10


img = plt.imshow(winners, interpolation='none', cmap='hot')
plt.colorbar(img)

#plt.subplot(122)
plt.figure(2)
plt.plot(u,c,'ro')
plt.grid(True)
plt.xticks(np.arange(min(u)-1, max(u), 10.0))
plt.ylabel('activations on grid')
plt.xlabel('place cell number')

plt.show()
