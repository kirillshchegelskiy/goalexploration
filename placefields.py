import numpy as np
import matplotlib.pyplot as plt
import scipy.misc
from PIL import Image

gridnum = 127
placecells = 50

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

#plt.subplot(121)
plan = plt.imread('/home/kirill/Thesis/Stage/worlds/bitmaps/circles_filled.png')
plan1bit = np.mean(plan, axis=2)
plan_resized = scipy.misc.imresize(plan1bit, (gridnum,gridnum))

#plt.imshow(plan, cmap='gray')
#print np.mean(plan, axis=2)
winners[np.where(plan_resized==0)]=-10
u,c = np.unique(winners, return_counts=True)
c[0]=0 

img = plt.imshow(winners, interpolation='none', cmap='hot')
plt.colorbar(img)

#plt.subplot(122)
plt.figure(2)
plt.plot(u,c,'ro')
plt.grid(True)
plt.xticks(np.arange(min(u)-1, max(u), 10.0))
plt.xlim([0,placecells])
plt.ylabel('activations on grid')
plt.xlabel('place cell number')

plt.show()
