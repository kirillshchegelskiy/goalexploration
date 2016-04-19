import numpy as np
import matplotlib.pyplot as plt

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

print winners

plt.figure(1)

plt.subplot(121)
plan = plt.imread('/home/kirill/Thesis/Stage/worlds/bitmaps/cave.png')
plt.imshow(plan, cmap='gray')

plt.subplot(122)
img = plt.imshow(winners, interpolation='none')
plt.colorbar(img)

plt.show()
