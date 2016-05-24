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


#plt.subplot(121)
plan = plt.imread('/home/kirill/Thesis/Stage/worlds/bitmaps/cave_filled.png')
#plan = np.mean(plan, axis=2) #uncomment for circles
plan_resized = scipy.misc.imresize(plan, (gridnum,gridnum))

#plt.imshow(plan, cmap='gray')
#print np.mean(plan, axis=2)
print "Obstacle area: ", np.shape(np.where(plan_resized==0))
winners[np.where(plan_resized==0)]=-5
u,c = np.unique(winners, return_counts=True)
c[0]=0 

goals = u[np.where(c>float(gridnum)*float(gridnum)/20.0)]
print "Threshold: ", float(gridnum)*float(gridnum)/20.0
print "Goals array: ", goals
np.savetxt("goals.txt", goals, fmt='%d')

pfs_counts = c[np.where(c>float(gridnum)*float(gridnum)/20.0)]

print "PF sizes: ", pfs_counts
print "Mean value of Place Field sizes: ", np.mean(pfs_counts)
print "Variance (sqrt) of Place Field sizes: ", np.sqrt(np.var(pfs_counts))

plt.figure(1)
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
plotstring = "Mean: {0:.2f} \n Variance: {1:.2f}".format(np.mean(pfs_counts), np.sqrt(np.var(pfs_counts)))
plt.text(40, 1000, plotstring, horizontalalignment='center', verticalalignment='center', size=18)
plt.show()
