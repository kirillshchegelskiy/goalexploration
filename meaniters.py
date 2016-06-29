import numpy as np
import matplotlib.pyplot as plt

maxruns = 1000

f = open('naviters.txt','r')
filelines = f.readlines()
f.close()

iters = np.array(filelines)
iters = iters.astype(np.int)

print "Mean iterations number: ", np.mean(iters)
