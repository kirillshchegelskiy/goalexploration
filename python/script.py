import matplotlib.pyplot as plt
import numpy as np

arr = imread('cave.png')
for i in range(0,499,20): 
    arr[i:i+2,:] = 0
    arr[:,i:i+2] = 0
 
plt.imshow(arr, cmap = 'gray')
plt.show()


