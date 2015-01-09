#!/opt/epd-7.0-2-rh5-x86_64/bin/python

import numpy as np
import tables as tb
import random as rnd
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

t=time.time()
print "Start"
colors = ['black', 'blue', 'green', 'red', 'yellow', 'magenta']
fig = plt.figure()
ax = Axes3D(fig)

for i in range(0, 1000):
    filename = '../hdf5_sorted/mill2sort-'+str(i)
    h5 = tb.openFile(filename, 'r')
    data = h5.root.data.read()
    h5.close()
    print np.amin(data[:,0])
    print np.amax(data[:,0])

print "Done in ", time.time()-t
