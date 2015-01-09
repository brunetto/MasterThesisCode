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

for i in range(0, 1000, 100):
    print "Start loop ", i
    filename = '../hdf5_sorted/mill2sort-'+str(i)
    h5 = tb.openFile(filename, 'r')
    data = h5.root.data.read()
    h5.close()
    data = data[rnd.sample(xrange(data.shape[0]),100), :]
    ax.scatter(data[:,0], data[:,1], data[:,2], color=colors[i%6])    

plt.savefig('mill2_100_x_blocco_10_blocchi')
print "Done in ", time.time()-t
