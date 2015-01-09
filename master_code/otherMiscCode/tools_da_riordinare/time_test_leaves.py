#!/usr/bin/env python

import numpy as np
import kd3
import time

def timeami(filename):
	snap = read_snapshots.read_snapshot(filename)
	pos = snap['pos']
	del snap
	size = pos.shape[0]
	for i in [10, 100, 200, 300, 400, 500, 1000, 2000, 5000, 8000, 10000, 30000, 50000, 60000, 70000, 100000, 200000]:
		t1 = time.time()
		print "Chunk dimension ", i
		ds = kd3.minkowski_distance_p(pos[0:i,:][:,np.newaxis,:], pos[0:i, :][np.newaxis,:,:]).ravel()
		print "Time to calculate the distance ", time.time()-t1, " tot ", (time.time()-t1)*size
		t2 = time.time()
		ds.sort()
		print "Time to sort ", time.time()-t2, " tot ", (time.time()-t2)*size
		print "Time for all ", time.time()-t1, " tot ", (time.time()-t1)*size
