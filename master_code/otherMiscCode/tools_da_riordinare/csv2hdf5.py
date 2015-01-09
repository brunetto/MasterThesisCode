#!/usr/bin/env python

import numpy as np
import tables as tb
import time

t = time.time()

fofhdf5 = tb.openFile('mill2_fof_snap67.h5', 'w')

fof_data = np.genfromtxt('fof0.csv', dtype=([('fofId', 'i8'), ('np', 'i4'), ('mass', 'f4'), ('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('ix', 'i4'), ('iy', 'i4'), ('iz', 'i4'), ('m_crit_200', 'f4'), ('r_crit_200', 'f4'), ('m_mean_200', 'f4'), ('r_meam_200', 'f4'), ('m_tophat', 'f4'), ('r_tophat', 'f4'), ('numSubs', 'i4')]), comments='#', delimiter=',', skiprows=26)


table = fofhdf5.createTable(fofhdf5.root, description=fof_data, name='fof_data_snap67', title="fof_data_snap67", expectedrows=11697806)

for i in range(1, 20):
    fof_data = np.genfromtxt('fof'+str(i)+'.csv', dtype=([('fofId', 'i8'), ('np', 'i4'), ('mass', 'f4'), ('x', 'f4'), ('y', 'f4'), ('z', 'f4'), ('ix', 'i4'), ('iy', 'i4'), ('iz', 'i4'), ('m_crit_200', 'f4'), ('r_crit_200', 'f4'), ('m_mean_200', 'f4'), ('r_meam_200', 'f4'), ('m_tophat', 'f4'), ('r_tophat', 'f4'), ('numSubs', 'i4')]), comments='#', delimiter=',', skiprows=26)
    
    table.append(fof_data)
    table.flush()
    print "Loop ", i, " done."


fofhdf5.close()
print "Done in ", time.time()-t, "secondi !!!"
