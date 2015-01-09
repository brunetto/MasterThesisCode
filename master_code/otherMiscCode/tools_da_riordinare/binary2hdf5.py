#!/usr/bin/env python

import time
import kd3hdf5

t = time.time()

def bin2hdf5(bin_file, h5_file, tree = False):
    import tables as tb
    import modified_read_snapshots as rs

    snap = rs.read_snapshot(bin_file)
    h5f = kd3hdf5.KDTree(h5_file, 'w')
    h5f.data_store(snap['pos'])
    if tree == True:
        h5f.tree_build
    h5f.close()

def main():
    print "start"
    for i in range(0, 512):
    #for i in [0, 10, 100, 200, 511]:
        t2 = time.time()
        print "Loop ", i
        t3 = time.time()
        bin2hdf5('../binary/snap_newMillen_subidorder_067.'+str(i), '../hdf5/data_'+str(i))
        print "Loop ", i, " finished in ", time.time()-t2

    print "That's all folks, in ", time.time()-t, "!!!"

if __name__ == "__main__":
    main()
