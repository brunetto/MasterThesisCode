#! /usr/bin/env python

import numpy as np
import tables as tb

slice_dim = 0.1    #Mpc/h
all_gif = tb.openFile('../gif2_sorted.h5', 'r')
all_data = all_gif.root.data.read()
b_size = 110    #Mpc/h
slice_num = int(b_size / slice_dim)

def breakgif(i, slice_dim, slice_num):
    start_dim = i*slice_dim
    stop_dim = start_dim + slice_dim
    print "Loop ", i, " of ", slice_num , " between ", start_dim, " and ", stop_dim 
    pos_greater = all_data[start_dim < all_data[:,0]]
    pos = pos_greater[pos_greater[:,0] <= stop_dim]
    h5=tb.openFile('../gif2_'+str(i)+'.h5', 'w')
    h5.createArray(h5.root, 'data', pos, title='gif2_sorted_rebinned_pos')
    h5.flush()
    h5.close()



for i in xrange(slice_num):
    breakgif(i, slice_dim, slice_num)

print "Done."
