#! /use/bin/env python

import numpy as np
import tables as tb
from enthought.mayavi import mlab

pippo = tb.openFile('prova.h5', 'r')
maxes = pippo.root.tree._v_attrs.maxes
mins = pippo.root.tree._v_attrs.mins

x_M = maxes[0]
y_M = maxes[1]
z_M = maxes[2]
x_m = mins[0]
y_m = mins[1]
z_m = mins[2]

floor_x = np.array([x_m, x_M, x_M, x_m, x_m])
floor_y = np.array([y_m, y_m, y_M, y_M, y_m])
floor_z = np.array([z_m, z_m, z_m, z_m, z_m])



roof_x = np.array([x_m, x_M, x_M, x_m, x_m])
roof_y = np.array([y_m, y_m, y_M, y_M, y_m])
roof_z = np.array([z_M, z_M, z_M, z_M, z_M])



edge1x = np.array([x_m, x_m])
edge2x = np.array([x_M, x_M])
edge3x = np.array([x_M, x_M])
edge4x = np.array([x_m, x_m])


edge1y = np.array([y_m, y_m])
edge2y = np.array([y_m, y_m])
edge3y = np.array([y_M, y_M])
edge4y = np.array([y_M, y_M])

edge1z = np.array([z_m, z_M])
edge2z = np.array([z_m, z_M])
edge3z = np.array([z_m, z_M])
edge4z = np.array([z_m, z_M])



mlab.plot3d(floor_x, floor_y, floor_z)
mlab.plot3d(roof_x, roof_y, roof_z)

mlab.plot3d(edge1x, edge1y, edge1z)
mlab.plot3d(edge2x, edge2y, edge2z)
mlab.plot3d(edge3x, edge3y, edge3z)
mlab.plot3d(edge4x, edge4y, edge4z)

idxs = np.array([])
for leaf in self.h5file.walkNodes(node, classname='Array'):
    idxs = np.hstack((idxs, leaf.read()))

points = self.h5file.root.data.read()[idxs.astype(int)]
mlab.points3d(points[:,0], points[:,1],points[:,2])

#####################################################################
import numpy as np
import read_snapshots as rs
from enthought.mayavi import mlab
import random as rnd


#BOX

x_M = 100
y_M = 100
z_M = 100
x_m = 0
y_m = 0
z_m = 0

floor_x = np.array([x_m, x_M, x_M, x_m, x_m])
floor_y = np.array([y_m, y_m, y_M, y_M, y_m])
floor_z = np.array([z_m, z_m, z_m, z_m, z_m])



roof_x = np.array([x_m, x_M, x_M, x_m, x_m])
roof_y = np.array([y_m, y_m, y_M, y_M, y_m])
roof_z = np.array([z_M, z_M, z_M, z_M, z_M])



edge1x = np.array([x_m, x_m])
edge2x = np.array([x_M, x_M])
edge3x = np.array([x_M, x_M])
edge4x = np.array([x_m, x_m])


edge1y = np.array([y_m, y_m])
edge2y = np.array([y_m, y_m])
edge3y = np.array([y_M, y_M])
edge4y = np.array([y_M, y_M])

edge1z = np.array([z_m, z_M])
edge2z = np.array([z_m, z_M])
edge3z = np.array([z_m, z_M])
edge4z = np.array([z_m, z_M])



mlab.plot3d(floor_x, floor_y, floor_z)
mlab.plot3d(roof_x, roof_y, roof_z)

mlab.plot3d(edge1x, edge1y, edge1z)
mlab.plot3d(edge2x, edge2y, edge2z)
mlab.plot3d(edge3x, edge3y, edge3z)
mlab.plot3d(edge4x, edge4y, edge4z)

snap_1=rs.read_snapshot('snap_newMillen_subidorder_067.0')
pos_1 = snap_1['pos'][rnd.sample(xrange(snap_1['pos'].shape[0]),1000), :]
del snap_1

snap_2=rs.read_snapshot('snap_newMillen_subidorder_067.100')
pos_2 = snap_2['pos'][rnd.sample(xrange(snap_2['pos'].shape[0]),1000), :]
del snap_2

snap_3=rs.read_snapshot('snap_newMillen_subidorder_067.200')
pos_3 = snap_3['pos'][rnd.sample(xrange(snap_3['pos'].shape[0]),1000), :]
del snap_3

snap_4=rs.read_snapshot('snap_newMillen_subidorder_067.300')
pos_4 = snap_4['pos'][rnd.sample(xrange(snap_4['pos'].shape[0]),1000), :]
del snap_4

snap_5=rs.read_snapshot('snap_newMillen_subidorder_067.511')
pos_5 = snap_5['pos'][rnd.sample(xrange(snap_5['pos'].shape[0]),1000), :]
del snap_5


mlab.points3d(pos_1[:,0], pos_1[:,1],pos_1[:,2], scale_factor=2)
mlab.points3d(pos_2[:,0], pos_2[:,1],pos_2[:,2], scale_factor=2)
mlab.points3d(pos_3[:,0], pos_3[:,1],pos_3[:,2], scale_factor=2)
mlab.points3d(pos_4[:,0], pos_4[:,1],pos_4[:,2], scale_factor=2)
mlab.points3d(pos_5[:,0], pos_5[:,1],pos_5[:,2], scale_factor=2)

h5file = tb.openFile('../mill2_web/mill2_fof_snap67.h5', mode = "r")
pos=h5file.root.all_fof_pos.read()
fof = pos['pos'][rnd.sample(xrange(pos.shape[0]),200000), :]
mlab.points3d(fof[:,0], fof[:,1],fof[:,2], scale_factor=2)


********************
import matplotlib.pyplot as plt
import random as rnd
from mpl_toolkits.mplot3d import Axes3D

pos=pos[rnd.sample(xrange(pos.shape[0]),10000), :]
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(pos[:,0], pos[:,1], pos[:,2], color=colours[i%6])

