#!/usr/bin/env python

import time
import tables as tb
import numpy as np
import sys

# Global variables for patching the tree.
n_traverse = 0

####################################################################################
#                                        kd3                                       #                    
#                                                                                  #
####################################################################################

# Modified from the original by Anne
# Copyright Anne M. Archibald 2008
# Released under the scipy license

from heapq import heappush, heappop
import scipy.sparse

def minkowski_distance_p(x,y,p=2):
    """Compute the pth power of the L**p distance between x and y

    For efficiency, this function computes the L**p distance but does
    not extract the pth root. If p is 1 or infinity, this is equal to
    the actual L**p distance.  The python broadcasting rules guarantee
    that the distance is compute among alle the possible couples.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if p==np.inf:
        return np.amax(np.abs(y-x),axis=-1)
    elif p==1:
        return np.sum(np.abs(y-x),axis=-1)
    else:
        return np.sum(np.abs(y-x)**p,axis=-1)
    
def minkowski_distance(x,y,p=2):
    """Compute the L**p distance between x and y"""
    x = np.asarray(x)
    y = np.asarray(y)
    if p==np.inf or p==1:
        return minkowski_distance_p(x,y,p)
    else:
        return minkowski_distance_p(x,y,p)**(1./p)

def min_distance_rectangle(tree_1, tree_2, p=2.):
    """Compute the minimum distance between points in the two hyperrectangles."""
    tree_1_maxes = np.maximum(tree_1._v_attrs.maxes,tree_1._v_attrs.mins).astype(np.float)
    tree_1_mins = np.minimum(tree_1._v_attrs.maxes,tree_1._v_attrs.mins).astype(np.float)
    tree_2_maxes = np.maximum(tree_2._v_attrs.maxes,tree_2._v_attrs.mins).astype(np.float)
    tree_2_mins = np.minimum(tree_1._v_attrs.maxes,tree_2._v_attrs.mins).astype(np.float)
    return minkowski_distance(0, np.maximum(0,np.maximum(tree_1_mins-tree_2_maxes,tree_2_mins-tree_1_maxes)),p)

def max_distance_rectangle(tree_1, tree_2, p=2.):
    """Compute the maximum distance between points in the two hyperrectangles."""
    tree_1_maxes = np.maximum(tree_1._v_attrs.maxes,tree_1._v_attrs.mins).astype(np.float)
    tree_1_mins = np.minimum(tree_1._v_attrs.maxes,tree_1._v_attrs.mins).astype(np.float)
    tree_2_maxes = np.maximum(tree_2._v_attrs.maxes,tree_2._v_attrs.mins).astype(np.float)
    tree_2_mins = np.minimum(tree_1._v_attrs.maxes,tree_2._v_attrs.mins).astype(np.float)
    return minkowski_distance(0, np.maximum(tree_1_maxes-tree_2_mins,tree_2_maxes-tree_1_mins),p)

class KDTree(object):
    """kd-tree for quick nearest-neighbor lookup

    This class provides an index into a set of k-dimensional points
    which can be used to rapidly look up the nearest neighbors of any
    point.

    The algorithm used is described in Maneewongvatana and Mount 1999.
    The general idea is that the kd-tree is a binary trie, each of whose
    nodes represents an axis-aligned hyperrectangle. Each node specifies
    an axis and splits the set of points based on whether their coordinate
    along that axis is greater than or less than a particular value.

    During construction, the axis and splitting point are chosen by the
    "sliding midpoint" rule, which ensures that the cells do not all
    become long and thin.

    The tree can be queried for the r closest neighbors of any given point
    (optionally returning only those within some maximum distance of the
    point). It can also be queried, with a substantial gain in efficiency,
    for the r approximate closest neighbors.

    For large dimensions (20 is already large) do not expect this to run
    significantly faster than brute force. High-dimensional nearest-neighbor
    queries are a substantial open problem in computer science.

    The tree also supports all-neighbors queries, both with arrays of points
    and with other kd-trees. These do use a reasonably efficient algorithm,
    but the kd-tree is not necessarily the best data structure for this
    sort of calculation.
    """

    def __init__(self, filename, mode):
        if mode == 'read' or mode == 'r':
            self.h5file = tb.openFile(filename, mode = "r")
        elif mode == 'append' or mode == 'a':
            self.h5file = tb.openFile(filename, mode = "a")
        elif mode == 'build' or mode == 'w' or mode == 'write':
            self.h5file = tb.openFile(filename, mode = "w")
    
    def data_store(self, data):
        t = time.time()
        self.h5file.createArray(self.h5file.root, 'data', np.asarray(data), title='data')
        self.h5file.root.data._v_attrs.n_elements = data.shape[0]
        self.h5file.root.data._v_attrs.m_dimensions = data.shape[1]
        self.h5file.root.data._v_attrs.maxes = np.amax(data,axis=0)   # maxes and mins for each coord
        self.h5file.root.data._v_attrs.mins = np.amin(data,axis=0)
        print self.h5file.root.data._v_attrs.n_elements, " Stored in ", time.time()-t, "seconds."
        t = time.time()
        self.h5file.root.data.flush()
        print time.time()-t, " seconds to commit changes."

    def c_data_store(self, data):
        t = time.time()
        self.h5file.createCArray(self.h5file.root, 'data', np.asarray(data), title='data')
        self.h5file.root.data._v_attrs.n_elements = data.shape[0]
        self.h5file.root.data._v_attrs.m_dimensions = data.shape[1]
        self.h5file.root.data._v_attrs.maxes = np.amax(data,axis=0)   # maxes and mins for each coord
        self.h5file.root.data._v_attrs.mins = np.amin(data,axis=0)
        print self.h5file.root.data._v_attrs.n_elements, " Stored in ", time.time()-t, "seconds."
        t = time.time()
        self.h5file.root.data.flush()
        print time.time()-t, " seconds to commit changes."

    def tree_build(self, leafsize = 300):
        """Construct a kd-tree.
    
        Parameters:
        ===========
    
        data : array-like, shape (n,k)
            The data points to be indexed. This array is not copied, and
            so modifying this data will result in bogus results.
        leafsize : positive integer
            The number of points at which the algorithm switches over to
            brute-force.
        """
    
        print "Start building the tree..."
        t = time.time()
        #self.tree = self.h5file.createGroup(self.h5file.root, 'tree')
        self.h5file.root._v_attrs.nodetag = 0
        self.leafsize = self.h5file.root._v_attrs.leafsize = int(leafsize)     # elements in the leaves
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")
        self.__build(self.h5file.root, np.arange(self.h5file.root.data._v_attrs.n_elements), 
                     self.h5file.root.data._v_attrs.maxes, self.h5file.root.data._v_attrs.mins, 
                     'tree')   # run the build for the tree
        # Number of nodes.
        self.h5file.root.tree._v_attrs.n_nodes = self.h5file.root._v_attrs.nodetag
        self.h5file.createHardLink(self.h5file.root.tree, 'data', self.h5file.root.data)
        ### Check the definition of the most bounded particle.
        # Create a leafnode at h5file.root named h5file.root.mbp
        self.__leafnode(self.h5file.root, [0], 'mbp', self.h5file.root.tree.data.read()[0], self.h5file.root.tree.data.read()[0])

        print "Tree built in ", time.time()-t, " seconds, for ", self.h5file.root.tree._v_attrs.children, " particles."
        
        t = time.time()
        self.h5file.flush()
        print time.time()-t, " seconds to commit the changes."
        
        #return self.h5file
    
    def __leafnode(self, parent, idx, name, maxes, mins):
        self.h5file.root._v_attrs.nodetag += 1
        leaf = self.h5file.createArray(parent, name, idx, title='leafnode_'+str(self.h5file.root._v_attrs.nodetag))
        leaf._v_attrs.tag = self.h5file.root._v_attrs.nodetag    # unique (in this tree) index for the leaf        
        leaf._v_attrs.children = len(idx)   # number of elements in the leaf
        leaf._v_attrs.type = 'leaf'
        leaf._v_attrs.maxes = maxes
        leaf._v_attrs.mins = mins
        #leaf.flush()
        return leaf
    
    def __innernode(self, parent, split_dim, split, name, children, maxes, mins):
        self.h5file.root._v_attrs.nodetag +=1
        inner = self.h5file.createGroup(parent, name, title='innernode_'+str(self.h5file.root._v_attrs.nodetag))
        inner._v_attrs.split_dim = split_dim   # splitting dimension
        inner._v_attrs.split = split   # where the split happened
        inner._v_attrs.children = children    # number of elements in the node
        inner._v_attrs.tag = self.h5file.root._v_attrs.nodetag   #unique (in this tree) index of the node        
        inner._v_attrs.type = 'inner'
        inner._v_attrs.maxes = maxes
        inner._v_attrs.mins = mins
        return inner
    
    def __build(self, parent, idx, maxes, mins, name):
        if len(idx)<=self.leafsize:   
            self.__leafnode(parent, idx, name, maxes, mins)
        else:
            data = self.h5file.root.data.read()[idx]   # extract the needed from all the data
            #maxes = np.amax(data,axis=0)
            #mins = np.amin(data,axis=0)
    
            # Calculate for each axis the distance between the max
            # and the min, then argmax return the axes with the max
            # distance, so "d" is the dimension with the max range.
            d = np.argmax(maxes-mins)   
            maxval = maxes[d]   # coord of the max in the axis with max range
            minval = mins[d]   # coord of the min  in the axis with max range
            if maxval==minval:
                # all points are identical (zero dimension); warn user?
                self.__leafnode(parent, idx, name, maxes, mins)
    
            data = data[:,d] # take all the coord in the "max range" dimension
    
            # Sliding midpoint rule; see Maneewongvatana and Mount 1999
            # for arguments that this is a good idea.
            split = (maxval+minval)/2   # coord where to split (=half of the range)
    
            # Tuple containing the indexes of the non zero elements in
            # the 0-th dimension corresponding to the condition.
            less_idx = np.nonzero(data<=split)[0]   
            greater_idx = np.nonzero(data>split)[0]   
            if len(less_idx)==0:
                split = np.amin(data)   
                less_idx = np.nonzero(data<=split)[0]
                greater_idx = np.nonzero(data>split)[0]
            if len(greater_idx)==0:
                split = np.amax(data)
                less_idx = np.nonzero(data<split)[0]
                greater_idx = np.nonzero(data>=split)[0]
            if len(less_idx)==0:
                # _still_ zero? all must have the same value
                assert np.all(data==data[0]), "Troublesome data array: %s" % data
                split = data[0]
                less_idx = np.arange(len(data)-1)
                greater_idx = np.array([len(data)-1])
    
            lessmaxes = np.copy(maxes)   # lower set maxes
            lessmaxes[d] = split   # lower set ""max_range" axis" max is the split coord
            greatermins = np.copy(mins)   # higher set mins
            greatermins[d] = split   # higher set "max_range" axis min is the split coord
            inn = self.__innernode(parent, d, split, name, idx.size, maxes, mins)
            self.__build(inn._v_pathname, idx[less_idx], lessmaxes, mins, 'less') 
            self.__build(inn._v_pathname, idx[greater_idx], maxes, greatermins, 'greater')


    def __count_neighbors(self, tree_1, tree_2, r, strategy = 'log_nosqrt_sort', self_corr = False, p=2.):
        """Count how many nearby pairs can be formed.
        
        Count the number of pairs (x1,x2) can be formed, with x1 drawn
        from self and x2 drawn from other, and where distance(x1,x2,p)<=r.
        This is the "two-point correlation" described in Gray and Moore 2000,
        "N-body problems in statistical learning", and the code here is based
        on their algorithm.
        
        Parameters
        ==========
    
        other : KDTree
    
        r : float or one-dimensional array of floats The radius to
            produce a count for. Multiple radii are searched with a
            single tree traversal.  
        p : float, 1<=p<=infinity 
            Which Minkowski p-norm to use
    
        Returns
        =======
    
        result : integer or one-dimensional array of integers
            The number of pairs. Note that this is internally stored in a numpy int,
            and so may overflow if very large (two billion).
        """
        
        # Self-corr check
    
        global n_traverse
        n_traverse = 0
    
        if self_corr:
            print "Self corr ", self_corr
            tot_traverse = tree_1._v_attrs.n_nodes*(tree_1.tree._v_attrs.n_nodes-1)/2
        else:
            tot_traverse = tree_1._v_attrs.n_nodes*tree_2._v_attrs.n_nodes
    
        def traverse(node1, node2, rad_idx):
            global n_traverse
            n_traverse += 1
            if (n_traverse % 1000 == 0):
                print "traverse ", n_traverse, " di ", tot_traverse, " with ", node1._v_attrs.tag, node2._v_attrs.tag
            
            min_r = min_distance_rectangle(node1, node2, p)   # min dist between the two nodes
            max_r = max_distance_rectangle(node1, node2, p)   # min dist between the two nodes
            
            # Indexes of the radii enterely including the nodes
            # Before there was also other strategies
            if strategy == 'log_nosqrt_sort':
                 #print "calcola i raggi che includono i nodi"
                 included = r[rad_idx] > 2*np.log10(max_r)
            else:
                print "Something wrong checking the entirely including radii!!!"
                exit()
    
            # If self-corr and not yet checked nodes: sum the number of couples in the nodes.
            if (self_corr and (node1._v_attrs.tag < node2._v_attrs.tag)):    
                result[rad_idx[included]] += node1._v_attrs.children*node2._v_attrs.children   
            # If self-corr and identical nodes, add half of the couples.
            elif (self_corr and (node1._v_attrs.tag == node2._v_attrs.tag)):    
                result[rad_idx[included]] += node1._v_attrs.children*node2._v_attrs.children/2   
            # If self-corr and yet checked, drop. (now redundant)
            elif (self_corr and (node1._v_attrs.tag > node2._v_attrs.tag)):   
                pass
            # If not sef-corr, add all the couples.        
            else:
                result[rad_idx[included]] += node1._v_attrs.children*node2._v_attrs.children   
    
            if np.all(result>=0) == False:
                print "Argh!!! Negative count adding entirely included nodes!!!"
                print "result ", result
                print "nodes tag ", node1._v_attrs.tag, node2._v_attrs.tag
                print "min_r, max_r ", min_r, max_r
                print "number of particles to add ", node1._v_attrs.children*node2._v_attrs.children
                exit()
    
            # Idxs of the radii intersecting the nodes.
            if strategy == 'log_nosqrt_sort':
                if min_r == 0:
                    rad_idx = rad_idx[(0 <=r[rad_idx]) & (r[rad_idx]<=2*np.log10(max_r))]
                else:
                    rad_idx = rad_idx[(2*np.log10(min_r)<=r[rad_idx]) & (r[rad_idx]<=2*np.log10(max_r))]
            else:
                print "Something wrong checking the intersecting radii!!!"
                exit()
    
            # No radii intersecting the nodes.
            if len(rad_idx)==0:
                return
                    
            # If the first node is a leaf
            if node1._v_attrs.type == 'leaf':
                # and also the second.
                if node2._v_attrs.type == 'leaf':
    
                    ### Open leaves and count couples ###
    
                    #Before there was also other strategies
                    if strategy == 'log_nosqrt_sort':
                        # Calculate all the possible distances.                                       
                        ds = minkowski_distance_p(tree_1.data.read()[node1.read()][:,np.newaxis,:],                 
                                                tree_2.data.read()[node2.read()][np.newaxis,:,:],                
                                                p).ravel()                                            
                        ds.sort()   # sorting all the distances                                                   
                                                                                                      
                        # If is self-corr and not already checked.                                    
                        if (self_corr and (node1._v_attrs.tag <= node2._v_attrs.tag)):                             
                            # Self-corr and identical leaves: half of the result.                     
                            if node1._v_attrs.tag == node2._v_attrs.tag:                        
                                result[rad_idx] += (np.searchsorted(ds,r[rad_idx],side='right'))/2
                          
                            # Self-corr different leaves.                                             
                            else:                                                                     
                                result[rad_idx] += (np.searchsorted(ds,r[rad_idx],side='right'))
    
                        # If not self-corr.                                                           
                        elif (self_corr != True):                                               
                            result[rad_idx] += (np.searchsorted(ds,r[rad_idx],side='right'))
    
                        if np.all(result>=0) == False:
                            print "Argh!!! Negative count opening leaves"
                            print "result ", result
                            print "nodes tag ", node1._v_attrs.tag, node2._v_attrs.tag
                            print "min_r, max_r ", min_r, max_r
                            exit()
    
    ###################################
                        #print "Traverse result ", result
                            
                # First node is a leaf but second is not.
                else:
                    if (self_corr and (node1._v_attrs.tag > node2.less._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.tag, node2.less.tag
                        pass
                    else:
                        traverse(node1, node2.less, rad_idx)
    
                    if (self_corr and (node1._v_attrs.tag > node2.greater._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.tag, node2.greater.tag
                        pass
                    else:
                        traverse(node1, node2.greater, rad_idx)
            # If first node is not a leaf 
            else:
                # but second node is
                if node2._v_attrs.type == 'leaf': 
    
                    if (self_corr and (node1.less._v_attrs.tag > node2._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.less.tag, node2.tag
                        pass
                    else:
                        traverse(node1.less, node2, rad_idx)
    
                    if (self_corr and (node1.greater._v_attrs.tag > node2._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.tag
                        pass
                    else:
                        traverse(node1.greater, node2, rad_idx)
                # Second node is not a leaf
                else:
                    if (self_corr and (node1.less._v_attrs.tag > node2.less._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.less.tag, node2.less.tag
                        pass
                    else:
                        traverse(node1.less, node2.less, rad_idx)
    
                    if (self_corr and (node1.less._v_attrs.tag > node2.greater._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.tag
                        pass
                    else:
                        traverse(node1.less, node2.greater, rad_idx)
    
                    if (self_corr and (node1.greater._v_attrs.tag > node2.less._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.less.tag
                        pass
                    else:
                        traverse(node1.greater, node2.less, rad_idx)
    
                    if (self_corr and (node1.greater._v_attrs.tag > node2.greater._v_attrs.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.greater.tag
                        pass
                    else:
                        traverse(node1.greater, node2.greater, rad_idx)
            
        ############## Count_neighbours "main". ##############

        if np.shape(r) == ():
            r = np.array([r])
            result = np.zeros(1,dtype=long)
            traverse(tree_1, tree_2, np.arange(1))
            return result[0]
        elif len(np.shape(r))==1:
            r = np.asarray(r)
            n, = r.shape
            result = np.zeros(n,dtype=long)
            traverse(tree_1, tree_2, np.arange(n))
            return result
        else:
            raise ValueError("r must be either a single value or a one-dimensional array of values")
    

    def profile(self, r, strategy = 'log_nosqrt_sort', self_corr = False, p=2.):
        tree_1 = self.h5file.root.mbp
        tree_2 = other.h5file.root.tree
        t = time.time()
        counts = self.__count_neighbors(tree_1, tree_2, r, strategy = 'log_nosqrt_sort', self_corr = False, p=2.)
        print "Counting done in ", time.time()-t
        return counts
    
    def correlation(self, other, r, strategy = 'log_nosqrt_sort', self_corr = False, p=2.):
        tree_1 = self.h5file.root.tree
        tree_2 = other.h5file.root.tree
        t = time.time()
        counts = self.__count_neighbors(tree_1, tree_2, r, strategy = 'log_nosqrt_sort', self_corr = False, p=2.)
        print "Counting done in ", time.time()-t
        return counts

    def plot_node_bound(self, node):
        from enthought.mayavi import mlab
        maxes = node._v_attrs.maxes
        mins = node._v_attrs.mins

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

    def plot_node_points(self, node):
        from enthought.mayavi import mlab

        idxs = np.array([])
        for leaf in self.h5file.walkNodes(node, classname='Array'):
            idxs = np.hstack((idxs, leaf.read()))

        points = self.h5file.root.data.read()[idxs.astype(int)]
        mlab.points3d(points[:,0], points[:,1],points[:,2])

    def close(self):
        self.h5file.close()





