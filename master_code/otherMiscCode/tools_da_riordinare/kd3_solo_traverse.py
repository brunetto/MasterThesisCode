#!/usr/bin/env python

import f_dist
import numpy as np
import sys
import time
import logging

klogger = logging.getLogger("Main_log.modules.kd3")

# Global variables for patching the tree.
n_traverse = 0
n_leaves = 0
n_nodes = 0
nodetag = 0
sort_times = 0
dist_times = 0
leaves = 0

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

def minkowski_distance_p_f(x,y,p=2):
    """Fortan function to calculate the distnces.. not working for now
    """
    result = f_dist.f_dist(x.shape[0], y.shape[0], x[:,0], x[:,1], x[:,2], y[:,0], y[:,1], y[:,2], np.empty((x.shape[0]*y.shape[0], 1))) 
    return result
   
def minkowski_distance(x,y,p=2):
    """Compute the L**p distance between x and y"""
    x = np.asarray(x)
    y = np.asarray(y)
    if p==np.inf or p==1:
        return minkowski_distance_p(x,y,p)
    else:
        return minkowski_distance_p(x,y,p)**(1./p)

class Rectangle(object):
    """Hyperrectangle class.

    Represents a Cartesian product of intervals.
    """
    def __init__(self, maxes, mins):
        """Construct a hyperrectangle."""
        self.maxes = np.maximum(maxes,mins).astype(np.float)
        self.mins = np.minimum(maxes,mins).astype(np.float)

        # Number of minimi = number of dimensions of the hyper-rectangle.
        self.m, = self.maxes.shape

    def __repr__(self):
        return "<Rectangle %s>" % zip(self.mins, self.maxes)

    def volume(self):
        """Total volume."""
        return np.prod(self.maxes-self.mins)

    def split(self, d, split):
        """Produce two hyperrectangles by splitting along axis d.

        In general, if you need to compute maximum and minimum
        distances to the children, it can be done more efficiently
        by updating the maximum and minimum distances to the parent.
        """ # FIXME: do this
        mid = np.copy(self.maxes)
        mid[d] = split     #???
        less = Rectangle(self.mins, mid)
        mid = np.copy(self.mins)
        mid[d] = split    #???
        greater = Rectangle(mid, self.maxes)
        return less, greater

    def min_distance_rectangle(self, other, strategy, p=2.):
        """Compute the minimum distance between points in the two hyperrectangles."""
        return minkowski_distance_p(0, np.maximum(0,np.maximum(self.mins-other.maxes,other.mins-self.maxes)),p)

    def max_distance_rectangle(self, other, strategy, p=2.):
        """Compute the maximum distance between points in the two hyperrectangles."""
        return minkowski_distance_p(0, np.maximum(self.maxes-other.mins,other.maxes-self.mins),p)

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
    # Unique (inside the tree) index for pruning in case of
    # self-correlation.
        
    def __init__(self, data, leafsize=300):
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

        self.data = np.asarray(data)   # put the data in an array if not yet done
        self.n, self.m = np.shape(self.data)   # n = number of elements, m = dimensions
        self.leafsize = int(leafsize)     # elements in the leaves
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")
        self.maxes = np.amax(self.data,axis=0)   # maxes and mins for each coord
        self.mins = np.amin(self.data,axis=0)
        
        global nodetag
        nodetag = 0  
        global n_leaves
        n_leaves = 0
        global n_node
        n_node = 0
        
        t = time.time()
        self.tree = self.__build(np.arange(self.n), self.maxes, self.mins)   # run the build for the tree
        klogger.info("%s seconds to build a %s particles tree.", time.time()-t, self.n)

        # Fix the nested class bug.
        klogger.info("Fix the nested class bug.")
        setattr(sys.modules[__name__], 'node', KDTree.node)
        setattr(sys.modules[__name__], 'innernode', KDTree.innernode)
        setattr(sys.modules[__name__], 'leafnode', KDTree.leafnode)

    def __sizeof__(self):
        return object.__sizeof__(self) + sum(sys.getsizeof(v) for v in self.__dict__.values())

    class node(object):    # metaclass
        pass
    class leafnode(node):
        def __init__(self, idx):
            self.idx = idx   # indexes of the elements in the leaf, unique but non contiguos
            self.children = len(idx)   # number of elements in the leaf
            global nodetag
            nodetag += 1
            global n_leaves
            n_leaves += 1
            self.tag = nodetag    # unique (in this tree) index for the leaf
            #print "leafnode check assignement: ", self.tag
        def __sizeof__(self):
            return object.__sizeof__(self) + sys.getsizeof(self.idx[0])*len(self.idx)

    class innernode(node):
        def __init__(self, split_dim, split, less, greater):
            self.split_dim = split_dim   # splitting dimension
            self.split = split   # where the split happened
            self.less = less   # address of the lower half
            self.greater = greater   # address of the higher half
            self.children = less.children+greater.children    # number of elements in the node
            global nodetag
            nodetag += 1
            global n_node
            n_node += 1
            self.tag = nodetag   #unique (in this tree) index of the node
            #print "innernode check assignement: ", self.tag
        def __sizeof__(self):
            return object.__sizeof__(self) + sum(sys.getsizeof(v) for v in self.__dict__.values()) 

    def __build(self, idx, maxes, mins):
        if len(idx)<=self.leafsize:   
            return KDTree.leafnode(idx)
        else:
            data = self.data[idx]   # extract the needed from all the data
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
                return KDTree.leafnode(idx)

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
            lessmaxes[d] = split   # lower set "max_range" axis max is the split coord
            greatermins = np.copy(mins)   # higher set mins
            greatermins[d] = split   # higher set "max_range" axis min is the split coord
            return KDTree.innernode(d, split,
                                    self.__build(idx[less_idx],lessmaxes,mins),
                                    self.__build(idx[greater_idx],maxes,greatermins))

    def count_neighbors(self, other, r, strategy = 'log_nosqrt_sort', self_corr = False, p=2.):
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
        global sort_times
        sort_times = 0
        global leaves
        leaves = 0
        global dist_times
        dist_times = 0


        if self_corr:
            klogger.info("Self corr %s", self_corr)
            tot_traverse = self.tree.tag*(self.tree.tag-1)/2
        else:
            tot_traverse = self.tree.tag*other.tree.tag

        klogger.info("Start counting: total traverse %s, ")

        def traverse(node1, rect1, node2, rect2, rad_idx):
            global n_traverse
            global sort_times
            global leaves 
            global dist_times

            n_traverse += 1
            if (n_traverse % 100000 == 0):
                klogger.info("traverse %s of %s with %s, %s", n_traverse, tot_traverse, node1.tag, node2.tag)
                
            min_r = rect1.min_distance_rectangle(rect2, strategy, p)   # min dist between the two nodes
            max_r = rect1.max_distance_rectangle(rect2, strategy, p)   # min dist between the two nodes

            # Indexes of the radii enterely including the nodes
            # Before there was also other strategies
            if strategy == 'log_nosqrt_sort' or strategy == 'log_nosqrt_nosort':
#                print "calcola i raggi che includono i nodi"
                included = r[rad_idx] > max_r
            else:
                klogger.error("Something wrong checking the entirely including radii!!!")
                exit()

            # If self-corr and not yet checked nodes: sum the number of couples in the nodes.
            if (self_corr and (node1.tag < node2.tag)):    
                pass
            # If self-corr and identical nodes, add half of the couples.
            elif (self_corr and (node1.tag == node2.tag)):    
                pass
            # If self-corr and yet checked, drop. (now redundant)
            elif (self_corr and (node1.tag > node2.tag)):   
                pass
            # If not sef-corr, add all the couples.        
            else:
                pass

            if np.all(result>=0) == False:
                klogger.error("Argh!!! Negative count adding entirely included nodes!!!")
                klogger.error("result %s", result)
                klogger.error("nodes tag %s", node1.tag, node2.tag)
                klogger.error("min_r, max_r %s", min_r, max_r)
                klogger.error("number of particles to add %s", node1.children*node2.children)
                exit()
 
            # Idxs of the radii intersecting the nodes.

            if strategy == 'log_nosqrt_sort' or strategy == 'log_nosqrt_nosort':
#                if min_r == 0:
#                    rad_idx = rad_idx[(0 <=r[rad_idx]) & (r[rad_idx]<=max_r)]
#                else:
                rad_idx = rad_idx[(min_r<=r[rad_idx]) & (r[rad_idx]<=max_r)]
            else:
                klogger.error("Something wrong checking the intersecting radii!!!")
                exit()

            # No radii intersecting the nodes.
            if len(rad_idx)==0:
                return
                    
            # If the first node is a leaf
            if isinstance(node1,KDTree.leafnode):
                # and also the second.
                if isinstance(node2,KDTree.leafnode):
                    ### Open leaves and count couples ###
                    #Before there was also other strategies
                    if strategy == 'log_nosqrt_sort':
                        # Calculate all the possible distances.
                        t = time.time()
                        pass
                        sort_times += time.time()-t
                        leaves +=1 
                        # If is self-corr and not already checked.                                    
                        if (self_corr and (node1.tag <= node2.tag)):                             
                            # Self-corr and identical leaves: half of the result.                     
                            if node1.tag == node2.tag:                        
                                pass
                          
                            # Self-corr different leaves.                                             
                            else:                                                                     
                                pass
  
                        # If not self-corr.                                                           
                        elif (self_corr != True):                                               
                            pass

                        if np.all(result>=0) == False:
                            klogger.error("Argh!!! Negative count opening leaves")
                            klogger.error("result %s", result)
                            klogger.error("nodes tag %s, %s", node1.tag, node2.tag)
                            klogger.error("min_r, max_r %s, %s", min_r, max_r)
                            exit()

#######################################

                # First node is a leaf but second is not.
                else:
                    less, greater = rect2.split(node2.split_dim, node2.split)

                    if (self_corr and (node1.tag > node2.less.tag)):
                        #print "Pruning traverse on ", node1.tag, node2.less.tag
                        pass
                    else:
                        traverse(node1, rect1, node2.less, less, rad_idx)

                    if (self_corr and (node1.tag > node2.greater.tag)):
                        #print "Pruning traverse on ", node1.tag, node2.greater.tag
                        pass
                    else:
                        traverse(node1, rect1, node2.greater, greater, rad_idx)
            # If first node is not a leaf 
            else:
                # but second node is
                if isinstance(node2,KDTree.leafnode): 
                    less, greater = rect1.split(node1.split_dim, node1.split)

                    if (self_corr and (node1.less.tag > node2.tag)):
                        #print "Pruning traverse on ", node1.less.tag, node2.tag
                        pass
                    else:
                        traverse(node1.less, less, node2, rect2, rad_idx)

                    if (self_corr and (node1.greater.tag > node2.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.tag
                        pass
                    else:
                        traverse(node1.greater, greater, node2, rect2, rad_idx)
                # Second node is not a leaf
                else:
                    less1, greater1 = rect1.split(node1.split_dim, node1.split)
                    less2, greater2 = rect2.split(node2.split_dim, node2.split)

                    if (self_corr and (node1.less.tag > node2.less.tag)):
                        #print "Pruning traverse on ", node1.less.tag, node2.less.tag
                        pass
                    else:
                        traverse(node1.less,less1,node2.less,less2,rad_idx)

                    if (self_corr and (node1.less.tag > node2.greater.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.tag
                        pass
                    else:
                        traverse(node1.less,less1,node2.greater,greater2,rad_idx)

                    if (self_corr and (node1.greater.tag > node2.less.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.less.tag
                        pass
                    else:
                        traverse(node1.greater,greater1,node2.less,less2,rad_idx)

                    if (self_corr and (node1.greater.tag > node2.greater.tag)):
                        #print "Pruning traverse on ", node1.greater.tag, node2.greater.tag
                        pass
                    else:
                        traverse(node1.greater,greater1,node2.greater,greater2,rad_idx)
            
        ############## Count_neighbours "main". ##############
        R1 = Rectangle(self.maxes, self.mins)   # first tree rectangle
        R2 = Rectangle(other.maxes, other.mins)   # second tree rectangle
        
        if np.shape(r) == ():
            t1 = time.time()
            r = np.array([r])
            result = np.zeros(1,dtype=long)
            traverse(self.tree, R1, other.tree, R2, np.arange(1))
            t2 = time.time()
            trav_stats = {}
            trav_stats['Numer of particles'] = self.n
            trav_stats['Total numer of nodes'] = nodetag
            trav_stats['Number of leaves'] = n_leaves
            trav_stats['Number of innernodes'] = n_nodes
            trav_stats['Leafsize'] = self.leafsize
            trav_stats['Time for traverse'] = t2-t1
            trav_stats['sort_times'] = sort_times
            trav_stats['dist_times'] = dist_times
            trav_stats['n_opened_leaves'] = leaves
            trav_stats['Number of traverse'] = n_traverse
            trav_stats['Expected traverse'] = tot_traverse
            trav_stats['Traverse speed'] = n_traverse / (t2-t1)
            array_stats = np.zeros(10)
            array_stats[0] = self.n
            array_stats[1] = self.leafsize
            array_stats[2] = n_traverse
            array_stats[3] = tot_traverse
            array_stats[4] = leaves
            array_stats[5] = n_leaves**2
            array_stats[6] = dist_times
            array_stats[7] = sort_times
            array_stats[8] = t2-t1
            array_stats[9] = n_traverse/(t2-t1)
            klogger.info("Done with %s ans %s particles and leafsize %s: %s traverse of %s in %s seconds at %s nodes per second.", 
                         self.n, self.leafsize, n_traverse, tot_traverse, t2-t1, n_traverse/(t2-t1))
            return result[0], trav_stats, array_stats
        elif len(np.shape(r))==1:
            t1 = time.time()
            r = np.asarray(r)
            n, = r.shape
            result = np.zeros(n,dtype=long)
            traverse(self.tree, R1, other.tree, R2, np.arange(n))
            t2 = time.time()
            trav_stats = {}
            trav_stats['Numer of particles'] = self.n
            trav_stats['Total numer of nodes'] = nodetag
            trav_stats['Number of leaves'] = n_leaves
            trav_stats['Number of innernodes'] = n_nodes
            trav_stats['Leafsize'] = self.leafsize
            trav_stats['Time for traverse'] = t2-t1
            trav_stats['sort_times'] = sort_times
            trav_stats['dist_times'] = dist_times
            trav_stats['n_opened_leaves'] = leaves
            trav_stats['Number of traverse'] = n_traverse
            trav_stats['Expected traverse'] = tot_traverse
            trav_stats['Traverse speed'] = n_traverse / (t2-t1)
            array_stats = np.zeros(10)
            array_stats[0] = self.n
            array_stats[1] = self.leafsize
            array_stats[2] = n_traverse
            array_stats[3] = tot_traverse
            array_stats[4] = leaves
            array_stats[5] = n_leaves**2
            array_stats[6] = dist_times
            array_stats[7] = sort_times
            array_stats[8] = t2-t1
            array_stats[9] = n_traverse/(t2-t1)
            klogger.info("Done with %s and %s particles and leafsize %s: %s traverse of %s in %s seconds at %s nodes per second.", 
                         self.n, other.n, self.leafsize, n_traverse, tot_traverse, t2-t1, n_traverse/(t2-t1))
            return result, trav_stats, array_stats
        else:
            raise ValueError("r must be either a single value or a one-dimensional array of values")

