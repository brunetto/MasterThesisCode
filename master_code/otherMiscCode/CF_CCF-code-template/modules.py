#!/usr/bin/env python

import f_dist    # to prevent segmentation fault

import sys
import numpy as np
import random as rnd
import tables as tb
import time
import matplotlib
import logging
matplotlib.use('Agg') # define display on the server
import matplotlib.pyplot as plt
from pylab import *

import kd3

# Start logger.
mlogger = logging.getLogger("Main_log.modules")

###################################################################################
#                             Config file reading                                 #
###################################################################################

def parse_config(filename):
    """Read the config file and store the variables into a dictionary.
    Thanks to http://www.decalage.info

    """
    mlogger.info("Reading config file.")

    COMMENT_CHAR = '#'
    OPTION_CHAR =  '='

    options = {}
    f = open(filename)
    for line in f:
        # First, remove comments:
        if COMMENT_CHAR in line:
            # split on comment char, keep only the part before
            line, comment = line.split(COMMENT_CHAR, 1)
        # Second, find lines with an option=value:
        if OPTION_CHAR in line:
            # split on option char:
            option, value = line.split(OPTION_CHAR, 1)
            # strip spaces:
            option = option.strip()
            value = value.strip()
            # store in dictionary:
            options[option] = value
    f.close()
    return options


###################################################################################
#                             Variables container                                 #
###################################################################################

class Bunch(object):
    """Create a dictionary object containing all the initial variables.

    """

    def __init__(self, d=None):
        if d is not None: self.__dict__.update(d)

def var(filename):
    """Read the initial variables and return a dictionary object.

    Parameters
    ==========

    filename = string, the name of the config file
    
    Returns
    =======

    Bunch(locals()) = dictionary object containing the local variables
 
    """

    # Read the config file and create a dictionary.
    options = parse_config(filename)
    
    # Common variables.
    mlogger.info("Defining common parameters.")

    file_1 = options['file_1']
    file_2 = options['file_2']
    n_sel =  int(options['n_sel'])
    path = options['path']
    r_step = int(options['r_step'])
    leafsize = int(options['leafsize'])
    strategy = options['strategy']
    r_min = float(options['r_min'])
    r_max = float(options['r_max'])
    log_file = options['log_file']
    test = None
    console = False

    return Bunch(locals())


####################################################################################
#                                      functions                                   #
####################################################################################

def hdf5_data(path, filename, n_sel=None):
    """Open the hdf5 files and return the positions array.
    
    Parameters
    ==========

    filename = string, the name of the hdf5 file containing the "data" array
    path = string, the name of the path to the file
    n_sel = the number of particles to be selected from the set (0 or None means all)

    Returns
    =======

    pos = array of floats containing the particles coordinates

    """
    t = time.time()
    if path is not None:
        filename=path+filename

    h5f = tb.openFile(filename, mode = "r")
    if n_sel == 0 or n_sel == None:
        pos = h5f.root.data.read()
    else:
        mlogger.info("Selected %s particles.", n_sel)
        size = h5f.root.data.read().shape[0]
        pos = h5f.root.data.read()[rnd.sample(xrange(size),n_sel), :]
    h5f.close()
    mlogger.info("Time to open %s is %s ", filename, time.time()-t)
    return pos

def hdf5_data_slice(path, filename, p_1_min, p_1_max, r_max, n_sel=None):
    """Open the hdf5 files and return the positions array of a
    selected slice.
    
    Parameters
    ==========

    filename = string, the name of the hdf5 file containing the "data" array
    path = string, the name of the path to the file
    p_1_min, p_1_max = minimum and maximum x  position of the other set
    r_max = maximum radius for the correlation
    n_sel = the number of particles to be selected from the set (0 or None means all)

    Returns
    =======

    pos = array of floats containing the particles coordinates

    """
    t = time.time()

    # Define the physical limits of the set.
    xstart = p_1_min - r_max - 0.001
    xstop = p_1_max + r_max + 0.001

    if path is not None:
        filename=path+filename

    h5f = tb.openFile(filename, mode = "r")
    
    # Retrieve the indexes for the slicing.
    indexes = h5f.root.indexes.read()
    start = indexes[np.maximum(np.searchsorted(indexes[:,0], xstart, side='left')-1, 0), 1]
    stop = indexes[np.minimum(np.searchsorted(indexes[:,0], xstop, side='right')+1, indexes.shape[0]-1), 2]

    mlogger.info("Slicing second set from %s to %s", xstart, xstop)

    if n_sel == 0 or n_sel == None:
        pos = h5f.root.data.read(start,stop)
    else:
        mlogger.info("Selected %s particles.", n_sel)
        size = h5f.root.data.read().shape[0]
        pos = h5f.root.data.read(start,stop)[rnd.sample(xrange(size),n_sel),:]
    h5f.close()
    mlogger.info("Time to open %s is %s ", filename, time.time()-t)
    return pos

def gif2(filename, path):
    """Open the gif2 fortran binary files and return the positions array.
    
    Parameters
    ==========

    filename = string, the name of the bin file containing the "data" array
    path = string, the name of the path to the file

    Returns
    =======

    array of floats containing the particles coordinates

    """
    import read_gif2 as rg
    return rg.read_gif2_file(path+filename)

def random_data(n_rand, dim_size, offset=[0,0,0]):
    """Create a (n,3) array with the spatial coordinates of
       random particles.
    
    Parameters
    ==========

    n_rand = integer, number of random particles 
    dim_size = dimension of the box
    offset = offset to be applied to the coordinates (in case of simulation slicing)
    
    Returns
    =======

    rand = array of floats containing the random particles coordinates

    """

    t = time.time()
    rand = np.zeros((n_rand, 3))

    # Random coordinates.
    rand[:,0] = np.random.random_sample(rand[:,0].size) * dim_size[0] + offset[0]
    rand[:,1] = np.random.random_sample(rand[:,0].size) * dim_size[1] + offset[1]
    rand[:,2] = np.random.random_sample(rand[:,0].size) * dim_size[2] + offset[2]
    mlogger.info("Time to create the random particles %s", time.time()-t)
    return rand


def binning(r_min, r_max, r_step, strategy='log_nosqrt_sort'):
    """ Create the arrays with the radii values and the shells limits.
    
     np.linspace(np.log10(r_min), np.log10(r_max), 2*r_step+1)) create
     an array with 2*r_step+1 elements, the last is 2*r_step

    Parameters
    ==========

    r_min = float, lower limit for the correlation 
    r_max = float, upper limit for the correlation 
    r_step = integer, number of bins
    strategy = strategy for the counting (affects the binning)
    
    Returns
    =======

    shell, r = shells limits and radii (corrected by sqrt) values

    """

    t = time.time()
    if strategy == 'log_nosqrt_sort':
        # Squared logaritmically spaced linear bins, exponential spaced in linear scale.
        r_arr = pow(10, 2*np.linspace(np.log10(r_min), np.log10(r_max), 2*r_step+1))
    elif strategy == 'log_nosqrt_nosort':
        r_arr = 2*np.linspace(np.log10(r_min), np.log10(r_max), 2*r_step+1)
    else:
        #mlogger.error("Error creating the radii array!!!")
        exit()
    
    # Number of radii: r_step-1.
    xrange_limit = int(np.floor(r_arr.size / 2))
        
    r = np.zeros(xrange_limit)
    shell = np.zeros(xrange_limit+1)
    
    # Separate the radii from the shells limits.
    if strategy == 'log_nosqrt_sort':
        for i in xrange(xrange_limit):
            r[i] = np.sqrt(r_arr[2*i+1])
            shell[i] = r_arr[2*i]
    elif strategy == 'log_nosqrt_nosort':
        for i in xrange(xrange_limit):
            r[i] = r_arr[2*i+1]/2
            shell[i] = r_arr[2*i]
        r = pow(10, r)
    else:
        #mlogger.error("Wrong strategy in binning!!")
        exit()

    shell[r_step] = r_arr[2*r_step]
    #mlogger.info("Time to create shell and radii %s", time.time()-t)
    return shell, r


def tree_build(data, leafsize):
    """ Create the two kdtrees.

    Parameters
    ==========

    data = float arrays containing the particles positions
    leafsize = integer, max number of elements in each leaf

    Returns
    =======

    tree = kdtree of the two sets of particles
    t2-t1 = time needed for the build
    
    """

    # Create the trees.
    t1 = time.time()
    tree  = kd3.KDTree(data, leafsize)
    t2 = time.time()
    mlogger.info("Time to build the tree %s", t2-t1)
    mlogger.info("Size of the tree %s", sys.getsizeof(tree))
    return tree, t2-t1

def count(tree_1, tree_2, shell, self_corr, strategy = 'log_nosqrt_sort'):
    """Count the couples.

    Parameters
    ==========

    tree_1, tree_2 = the two trees 

    shell = float array of the shell bins

    self_corr = bool, if to self-correlate

    strategy = string, counting algorithm

    Returns
    =======

    data_shell_counts = array of long64, the count in shell result
    trav_stats = dict, statistics of the count
    stats_array = array of statistics

    """

    # Counts the particles pairs with distance <= radii.
    t = time.time()
    data_counts, trav_stats, stats_array = tree_1.count_neighbors(tree_2, shell, strategy, self_corr)

    # Subtract subsequent spheres to obtain the counts in shells.
    data_shell_counts = np.zeros(shell.size)
    #print "spheres data counts ", data_counts 
    for i in xrange(shell.size-1):
        data_shell_counts[i] = data_counts[i+1] - data_counts[i]

    #print "shells data counts ", data_shell_counts
    mlogger.info("Time for the couples count %s", time.time()-t)  
    return data_shell_counts, trav_stats, stats_array

def correlate(particles_counts, mixed_counts_1, mixed_counts_2, random_counts):
    """Calculate the correlation with the LS estimator.

    Parameters
    ==========

    particles_counts, mixed_counts_1, mixed_counts_2, random_counts = the counts

    Returns
    =======

    TPCF, TPCF_err = correlation function with errors
    t = time needed

    """
    
    t1 = time.time()
    two_pt_corr = (particles_counts-mixed_counts_1-mixed_counts_2+random_counts) / random_counts
    mask = (two_pt_corr <= 0) | (np.invert(np.isfinite(two_pt_corr)))                           
    TPCF = np.ma.masked_where(mask, two_pt_corr)
    TPCF_err = 1/np.sqrt(TPCF)
    t = time.time()-t1
    mlogger.info("Time for the correlation %s", t)
    return TPCF, TPCF_err, t

def plot_correlation(r, TPCF, err):
    """Plot the results

    Parameters
    ==========

    r = the radii
    TPCF = the correlation function
    err = error array

    """
    
    t = time.time()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Two point correlation function')
    ax.set_xlabel('Radius [Mpc/h]')
    ax.set_ylabel('TPCF')
    ax.errorbar(r, TPCF, err)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)
    plt.savefig('serial_TPCF.png')
    mlogger.info("Time to plot the result %s", time.time()-t)

###################################################################################
#                                 Tests                                           #
###################################################################################

def dist_test(v, leafsizes=None):
    """Test for the time needed for calculate the distances in one
    leaf as a function of the leafsize.

    Parameters
    ==========

    v = class containing come parameters
    leafsizes = list of the different leafsizes to be tested

    Returns
    =======
    
    times = array of times

    """
    mlogger.info("Start distances test...")
    if leafsizes == None:
        #chunks = [300, 600, 1200, 2400, 4800, 9600, 19200]
        chunks = [400, 600, 1200, 2400, 4800]
    else:
        chunks = leafsizes

    pos = hdf5_data(v.file_2, v.n_sel)
    times = np.zeros((len(chunks), 2))
    for i in chunks:
        mlogger.info("Chunk dimension %s", i)
        t1 = time.time()
        result = f_dist.f_dist(i, i, pos[:i, 0], pos[:i, 1], pos[:i, 2], pos[:i, 0], pos[:i, 1], pos[:i, 2], np.empty((i*i)))
        times[chunks.index(i), 0] = time.time()-t1
        t1 = time.time()
        result.sort
        times[chunks.index(i), 1] = time.time()-t1
        if i == 600:
            mlogger.error("Time for sorting %s", times[chunks.index(i), 1])
    mlogger.info("Chunks %s", chunks)
    mlogger.info("Times %s", times)
    #print_stats_array(v.log_file+"_dist", times, 2)
    return times

def tree_test_build(data, leafsize):
    """ Create the two kdtrees for traverse testing as
    function of the leafsize.

    Parameters
    ==========

    data = float arrays containing the particles positions
    leafsizes = list of the different leafsizes to be tested

    Returns
    =======

    tree = kdtree of the two sets of particles
    t2-t1 = time needed
    
    """
    import kd3_solo_traverse as kd3
    # Create the trees.
    t1 = time.time()
    tree  = kd3.KDTree(data, leafsize)
    t2 = time.time()
    mlogger.info("Time for build the tree %s", t2-t1)
    mlogger.info("Size of the tree %s", sys.getsizeof(tree))
    return tree, t2-t1


def trav_test(v, dist=False):
    """Perform a test joining the distances calculation times with the
    traverse times as functions of the leafsizes.

    Parameters
    ==========

    v = class containing come parameters
    dist = bool, if perform the dist test

    Returns
    =======

    stats_array = array of statistics

    """
    mlogger.info("Start traverse test...")
    #tot_part = [1000000]
    tot_part = [0]
    chunks = [200, 300, 500]    
    mlogger.info("Total number of particles and chunk dimensions lists: %s, %s", tot_part, chunks)
    stats_array = np.zeros((len(tot_part)*len(chunks),10))
    if dist == True:
        mlogger.info("Starting dist test...")
        dist_array = dist_test(v, chunks)
    for n_sel in range(len(tot_part)):
        for leafsize in range(len(chunks)):
            mlogger.info("Start with particles %s and leafsize %s", tot_part[n_sel], chunks[leafsize])
            solo_trav = {}
            # Read particles.
            p_1 = hdf5_data(v.file_1, tot_part[n_sel])
            p_2 = hdf5_data(v.file_2, tot_part[n_sel])
            # Generate binning.
            shell, r = binning(v.r_min, v.r_max, v.r_step, v.strategy)
            # Build trees.
            d_tree_1, solo_trav['d_tree_1 build time'] = tree_test_build(p_1, chunks[leafsize])
            d_tree_2, solo_trav['d_tree_2 build time'] = tree_test_build(p_2, chunks[leafsize])
            del p_1, p_2
            # Counting pairs.
            mlogger.info("Start counting...")
            particles_counts, solo_trav['data'], stats_array[len(chunks)*n_sel+leafsize, :] = count(d_tree_1, d_tree_2, shell, False, v.strategy)
            mlogger.info("Total particles and chunk dimension: %s, %s", tot_part[n_sel], chunks[leafsize])
            for i in solo_trav.keys():
                mlogger.info("%s: %s", i, solo_trav[i])
    stats_array[:, 6:8] = dist_array
    # Add the expected dist time and expected total tim
    exp_dist_time = stats_array[:,4]*stats_array[:,6]
    exp_sort_time = stats_array[:,4]*stats_array[:,7]
    stats_array = np.hstack((stats_array, exp_dist_time[:, np.newaxis], exp_sort_time[:, np.newaxis]))
    exp_tot_time = exp_dist_time+exp_sort_time+stats_array[:,8]
    stats_array = np.hstack((stats_array, exp_tot_time[:, np.newaxis]))
    print_stats_array(v.log_file+"_trav", stats_array)
    return stats_array

def print_stats_array(logfile, stats_array, header=1):   
    """Print the statistics in a file.

    Parameters
    ==========

    logfile = the name of the file where to save
    stats_array = the array containing the statistics
    header = type of header
    
    """
    file_dat = open(logfile+".dat", 'w')
    if header == 1:
        header = """#TESTS ON PICO, rev25
#TRAVERSE
#tot_part, leafsize, trav_done, trav_expected, opened_leaves, expected_leaves, dist_time, sort_time, time[s], speed[node/second], expected_dist_time, expected_sort_time, expected_tot_time
"""
    else:
        header = """#TESTS ON PICO, rev25
#DIST+SORT
#dist_time, sort_time
"""        
    file_dat.write(header)
    np.savetxt(file_dat, stats_array)
    file_dat.flush()
    file_dat.close()

def plot_trav_dat(trav_file, 
             machine_name='Pico',
             ):
    """Plot the kd3 traverse data.

    Parameters
    ==========

    trav_file = the name of the file containing the statistics
    machine_name = name of the machine where the test was done

    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as font_manager

    trav_dat = np.genfromtxt(trav_file, dtype='f8', comments='#', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12))
    tot_part = trav_dat[:,0]
    leafsize = trav_dat[:,1]
    trav_done = trav_dat[:,2]
    trav_expected = trav_dat[:,3]
    opened_leaves = trav_dat[:,4]
    expected_leaves = trav_dat[:,5]
    dist_time = trav_dat[:,6]
    sort_time = trav_dat[:,7]
    trav_time = trav_dat[:,8]
    trav_speed = trav_dat[:,9]
    expected_dist_time = trav_dat[:,10]
    expected_sort_time = trav_dat[:,11]
    expected_tot_time = trav_dat[:,12]

    leaves_ratio = opened_leaves*1./expected_leaves

    mlogger.info("open figure")
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_title('Traverse stats @'+machine_name)
    ax1.set_xlabel('Leafsize')
    ax1.set_ylabel('Time [s]')
    #ax1.set_xticks(range(leafsize.size))
    #ax1.set_xticklabels(leafsize)
    ax1.set_yscale('linear')
    ax1.set_xscale('linear')
    ax1.grid(True)

    ax1.plot(leafsize, trav_time, color = 'blue', linestyle = '-.', marker = 'x', label = 'traverse time '+str(tot_part[0]))
    ax1.plot(leafsize, expected_sort_time, color = 'green', linestyle = ':', marker = 'x', label = 'expected sort time '+str(tot_part[0]))
    ax1.plot(leafsize, expected_dist_time, color = 'red', linestyle = '--', marker = 'x', label = 'expected dist time '+str(tot_part[0]))
    ax1.plot(leafsize, expected_tot_time, color = 'black', linestyle = '-', marker = 'x', label = 'expected tot time '+str(tot_part[0]))

    ax1.legend(loc='best', prop=font_manager.FontProperties(size=7))
    image_name = open(trav_file+"trav_stats_lin.png", 'w') 
    plt.savefig(image_name)

    ax1.set_yscale('log')
    ax1.set_xscale('log')

    image_name = open(trav_file+"trav_stats_log.png", 'w')
    plt.savefig(image_name)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_title('Traverse stats @'+machine_name)
    ax2.set_xlabel('Leafsize')
    ax2.set_ylabel('Opened leaves/Total leaves')
    #ax2.set_xticks(range(leafsize.size))
    #ax2.set_xticklabels(leafsize)
    ax2.set_yscale('linear')
    ax2.set_xscale('linear')
    ax2.grid(True)

    ax2.plot(leafsize, leaves_ratio, color = 'blue', linestyle = '-', marker = 'x', label = 'opened leaves ratio'+str(tot_part[0]))

    ax2.legend(loc='best', prop=font_manager.FontProperties(size=7))
    image_name = open(trav_file+"opened_leaves_ratio_lin.png",'w')
    plt.savefig(image_name)

    ax2.set_yscale('log')
    ax2.set_xscale('log')
    image_name = open(trav_file+"opened_leaves_ratio_log.png",'w')
    plt.savefig(image_name)

