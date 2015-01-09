#!/home/brunetto/enthougth/epd-7.0-2-rh5-x86_64/bin/python

#===============================================================================
#    Modules import.
#===============================================================================

#import f_dist    # to prevent segmentation fault

import os
import logging
import sys, argparse
import numpy as np
import time
from shutil import copy as shcopy

import modules as mod

#===============================================================================
#    MAIN                                   
#===============================================================================

def main(argv=None):

#===============================================================================
#    Parameters and arguments parsing
#===============================================================================

    tt_glob = time.time()

    # Parameters load from config file.
    v = mod.var('../config.txt')

    if argv is None:
        # Check for CLI variables.
        argv = sys.argv

        parser = argparse.ArgumentParser()
        parser.add_argument('-f1', '-file_1', '--file_1', action='store', dest='file_1', default=v.file_1, 
                    help='Filename of the first set')
        parser.add_argument('-f2', '-file_2', '--file_2', action='store', dest='file_2', default=v.file_2, 
                    help='Filename of the second set')
        parser.add_argument('-l', '--leafsize', action='store', dest='leafsize', type=int, default=v.leafsize, 
                    help='Max number of particles in a leaf')
        parser.add_argument('-n', '--num_part', action='store', dest='num_part', type=int, default=v.n_sel, 
                    help='Number of particles to be selected')
        parser.add_argument('-m', '-m_factor', '--m_factor', action='store', dest='m_factor', default=1, 
                    help='How many randoms respect to the data?')
        parser.add_argument('-s', '--strategy', action='store', dest='strategy', default=v.strategy, 
                    help='Counting strategy: log_nosqrt_sort, log_nosqrt_nosort')
        parser.add_argument('-t', '--test', action='store_true', dest='test', default=None, 
                    help='Test run.')
        parser.add_argument('-log', '-lg', '--log', action='store', dest='log_file', default=None, 
                    help='Log file basename')
        parser.add_argument('-c', '--console', action='store_true', dest='console', default=False, 
                    help='Add a console output')
        parser.add_argument('-sl', '--slicing', action='store_true', dest='slicing', default=False, 
                    help='Activate the slicing on the second set.')
        parser.add_argument('-mb', '--mass_bins', action='store', dest='mass_bin', default=None, 
                    help='Select mass bin for haloes.')
        parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

        cli = parser.parse_args()

        # Overwrite the config file parameters
        v.file_1 = cli.file_1
        v.file_2 = cli.file_2
        v.leafsize = cli.leafsize
        v.n_sel = cli.num_part
        m_factor = int(cli.m_factor)
        v.strategy = cli.strategy
        v.log_file = cli.log_file
        v.test = cli.test
        v.console = cli.console
        slicing = cli.slicing
        mass_bin = cli.mass_bin

    elif isinstance(argv, dict):
        # Reading variables passed to Main as a function (Guido docet
        # http://www.artima.com/weblogs/viewpost.jsp?thread=4829).
        for i in argv.keys():
            if i in ("-f1", "--file_1"): 
                v.file_1 = argv[i]
            elif i in ("-f2", "--file_2"): 
                v.file_2 = argv[i]
            elif i in ("-l", "--leafsize"): 
                v.leafsize = argv[i]
            elif i in ("-n", "--num_part"):
                v.n_sel = argv[i]
            elif i in ("-s", "--strategy"):
                v.strategy = argv[i]
            elif i in ("-log", "-lg", "--log"):
                v.strategy = argv[i]
            elif i in ("-t", "--test"):
                v.test = True
            elif i in ("-c", "--console"):
                v.console = argv[i]
            elif i in ("-sl", "--slicing"):
                slicing = argv[i]
            elif i in ("-mb", "--mass_bins"):
	        mass_bin = argv[i]
            else:
                print "Wrong parameter passed to main function, exit!!!"
                sys.exit(1)

#===============================================================================
#    Logger
#===============================================================================


    # Create logger.
    logger = logging.getLogger("Main_log")
    logging.captureWarnings(True)
    logger.setLevel(logging.DEBUG)
    # Create file handler which logs even debug messages.
    fh = logging.FileHandler(v.log_file+".log", 'w')
    fh.setLevel(logging.DEBUG)
    # Create formatter.
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    # Add the handlers to the logger.
    logger.addHandler(fh)
    if v.console is True:
        # Create console handler.
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    # Start logging.
    logger.info("Log started.")
    logger.info("Process PID is %s", os.getpid())
    logger.info("Copying proc status file at the beginning.")
    status_src = '/proc/'+str(os.getpid())+'/status'
    status_dest = '../logs/status'+str(os.getpid())+'@beginning'
    shcopy(status_src, status_dest)

#===============================================================================
#    Test
#===============================================================================

    # Check for test run.
    if v.test is not None:
        logger.info("Selected test %s", cli.test)
        if v.test == True:
            stats_array = mod.trav_test(v, True)
            logger.info("plotting test data...")
            mod.plot_trav_dat(v.log_file+"_trav"+".dat", v.log_file+"_dist"+".dat")
            logger.info("End test.")
            logger.info( "Wall time for all %s", time.time()-tt_glob)
            return 42
        else:
            logger.error("Something wrong in test selection.")
    else:
        logger.info("This is not a test!")

    # Dictionary for statistics.
    stats = {}

#===============================================================================
#    Data reading or creation
#===============================================================================

    # Retrieve the positions.
    p_1 = mod.hdf5_data(v.path, v.file_1, v.n_sel)

    # Binary gif files.
    #p_1 = mod.gif2(v.file_1)
    #p_2 = mod.gif2(v.file_2)

    # Random data for test.
    #p_1 = mod.random_data(10000, b_size_1, 0) 
    #p_2 = mod.random_data(10000, b_size_1, 0) 

    # Usefull spatial info about the first set.
    min_1 = np.amin(p_1, 0)
    max_1 = np.amax(p_1, 0)
    b_size_1 = max_1 - min_1
    offset_1 = np.array([min_1[0],0,0])
 
    if slicing == False:
        logger.info("slicing is off.")
        p_2 = mod.hdf5_data(v.path, v.file_2, v.n_sel)
        if mass_bin != None:
	  logger.info("but selecting mass bin %s", mass_bin)
          p_2 = mod.hdf5_data(v.path, v.file_2, None, mass_bin)
    elif slicing == True:
        logger.info("slicing is on.")
        p_2 = mod.hdf5_data_slice(v.path, v.file_2, min_1[0], 
                                  max_1[0], v.r_max, v.n_sel)
    else: 
        print "problem with slicing choice"

    # Usefull spatial info about the second set.    
    min_2 = np.amin(p_2, 0)
    max_2 = np.amax(p_2, 0)
    b_size_2 = max_2 - min_2
    offset_2 = np.array([min_2[0],0,0])

    logger.info("First set limits %s %s", min_1, max_1)
    logger.info("Second set limits %s %s", min_2, max_2)

    stats['Set 1'] = p_1.shape[0]
    stats['Set 2'] = p_2.shape[0]

    logger.info("path is %s", v.path)
    logger.info("filenames %s, %s", v.file_1, v.file_2)

    # Check for data self correlation.
    if p_1.size == p_2.size:
        self_corr = (p_1==p_2).all()
    else:
        self_corr = False
    
    logger.info("Data self correlation is %s", self_corr)
    
    logger.info("We are going to create random %s * dim_data.", m_factor)

    logger.info("Correlating on n_bin = %s", v.r_step)
    logger.info("From [Mpc] %s", v.r_min)
    logger.info("To [Mpc] %s", v.r_max)

    # Generate random.
    random_1 = mod.random_data(p_1.shape[0]*m_factor, b_size_1, offset_1) #fof
    if self_corr == True:
        random_2 = random_1
    else:
        random_2 = mod.random_data(p_2.shape[0]*m_factor, b_size_2, offset_2) #particles

    # Create the result files. 
    result_file = open(v.log_file+'-result.dat', 'a')

#===============================================================================
#    Binning
#===============================================================================


    # Generate binning.
    logger.info("Binning...")
    shell, r = mod.binning(v.r_min, v.r_max, v.r_step, v.strategy)

    # Save binning.
    result_file.write("Shells ")
    np.savetxt(result_file, shell[np.newaxis,:])
    result_file.write("Radii ")
    np.savetxt(result_file, r[np.newaxis,:])
    result_file.flush()

#===============================================================================
#    Trees build
#===============================================================================

    # Build trees.
    logger.info("Start building the trees...")
    logger.info("Leafsize %s", v.leafsize)
    d_tree_1, stats['d_tree_1 build time'] = mod.tree_build(p_1, v.leafsize)

    # Data trees.
    if self_corr == True:
        logger.info("Copying the first data tree into the second.")
        d_tree_2 = d_tree_1
        stats['d_tree_2 build time'] = stats['d_tree_1 build time']
    else:
        d_tree_2, stats['d_tree_2 build time'] = mod.tree_build(p_2, v.leafsize)

    # Random trees.
    r_tree_1, stats['r_tree_1 build time'] = mod.tree_build(random_1, v.leafsize)
    if self_corr == True:
        logger.info("Copying the first random tree into the second.")
        r_tree_2 = r_tree_1
        stats['r_tree_2 build time'] = stats['r_tree_1 build time']
    else:
        r_tree_2, stats['r_tree_2 build time'] = mod.tree_build(random_2, v.leafsize)

    del p_1, p_2, random_1, random_2 # da controllare, il random non andra` cancellato poi

    logger.info("Copying proc status file at mid.")
    status_src = '/proc/'+str(os.getpid())+'/status'
    status_dest = '../logs/status'+str(os.getpid())+'@mid'
    shcopy(status_src, status_dest)


    stats['Time elapsed before traversing'] = time.time()-tt_glob
    logger.info("Time elapsed before traversing %s", stats['Time elapsed before traversing'])

#===============================================================================
#    Counting
#===============================================================================

    # Data counting.
    trav_stats = {}
    stats_array = np.zeros((4, 10))
    logger.info("Starting data-data counts...")

    particles_counts, trav_stats['data'], stats_array[0,:]= mod.count(d_tree_1, d_tree_2, shell, self_corr, v.strategy)

    logger.info("Particles counts:")
    logger.info("DD %s", particles_counts)
    result_file.write("Particles_counts ")
    np.savetxt(result_file, particles_counts[np.newaxis,:])
    result_file.flush()

    # Random counting.
    logger.info("Starting random-random counts...")

    random_counts, trav_stats['random'], stats_array[1,:] = mod.count(r_tree_1, r_tree_2, shell, self_corr, v.strategy)

    logger.info("Random counts:")
    logger.info("RR %s", random_counts)
    result_file.write("Random_counts ")
    np.savetxt(result_file, random_counts[np.newaxis,:])
    result_file.flush()
    
    # Mixed counting 1 (needed for the LS estimator).
    logger.info("Starting mixed_1 counts...")

    mixed_counts_1, trav_stats['mixed 1'], stats_array[2,:] = mod.count(d_tree_1, r_tree_2, shell, False, v.strategy)

    logger.info("Mixed_count_1:")
    logger.info("M1 %s", mixed_counts_1)
    result_file.write("Mixed_counts_1 ")
    np.savetxt(result_file, mixed_counts_1[np.newaxis,:])
    result_file.flush()

    # Mixed counting 2 (needed for the LS estimator).    
    logger.info("Starting mixed_2 counts...")

    mixed_counts_2, trav_stats['mixed 2'], stats_array[3,:] = mod.count(r_tree_1, d_tree_2, shell, False, v.strategy)

    logger.info("Mixed_count_2:")
    logger.info("M2 %s", mixed_counts_2)
    result_file.write("Mixed_counts_2 ")
    np.savetxt(result_file, mixed_counts_2[np.newaxis,:])
    result_file.flush()
    
    # Correlate.
    #TPCF, TPCF_err, stats['Correlation time'] =  mod.correlate(particles_counts, mixed_counts_1, mixed_counts_2, random_counts)
    
    # Print statistics.
    logger.info("particles_counts %s", particles_counts)
    logger.info("random_counts %s", random_counts)
    logger.info("mixed_counts_1 %s", mixed_counts_1)
    logger.info("mixed_counts_2 %s", mixed_counts_2)
    #logger.info( "TPCF %s", TPCF)
    #logger.info("r %s", r)

    #result_file.write("TPCF ")
    #np.savetxt(result_file, TPCF[np.newaxis,:])
    result_file.flush()
    result_file.close()

    # Plot results.
    #mod.plot_correlation(r, TPCF, TPCF_err)

#===============================================================================
#    Info
#===============================================================================

    # Print some other statistics.
    logger.info("########## Stats ##########")
    logger.info("Strategy %s", v.strategy)
    logger.info("Set 1 dimensions %s",  stats['Set 1'])
    logger.info("Set 2 dimensions %s", stats['Set 2'])
    logger.info("Set 1 mean density %s",  (stats['Set 1']/
                                           (b_size_1[0]*b_size_1[1]*b_size_1[2])))
    logger.info("Set 2 mean density %s", (stats['Set 2']/
                                           (b_size_2[0]*b_size_2[1]*b_size_2[2])))
    logger.info("Leafsize %s", v.leafsize)
    logger.info('d_tree_1 build time %s', stats['d_tree_1 build time'])
    logger.info('d_tree_2 build time %s', stats['d_tree_2 build time'])
    logger.info('r_tree_1 build time %s', stats['r_tree_1 build time'])
    logger.info('r_tree_2 build time %s', stats['r_tree_2 build time'])    
    logger.info('Time elapsed before traversing %s', stats['Time elapsed before traversing'])
    logger.info("Data traverse:")
    logger.info('Time for traverse %s', trav_stats['data']['Time for traverse'])  
    logger.info('Number of traverse %s', trav_stats['data']['Number of traverse']) 
    logger.info('Expected traverse %s', trav_stats['data']['Expected traverse'])
    logger.info('Traverse speed %s', trav_stats['data']['Traverse speed'])     
    logger.info('n_opened_leaves %s', trav_stats['data']['n_opened_leaves'])
    logger.info('sort_times %s', trav_stats['data']['sort_times'])
    logger.info('dist_times %s', trav_stats['data']['dist_times'])
    logger.info( "Wall time for all %s", time.time()-tt_glob)
    logger.info("Writing %s", v.log_file+".dat")
    mod.print_stats_array(v.log_file, stats_array)

#===============================================================================
#    Quitting
#===============================================================================

    # Save proc status file.
    logger.info("Copying proc status file at the end.")
    status_src = '/proc/'+str(os.getpid())+'/status'
    status_dest = '../logs/status'+str(os.getpid())+'@end'
    shcopy(status_src, status_dest)

    # Closing.
    logger.info("That's all Folks!:P")
    logging.shutdown()
    return 42 #The Answer.


#===============================================================================
#    Main start
#===============================================================================

# Start main function.
if __name__ == "__main__":
    sys.exit(main())
