#! /use/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker 

class MyFormatter(ticker.LogFormatter): 
    """Non funziona... servirebbe a sistemare i label degli assi logaritmici."""
    def __call__(self, val, pos=None): 
        fx = int(np.floor(np.log(abs(val))/np.log(self._base) +0.5))  
        isDecade = self.is_decade(fx) 
        if not isDecade and self.labelOnlyBase: 
            return '' 
        if (fx%2)==1:
            return '' 
        return '%d'%fx 

def plot_kd3_tests(trav_file="prova_traverse.dat",
                   dist_file="prova_distanze.dat",
                   plot='all',
                   n_part_dim=8,
                   n_chunks=7,
                   delimiter=None,
                   tot_part_array=np.array([1, 2, 4, 8, 16, 32, 64, 128])*100000,
                   xscale='log', yscale='log', machine_name='Uno',
                   plot_name='traverse_vs_distances_yxlog'):
    """Plot the kd3 tests.
    """
 
    class MyFormatter(ticker.LogFormatter): 
        """Non funziona... servirebbe a sistemare i label degli assi logaritmici."""
        def __call__(self, val, pos=None): 
            fx = int(np.floor(np.log(abs(val))/np.log(self._base) +0.5))  
            isDecade = self.is_decade(fx) 
            if not isDecade and self.labelOnlyBase: 
                return '' 
            if (fx%2)==1:
                return '' 
            return '%d'%fx 

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as font_manager
    import matplotlib.ticker as ticker
 
    t_traverse = np.genfromtxt(trav_file, dtype='f8', comments='#', delimiter=delimiter, usecols=(4))
    t_traverse = t_traverse.reshape((n_part_dim,n_chunks))
    t_dist = np.genfromtxt(dist_file, dtype='f8', comments='#', delimiter=delimiter, usecols=(1))
    size = np.genfromtxt(dist_file, dtype='i4', comments='#', delimiter=delimiter, usecols=(0))
    tot_part = tot_part_array
    
    t_dist_tot = np.empty(t_traverse.shape)
    for i in range(t_traverse.shape[0]):
	t_dist_tot[i,:] = t_dist[:]*tot_part[i]
        
    somma = t_traverse + t_dist_tot
        
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_title('Confronto tempi per distanze fortran e traverse fortran @'+machine_name)
    ax.set_xlabel('Dimensione foglia')
    ax.set_ylabel('Tempo [s]')
    ax.set_xticks(range(size.size))
    ax.set_xticklabels(size)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.grid(True)#, which="both")
    #plt.grid(True,which="both",ls="-", color='0.65')
    colors = ['blue', 'red', 'green', 'black', 'magenta', 'yellow']

    for i in range(t_traverse.shape[0]):
        ax.plot(t_traverse[i,:], color = colors[i%6], linestyle = ':', marker = 'x', label = 'traverse '+str(tot_part[i]))
    
    for i in range(t_traverse.shape[0]):
        ax.plot(t_dist_tot[i,:], color = colors[i%6], linestyle = '--', marker = 'x', label = 'dist '+str(tot_part[i]))
    
    for i in range(t_traverse.shape[0]):
        ax.plot(somma[i,:], color = colors[i%6], linestyle = '-', marker = 'x', label = 'somma '+str(tot_part[i]))
    
    ax.legend(loc='best', prop=font_manager.FontProperties(size=7))
    plt.savefig(plot_name)



def old_1():
    tempi = np.genfromtxt('tempi.dat', dtype='float', comments='#', usecols=(0, 2, 3, 4, 5, 6, 7))
    
    ### Tempo totale per tutti i metodi per il numero di particelle
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    all_y = tempi[[11, 23, 35, 47], 1:7]
    ax.set_title('Tempi totali di esecuzione per i 6 metodi')
    ax.set_xlabel('Numero di particelle')
    ax.set_ylabel('Tempo [s]')
    ax.set_xticks(range(4))
    ax.set_xticklabels(['1000', '3200', '10000', '32000'])
    ax.plot(all_y[:,0], color = 'blue', linestyle = '-', marker = 'x', label = 'lin_nosqrt_sort')
    ax.plot(all_y[:,1], color = 'cyan', linestyle = '-', marker = 'x', label = 'lin_sqrt_nosort')
    ax.plot(all_y[:,2], color = 'green', linestyle = '-', marker = 'x', label = 'lin_sqrt_sort')
    ax.plot(all_y[:,3], color = 'black', linestyle = '-', marker = 'x', label = 'log_nosqrt_nosort')
    ax.plot(all_y[:,4], color = 'red', linestyle = '-', marker = 'x', label = 'log_sqrt_nosort')
    ax.plot(all_y[:,5], color = 'yellow', linestyle = '-', marker = 'x', label = 'log_sqrt_sort')
    ax.legend(loc='best')
    ax.grid(True)
    plt.savefig('uno')

    ### Tempo per metodo per tutte le funzioni
    fig2 = plt.figure()
    ax1 = fig2.add_subplot(231)
    ax1.set_title('Timing di lin_nosqrt_sort')
    ax1.set_xlabel('Funzione')
    ax1.set_ylabel('Tempo [s]')
    ax1.set_ylim(ymax=6000)
    ax1.set_xticks(range(12))
    ax1.set_xticklabels(['sh_file', 'pos_file', 'sets', 'random_c', 'tree1', 'tree2', 'count', 'tree1_r', 'tree2_r', 'count_r', 'corr', 'all'], rotation = '50')
    ax1.plot(tempi[0:12, 1], color = 'blue', linestyle = '-', marker = 'x', label = '1000 particelle')
    ax1.plot(tempi[12:24, 1], color = 'green', linestyle = '-', marker = 'x', label = '3200 particelle')
    ax1.plot(tempi[24:36, 1], color = 'black', linestyle = '-', marker = 'x', label = '10000 particelle')
    ax1.plot(tempi[36:48, 1], color = 'red', linestyle = '-', marker = 'x', label = '32000 particelle')
    ax1.legend(loc='best')
    ax1.grid(True)
    
    ax2 = fig2.add_subplot(232)
    ax2.set_title('Timing di lin_sqrt_nosort')
    ax2.set_xlabel('Funzione')
    ax2.set_ylabel('Tempo [s]')
    ax2.set_ylim(ymax=6000)
    ax2.set_xticks(range(12))
    ax2.set_xticklabels(['sh_file', 'pos_file', 'sets', 'random_c', 'tree1', 'tree2', 'count', 'tree1_r', 'tree2_r', 'count_r', 'corr', 'all'], rotation = '50')
    ax2.plot(tempi[0:12, 2], color = 'blue', linestyle = '-', marker = 'x', label = '1000 particelle')
    ax2.plot(tempi[12:24, 2], color = 'green', linestyle = '-', marker = 'x', label = '3200 particelle')
    ax2.plot(tempi[24:36, 2], color = 'black', linestyle = '-', marker = 'x', label = '10000 particelle')
    ax2.plot(tempi[36:48, 2], color = 'red', linestyle = '-', marker = 'x', label = '32000 particelle')
    ax2.legend(loc='best')
    ax2.grid(True)
    
    ax3 = fig2.add_subplot(233)
    ax3.set_title('Timing di lin_sqrt_sort')
    ax3.set_xlabel('Funzione')
    ax3.set_ylabel('Tempo [s]')
    ax3.set_xticks(range(12))
    ax3.set_ylim(ymax=6000)
    ax3.set_xticklabels(['sh_file', 'pos_file', 'sets', 'random_c', 'tree1', 'tree2', 'count', 'tree1_r', 'tree2_r', 'count_r', 'corr', 'all'], rotation = '50')
    ax3.plot(tempi[0:12, 3], color = 'blue', linestyle = '-', marker = 'x', label = '1000 particelle')
    ax3.plot(tempi[12:24, 3], color = 'green', linestyle = '-', marker = 'x', label = '3200 particelle')
    ax3.plot(tempi[24:36, 3], color = 'black', linestyle = '-', marker = 'x', label = '10000 particelle')
    ax3.plot(tempi[36:48, 3], color = 'red', linestyle = '-', marker = 'x', label = '32000 particelle')
    ax3.legend(loc='best')
    ax3.grid(True)
    
    ax4 = fig2.add_subplot(234)
    ax4.set_title('Timing di log_nosqrt_nosort')
    ax4.set_xlabel('Funzione')
    ax4.set_ylabel('Tempo [s]')
    ax4.set_xticks(range(12))
    ax4.set_ylim(ymax=6000)
    ax4.set_xticklabels(['sh_file', 'pos_file', 'sets', 'random_c', 'tree1', 'tree2', 'count', 'tree1_r', 'tree2_r', 'count_r', 'corr', 'all'], rotation = '50')
    ax4.plot(tempi[0:12, 4], color = 'blue', linestyle = '-', marker = 'x', label = '1000 particelle')
    ax4.plot(tempi[12:24, 4], color = 'green', linestyle = '-', marker = 'x', label = '3200 particelle')
    ax4.plot(tempi[24:36, 4], color = 'black', linestyle = '-', marker = 'x', label = '10000 particelle')
    ax4.plot(tempi[36:48, 4], color = 'red', linestyle = '-', marker = 'x', label = '32000 particelle')
    ax4.legend(loc='best')
    ax4.grid(True)
    
    ax5 = fig2.add_subplot(235)
    ax5.set_title('Timing di log_sqrt_nosort')
    ax5.set_xlabel('Funzione')
    ax5.set_ylabel('Tempo [s]')
    ax5.set_xticks(range(12))
    ax5.set_ylim(ymax=6000)
    ax5.set_xticklabels(['sh_file', 'pos_file', 'sets', 'random_c', 'tree1', 'tree2', 'count', 'tree1_r', 'tree2_r', 'count_r', 'corr', 'all'], rotation = '50')
    ax5.plot(tempi[0:12, 5], color = 'blue', linestyle = '-', marker = 'x', label = '1000 particelle')
    ax5.plot(tempi[12:24, 5], color = 'green', linestyle = '-', marker = 'x', label = '3200 particelle')
    ax5.plot(tempi[24:36, 5], color = 'black', linestyle = '-', marker = 'x', label = '10000 particelle')
    ax5.plot(tempi[36:48, 5], color = 'red', linestyle = '-', marker = 'x', label = '32000 particelle')
    ax5.legend(loc='best')
    ax5.grid(True)
    
    ax6 = fig2.add_subplot(236)
    ax6.set_title('Timing di log_sqrt_sort')
    ax6.set_xlabel('Funzione')
    ax6.set_ylabel('Tempo [s]')
    ax6.set_xticks(range(12))
    ax6.set_ylim(ymax=6000)
    ax6.set_xticklabels(['sh_file', 'pos_file', 'sets', 'random_c', 'tree1', 'tree2', 'count', 'tree1_r', 'tree2_r', 'count_r', 'corr', 'all'], rotation = '50')
    ax6.plot(tempi[0:12, 6], color = 'blue', linestyle = '-', marker = 'x', label = '1000 particelle')
    ax6.plot(tempi[12:24, 6], color = 'green', linestyle = '', marker = 'x', label = '3200 particelle')
    ax6.plot(tempi[24:36, 6], color = 'black', linestyle = '-', marker = 'x', label = '10000 particelle')
    ax6.plot(tempi[36:48, 6], color = 'red', linestyle = '-', marker = 'x', label = '32000 particelle')
    ax6.legend(loc='best')
    ax6.grid(True)
    
    fig2.set_size_inches(25, 12)
    
    plt.savefig('due')
    
def old_2():
    ## Foglie

    foglie = np.genfromtxt('foglie_old.dat', dtype='float', comments='#', usecols=(0, 1, 2))

    fig3 = plt.figure()
    ax7 = fig3.add_subplot(311)
    ax7.set_title('Tempi di esecuzione per dimensione delle foglie')
    ax7.set_xlabel('Dimensione foglie')
    ax7.set_ylabel('Tempo [s]')
    leafnumbers = ['50', '100', '150', '200', '250', '300', '325', '350', '375', '400', '425', '450', '475', '500', '525', '550', '575', '600', '625', '650', '675', '700', '750', '800']
    ax7.set_xticks(range(24))
    ax7.set_xticklabels(leafnumbers)
    ax7.plot(foglie[0:24,2], color = 'blue', linestyle = '-', marker = 'x', label = '1000 particelle')
    ax7.legend(loc='best')
    ax7.grid(True)
    
    ax8 = fig3.add_subplot(312)
    ax8.set_title('Tempi di esecuzione per dimensione delle foglie')
    ax8.set_xlabel('Dimensione foglie')
    ax8.set_ylabel('Tempo [s]')
    leafnumbers = ['50', '100', '150', '200', '250', '300', '325', '350', '375', '400', '425', '450', '475', '500', '525', '550', '575', '600', '625', '650', '675', '700', '750', '800']
    ax8.set_xticks(range(24))
    ax8.set_xticklabels(leafnumbers)
    ax8.plot(foglie[24:48,2], color = 'cyan', linestyle = '-', marker = 'x', label = '10000 particelle')
    ax8.legend(loc='best')
    ax8.grid(True)
    
    ax9 = fig3.add_subplot(313)
    ax9.set_title('Tempi di esecuzione per dimensione delle foglie')
    ax9.set_xlabel('Dimensione foglie')
    ax9.set_ylabel('Tempo [s]')
    leafnumbers = ['50', '100', '150', '200', '250', '300', '325', '350', '375', '400', '425', '450', '475', '500', '525', '550', '575', '600', '625', '650', '675', '700', '750', '800']
    ax9.set_xticks(range(24))
    ax9.set_xticklabels(leafnumbers)
    ax9.plot(foglie[48:62,2], color = 'green', linestyle = '-', marker = 'x', label = '100000 particelle')
    ax9.legend(loc='best')
    ax9.grid(True)
    
    fig2.set_size_inches(15, 20)
    
    plt.savefig('tre')



def old_trav_dist_1():

    ###############plot traverse e distanze

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    traverse_10000 = np.array([0.808617115021, 0.24688911438, 0.209056854248, 0.0821368694305, 0.0643391609192, 0.032219171524, 0.0151379108429, 0.00276589393616, 0.00111198425293, 0.000225782394409])
    traverse_100000 = np.array([56.9056639671, 10.7256379128, 6.06676197052, 4.58031392097, 3.08165502548, 1.40718317032, 0.30656003952, 0.0901081562042, 0.0680449008942, 0.0303068161011])
    traverse_1000000 = np.array([0, 0, 1305.386163, 848.837960005, 582.702934027, 192.339005947, 61.6228399277, 13.0040619373, 5.46887397766, 4.03734183311])
    size = np.array([100, 200, 300, 400, 500, 1000, 2000, 5000, 8000, 10000])
    
    dist = np.array([2.73203850e-03, 5.65719604e-03, 1.28698349e-03, 1.94501877e-03, 8.79216194e-03, 3.45869064e-02, 2.17735052e-01, 5.59701920e-01, 2.01534510e+00, 7.04916716e+00])
    dist_10000 = dist*(10000/size)
    dist_100000 = dist*(100000/size)
    dist_1000000 = dist*(1000000/size)
    
    #dist_2 = np.array([5.29909134e-03, 9.28592682e-03, 1.89995766e-03, 8.25595856e-03, 3.27820778e-02, 2.05492020e-01, 5.29361963e-01, 8.24096918e-01])
    #dist2_10000 = dist_2*(10000/size)
    #dist2_100000 = dist_2*(100000/size)
    #dist2_1000000 = dist_2*(1000000/size)
    
    somma_10000 = traverse_10000 + dist_10000
    somma_100000 = traverse_100000 + dist_100000
    somma_1000000 = traverse_1000000 + dist_1000000
    
    ax.set_title('Confronto calcolo della distanza e traverse')
    ax.set_xlabel('Dimensione foglia')
    ax.set_ylabel('Tempo [s]')
    ax.set_xticks(range(size.size))
    ax.set_xticklabels(size)
    
    ax.plot(traverse_10000, color = 'blue', linestyle = ':', marker = 'x', label = 'traverse 10000')
    ax.plot(traverse_100000, color = 'red', linestyle = ':', marker = 'x', label = 'traverse 100000')
    ax.plot(traverse_1000000, color = 'green', linestyle = ':', marker = 'x', label = 'traverse 1000000')
    ax.plot(dist_10000, color = 'blue', linestyle = '--', marker = 'x', label = 'estimated dist 10000')
    ax.plot(dist_100000, color = 'red', linestyle = '--', marker = 'x', label = 'estimated dist 100000')
    ax.plot(dist_1000000, color = 'green', linestyle = '--', marker = 'x', label = 'estimated dist 1000000')
    #ax.plot(dist2_10000, color = 'blue', linestyle = '--', marker = 'x', label = 'estimated dist 10000')
    #ax.plot(dist2_100000, color = 'red', linestyle = '--', marker = 'x', label = 'estimated dist 100000')
    #ax.plot(dist2_1000000, color = 'green', linestyle = '--', marker = 'x', label = 'estimated dist 1000000')
    ax.plot(somma_10000, color = 'blue', linestyle = '-', marker = 'x', label = 'somma 10000')
    ax.plot(somma_100000, color = 'red', linestyle = '-', marker = 'x', label = 'somma 100000')
    ax.plot(somma_1000000, color = 'green', linestyle = '-', marker = 'x', label = 'somma 1000000')

    ax.legend(loc='best', prop=font_manager.FontProperties(size=7))
    ax.grid(True)
    plt.savefig('traverse_vs_distances')
    ax.set_yscale('log')
    plt.savefig('traverse_vs_distances_ylog')
    ax.set_xscale('log')
    #ax.yaxis.set_major_formatter(MyFormatter()) 
    plt.savefig('traverse_vs_distances_yxlog')
    ax.set_yscale('linear')

def old_trav_dist_2():
    # Dipendenza dei traverse fatti dalle dimensioni della foglia e dalle particelle totali

    traverse = np.genfromtxt('prova_traverse.dat', dtype='f8', comments='#', delimiter=None, usecols=(2,3))
    traverse_ratio = traverse[:,0]*1./traverse[:,1]
    #tot_part = np.genfromtxt('prova_traverse.dat', dtype='f8', comments='#', delimiter=None, usecols=(0))
    #leafsize = np.genfromtxt('prova_traverse.dat', dtype='f8', comments='#', delimiter=None, usecols=(1))
    fig1=plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_title('Ratio of done traverse')
    ax.set_xlabel('leafsize')
    #leafsize = np.array(['300', '600', '1200', '2400', '4800', '9600', '19200'])
    size = np.array([300, 600, 1200, 2400, 4800, 9600, 19200])
    ax.set_xticks(range(size.size))
    ax.set_xticklabels(size)
    #ax.set_xticks(range(leafsize.size))
    #ax.set_xticklabels=(leafsize)
    ax.set_ylabel('done/expected')
    ax.plot(traverse_ratio[0:7], color = 'blue', linestyle = '-', marker = 'x', label = '100000')
    ax.plot(traverse_ratio[7:14], color = 'red', linestyle = '-', marker = 'x', label = '200000')
    ax.plot(traverse_ratio[14:21], color = 'green', linestyle = '-', marker = 'x', label = '400000')
    ax.plot(traverse_ratio[21:28], color = 'black', linestyle = '-', marker = 'x', label = '800000')
    ax.plot(traverse_ratio[28:35], color = 'blue', linestyle = ':', marker = 'x', label = '1600000')
    ax.plot(traverse_ratio[35:42], color = 'red', linestyle = ':', marker = 'x', label = '3200000')
    ax.plot(traverse_ratio[42:49], color = 'green', linestyle = ':', marker = 'x', label = '6400000')
    ax.plot(traverse_ratio[49:56], color = 'black', linestyle = ':', marker = 'x', label = '12800000')
    ax.legend(loc='best', prop=font_manager.FontProperties(size=7))
    ax.grid(True)
    plt.show()
    plt.savefig('ratio_of_done_traverse')
    ax.set_yscale('log')
    plt.savefig('ratio_of_done_traverse_ylog')
    ax.set_xscale('log')
    plt.savefig('ratio_of_done_traverse_yxlog')

def old_trav_dist_3():
    #plot di traverse e distanze su uno col programma finito con fortran su tutte le dist

    t_traverse = np.genfromtxt('prova_traverse_fortran.dat', dtype='f8', comments='#', delimiter=None, usecols=(4))
    t_traverse = t_traverse.reshape((8,7))
    t_dist = np.genfromtxt('prova_distanze.dat', dtype='f8', comments='#', delimiter=None, usecols=(1))
    size = np.genfromtxt('prova_distanze.dat', dtype='i4', comments='#', delimiter=None, usecols=(0))
    tot_part = np.array([1, 2, 4, 8, 16, 32, 64, 128])*100000
    
    t_dist_tot = np.empty(t_traverse.shape)
    for i in range(t_traverse.shape[0]):
	t_dist_tot[i,:] = t_dist[:]*tot_part[i]
        
    somma = t_traverse + t_dist_tot
        
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.set_title('Confronto tempi per distanze fortran e traverse fortran @Uno')
    ax.set_xlabel('Dimensione foglia')
    ax.set_ylabel('Tempo [s]')
    ax.set_xticks(range(size.size))
    ax.set_xticklabels(size)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True)#, which="both")
    #plt.grid(True,which="both",ls="-", color='0.65')
    
    ax.plot(t_traverse[0,:], color = 'blue', linestyle = ':', marker = 'x', label = 'traverse 100000')
    ax.plot(t_traverse[1,:], color = 'red', linestyle = ':', marker = 'x', label = 'traverse 200000')
    ax.plot(t_traverse[2,:], color = 'green', linestyle = ':', marker = 'x', label = 'traverse 400000')
    ax.plot(t_traverse[3,:], color = 'cyan', linestyle = ':', marker = 'x', label = 'traverse 800000')
    ax.plot(t_traverse[4,:], color = 'black', linestyle = ':', marker = 'x', label = 'traverse 1600000')
    ax.plot(t_traverse[5,:], color = 'magenta', linestyle = ':', marker = 'x', label = 'traverse 3200000')
    ax.plot(t_traverse[6,:], color = 'yellow', linestyle = ':', marker = 'x', label = 'traverse 6400000')
    ax.plot(t_traverse[7,:], color = 'blue', linestyle = ':', marker = 'x', label = 'traverse 12800000')
    
    ax.plot(t_dist_tot[0,:], color = 'blue', linestyle = '--', marker = 'x', label = 'dist 100000')
    ax.plot(t_dist_tot[1,:], color = 'red', linestyle = '--', marker = 'x', label = 'dist 200000')
    ax.plot(t_dist_tot[2,:], color = 'green', linestyle = '--', marker = 'x', label = 'dist 400000')
    ax.plot(t_dist_tot[3,:], color = 'cyan', linestyle = '--', marker = 'x', label = 'dist 800000')
    ax.plot(t_dist_tot[4,:], color = 'black', linestyle = '--', marker = 'x', label = 'dist 1600000')
    ax.plot(t_dist_tot[5,:], color = 'magenta', linestyle = '--', marker = 'x', label = 'dist 3200000')
    ax.plot(t_dist_tot[6,:], color = 'yellow', linestyle = '--', marker = 'x', label = 'dist 6400000')
    ax.plot(t_dist_tot[7,:], color = 'blue', linestyle = '--', marker = 'x', label = 'dist 12800000')
    
    ax.plot(somma[0,:], color = 'blue', linestyle = '-', marker = 'x', label = 'somma 100000')
    ax.plot(somma[1,:], color = 'red', linestyle = '-', marker = 'x', label = 'somma 200000')
    ax.plot(somma[2,:], color = 'green', linestyle = '-', marker = 'x', label = 'somma 400000')
    ax.plot(somma[3,:], color = 'cyan', linestyle = '-', marker = 'x', label = 'somma 800000')
    ax.plot(somma[4,:], color = 'black', linestyle = '-', marker = 'x', label = 'somma 1600000')
    ax.plot(somma[5,:], color = 'magenta', linestyle = '-', marker = 'x', label = 'somma 3200000')
    ax.plot(somma[6,:], color = 'yellow', linestyle = '-', marker = 'x', label = 'somma 6400000')
    ax.plot(somma[7,:], color = 'blue', linestyle = '-', marker = 'x', label = 'somma 12800000')
    
    ax.legend(loc='best', prop=font_manager.FontProperties(size=7))
    
    plt.show()
    
    plt.savefig('traverse_vs_fortran_distances_yxlog')

