import numpy as np
import matplotlib
matplotlib.use('Agg') # define display on the server
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

# How many times rand are greater than data.
mult_fact = 1
description = "gif2_rev52_01_1_random_10x"

data = np.genfromtxt('particles.dat')
rand = np.genfromtxt('random.dat')
m1 = np.genfromtxt('m1.dat')
m2 = np.genfromtxt('m2.dat')
r = np.genfromtxt('radii.dat')

if data.ndim != 1:
	data = data[:, 1:-1]
	rand = rand[:, 1:-1] / mult_fact**2
	m1 = m1[:, 1:-1] / mult_fact
	m2 = m2[:, 1:-1] / mult_fact
	r = r[0, 1:]
	print "summing data"
	print data.shape
	data = data.sum(0)
	print data.shape
	rand = rand.sum(0)
	m1 = m1.sum(0)
	m2 = m2.sum(0)
else:
	print "not summing"
	data = data[1:-1]
	rand = rand[1:-1] / mult_fact**2
	m1 = m1[1:-1] / mult_fact
	m2 = m2[1:-1] / mult_fact
	r = r[0,1:]

np.savetxt('r.dat', r)

print data.shape, r.shape

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Raw counts')
ax.set_xlabel('Radius [Mpc/h]')
ax.set_ylabel('Counts')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.plot(r, data, color = 'blue', linestyle = '-', marker = 'x', label = 'Data')
ax.plot(r, rand, color = 'black', linestyle = '-', marker = 'x', label = 'Random')
ax.plot(r, m1, color = 'green', linestyle = '-', marker = 'x', label = 'Mixed_1')
ax.plot(r, m2, color = 'red', linestyle = '-', marker = 'x', label = 'Mixed_2')
ax.legend(loc='best', prop=font_manager.FontProperties(size=7))
#plt.show()
plt.savefig('conteggi_raw_'+description)

TPCF = (data-m1-m2+rand)/(rand)
TPCF_err = np.sqrt(data*(1+data/rand))/rand
np.savetxt('TPCF.dat', TPCF)
np.savetxt('TPCF_err.dat', TPCF_err)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('LS two point correlation function')
ax.set_xlabel('Radius [Mpc/h]')
ax.set_ylabel('TPCF')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.errorbar(r, TPCF, TPCF_err, color = 'blue', linestyle = '-', marker = '.', label = 'TPCF')
ax.legend(loc='best', prop=font_manager.FontProperties(size=12))
#plt.show()
plt.savefig('TPCF_'+description)



#####################################
#da sistemare.. mancano gli ultimi bin

def bin(n_bin, data, rand, m1, m2):
	n_bin = 9
        #quanti ne metto assieme, meglio dispari cos\`i tengo il bin centrale come media
	n_tog = r.size/n_bin        

	r_reb = np.zeros(n_bin)
	p_reb = np.zeros(n_bin)
	rand_reb = np.zeros(n_bin)
	m1_reb = np.zeros(n_bin)
	m2_reb = np.zeros(n_bin)
	
	central_bin=3#contando da zero, quindi human readable sarebbe 1 in pi\`u

	for i in range(n_bin):
		r_reb[i] = r[i*n_tog+central_bin]
		print i*n_tog+central_bin
		data_reb[i] = data[i*n_tog:i*n_tog+n_tog-1].sum()
		m1_reb[i] = m1[i*n_tog:i*n_tog+n_tog-1].sum()
		m2_reb[i] = m2[i*n_tog:i*n_tog+n_tog-1].sum()
		rand_reb[i] = rand_div[i*n_tog:i*n_tog+n_tog-1].sum()

	return data_reb, rand_reb, m1_reb, m2_reb, r_reb
