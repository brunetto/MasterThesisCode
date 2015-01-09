import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

# How many times rand are greater than data.
mult_fact = 8

data = np.genfromtxt('particles.dat')[1:-1]
rand = np.genfromtxt('random.dat')[1:-1]
m1 = np.genfromtxt('m1.dat')[1:-1]
m2 = np.genfromtxt('m2.dat')[1:-1]
r = np.genfromtxt('radii.dat')[1:]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Raw counts')
ax.set_xlabel('Radius [Mpc/h]')
ax.set_ylabel('Counts')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.plot(r, data, color = 'blue', linestyle = '-', marker = 'x', label = 'Data1')
ax.plot(r, rand, color = 'black', linestyle = '-', marker = 'x', label = 'Random1')
ax.plot(r, m1, color = 'green', linestyle = '-', marker = 'x', label = 'Mixed_11')
ax.plot(r, m2, color = 'red', linestyle = '-', marker = 'x', label = 'Mixed_21')
ax.legend(loc='best', prop=font_manager.FontProperties(size=7))
plt.show()
plt.savefig('conteggi_raw_millimil_rev44')


TPCF = (data-m1-m2+rand)/(rand)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('LS two point correlation function')
ax.set_xlabel('Radius [Mpc/h]')
ax.set_ylabel('TPCF')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.plot(r, TPCF, color = 'blue', linestyle = '-', marker = 'x', label = 'TPCF')
ax.legend(loc='best', prop=font_manager.FontProperties(size=12))
plt.show()
plt.savefig('TPCF_millimil_confronto')



#####################################
#da sistemare.. mancano gli ultimi bin

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
	p_reb[i] = data[i*n_tog:i*n_tog+n_tog-1].sum()
	m1_reb[i] = m1[i*n_tog:i*n_tog+n_tog-1].sum()
	m2_reb[i] = m2[i*n_tog:i*n_tog+n_tog-1].sum()
	rand_reb[i] = rand_div[i*n_tog:i*n_tog+n_tog-1].sum()

TPCF = (p_reb-m1_reb-m2_reb+rand_reb)/(rand_reb)
TPCF_err = np.sqrt(TPCF)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Raw counts rebinned')
ax.set_xlabel('Radius [Mpc/h]')
ax.set_ylabel('Counts')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.plot(r_reb, p_reb, color = 'blue', linestyle = '-', marker = 'x', label = 'Data')
ax.plot(r_reb, rand_reb, color = 'black', linestyle = '-', marker = 'x', label = 'Random')
ax.plot(r_reb, m1_reb, color = 'green', linestyle = '-', marker = 'x', label = 'Mixed_1')
ax.plot(r_reb, m2_reb, color = 'red', linestyle = '-', marker = 'x', label = 'Mixed_2')
ax.legend(loc='best', prop=font_manager.FontProperties(size=12))
plt.show()
plt.savefig('conteggi_raw_rebinnati_con errori')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('LS two point correlation function rebinned')
ax.set_xlabel('Radius [Mpc/h]')
ax.set_ylabel('TPCF')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(True)
ax.errorbar(r_reb, TPCF, TPCF_err, color = 'blue', linestyle = '-', marker = 'x', label = 'TPCF')
ax.legend(loc='best', prop=font_manager.FontProperties(size=12))
plt.show()
plt.savefig('TPCF_rebinned')
