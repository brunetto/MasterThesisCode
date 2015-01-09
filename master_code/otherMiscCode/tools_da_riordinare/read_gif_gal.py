#! /usr/bin/env python

file = open('lcdm_galaxy_cat.z0.00', 'rb')
i = 0
for riga in file.readlines():
    parole = riga.split()
    if len(parole) == 11:
        print "iterazione ", i 
        print "aggiungo "
        x.append(parole[5])
        y.append(parole[6])
        z.append(parole[7])
    i+=1
