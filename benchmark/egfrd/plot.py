#!/usr/bin/env/python

import sys

import numpy
import imp
import scipy.io
from matplotlib.pylab import *


N_A = 6.0221367e23

def plot_data(N, T, fmt):
    T = numpy.array(T)

    mean = T.mean(1)
    std_err = T.std()/math.sqrt(len(T))

    #errorbar(N, mean, yerr=std_err, fmt=fmt)
    print N, mean
    loglog(N, mean, fmt)


imp.load_source('out_Vspatiocyte', 'out_Vspatiocyte.py')
imp.load_source('out_VspatiocyteOld', 'out_VspatiocyteOld.py')
imp.load_source('out_VspatiocytePrevious', 'out_VspatiocytePrevious.py')
from out_Vspatiocyte import *
from out_VspatiocyteOld import *
from out_VspatiocytePrevious import *
from out_C import *
from out_V import *
from out_N300 import *
from out_N3000 import *
from out_BD import *
#from out_BD2 import *

from run_all import Nv, Nc, V300, N300, V3000, N3000, VBD, NBD, TBD, TBD2


# (40e-18 ** (1/3.0))**2 / 1e-12
# = 11.69607095285148

X = numpy.array([30,100,300,1000,3000,10000,30000,100000,1e8])


axes([.12,.14,.86,.83])

#for i in range(len(Nv)):
plot_data(Nv, data_V,'kx')

loglog(X, 0.0004* X**(5.0/3), 'k--')

figtext(.25, .18, r'(2) V = 1 pL')
figtext(.82, .85, r'$t \ \propto \ N^{5/3}$', color='k')


#for i in range(len(Nc)):
plot_data(Nc, data_C,'kx')
loglog(X, 0.19*X, 'b-')


#plot_data(Nv, data_Vspatiocyte,'r.')
plot_data(Nv, data_VspatiocyteOld,'k.')
plot_data(Nv, data_VspatiocytePrevious,'k.')

figtext(.14, .4, r'(1) C = 50 nM')
figtext(.8, .59, r'$t \  \propto \ N$', color='k')

# plot BD data
#plot_data(NBD, data_BD2,'k.')
plot_data(NBD, data_BD,'k.')

# loglog(X, 2e6* X, 'b:') # 1e-6 tau
# loglog(X, 2e4* X, 'b:') # 1e-4 tau
loglog(X, 1.3e5* X, 'b:') # 1e-5 tau

#figtext(.2, .82, r'BD', color='k')
#figtext(.19, .64, r'BD (relaxed)', color='k')

figtext(.19, .64, r'BD', color='k')

#loglog(data1[0] , data1[1], 'o-', label='Vol. = 1e-15 L')
#loglog(data2[0] , data2[1], 'o-', label='# particles = 600')
#loglog(data3[0] , data3[1], 'o-', label='Conc. = 1e-6 M')

xlabel('N [# particles]', size=22)
#xlabel('Concentration [M]')
#ylabel('time [s]', size=22)
#legend()
# xlim(4,9e6)
# ylim(1.1,2e11)

xlim(40,9e7)
ylim(0,2e11)

xticks(size=18)
yticks([60,3600,3600*24,3600*24*30, 3600*24*30*12],
       ['minute', 'hour', 'day', 'month', 'year'], size=16)

#grid()

C300 = numpy.array(N300) / (numpy.array(V300)*N_A)
C3000 = numpy.array(N3000) / (numpy.array(V3000)*N_A)
print C300, C3000
# Cx3000=numpy.array([
# #    9.35e-11,
#     9.35e-10,
#     9.35e-9,
#     9.35e-8,#N=3000,V=40e-15
#     9.35e-7,#N=3000,V=40e-16
#     9.35e-6,#N=3000,V=40e-17
#     9.35e-5,#N=3000,V=40e-18
#     9.35e-4,
#     9.35e-4*3

#     ])

# Cx300=numpy.array([
#     9.35e-10,#N=300,V=40e-14
#     9.35e-9,#N=300,V=40e-15
#     9.35e-8,#16
#     9.35e-7,#17
#     9.35e-6,#18
#     9.35e-5,#19
#     9.35e-4,#20
#     9.35e-4*3
# #    3.74e-3,#1e-21
# #    9.35e-3,#4e-21
#     ])

#data_N3000 *= 11696
#data_N300 *= 11696

show()
#savefig('fig1.eps')




#>>> _gfrd.S_irr(.0001 * 1e-8**2/1e-12, 1e-8, 10 * 1e-8 * 1e-12, 1e-12, 1e-8)
#0.99116163945434221


