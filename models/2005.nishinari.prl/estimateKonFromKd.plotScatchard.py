import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from matplotlib.ticker import MaxNLocator
import glob
import scipy.constants
#rc('text', usetex=True)


labelFontSize = 18
tickFontSize = 18
legendFontSize = 18
lineFontSize = 18
size = 50
T_total = 686699.0
Kd = 100e-9*scipy.constants.N_A*1e+3*1.79e-17
Ka = 1/Kd
KT = []
y = []
for i in range(size):
  KT.append(i/float(size)*T_total)
  y.append(Ka*T_total-Ka*KT[i])

plot(KT, y, color='m')

x = [58542, 107917, 150125, 186623, 218500]
y = [58442.0/101, 107717.0/200, 149825.0/300, 186223.0/395, 218000.0/512]
plot(x, y, color='b')

data = []
for name in glob.glob('IterateLog107917p0.12.csv'):
  data.extend((genfromtxt(name, delimiter=',').T)[1])
data = numpy.array(data)
print "mean:", numpy.mean(data)
print "#tubulin:", (genfromtxt(name, delimiter=',').T)[2][0]

#xscale('log')
ax = gca()
ax.grid(color='b', linestyle='--')
xlabel(r'KT', fontsize=lineFontSize)
ylabel(r'KT/K', fontsize=lineFontSize)

show()

