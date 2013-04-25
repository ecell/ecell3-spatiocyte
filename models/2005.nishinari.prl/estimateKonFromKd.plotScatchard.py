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

x = [58542, 107917, 150125, 186623, 218500, 246582, 271511, 293792, 313826, 331939, 363415, 412336, 448633, 545241]
y = [58442.0/100, 107717.0/199, 149825.0/300, 186223.0/395, 218000.0/497, 245982.0/586, 270811.0/695, 292992.0/785, 312926.0/868, 330939.0/955, 363215.0/1160, 410736.0/1633, 446633.0/1889, 541241.0/4027]
plot(x, y, 'bo')

data = []
for name in glob.glob('IterateLog545241p0.128nr.csv'):
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

