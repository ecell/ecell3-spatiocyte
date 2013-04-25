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
T_total = 638703.0
Kd = 100e-9*scipy.constants.N_A*1e+3*1.79e-17
Ka = 1/Kd
KT = []
y = []
for i in range(size):
  KT.append(i/float(size)*T_total)
  y.append(Ka*T_total-Ka*KT[i])

plot(KT, y, color='m')

#xscale('log')
ax = gca()
ax.grid(color='b', linestyle='--')
xlabel(r'KT', fontsize=lineFontSize)
ylabel(r'KT/K', fontsize=lineFontSize)

show()

