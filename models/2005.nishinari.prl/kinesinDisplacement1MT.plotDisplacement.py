import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from matplotlib.ticker import MaxNLocator
import glob
#rc('text', usetex=True)


labelFontSize = 18
tickFontSize = 18
legendFontSize = 18
lineFontSize = 18
data = []
for name in glob.glob('LifetimeLog.csv*'):
  data.extend((genfromtxt(name, delimiter=',').T)[1])
#print data
data = numpy.array(data)*1e+6
data2 = []
for i in range(len(data)):
  if data[i] > 0.1: 
    data2.append(data[i])
data = numpy.array(data2)
n, bins, patches = hist(data, 45, color='r', normed=1, facecolor='green', alpha=0.75)
#n, bins, patches = hist(data, 50, color='r', histtype='step', normed=1, facecolor='green', alpha=0.75)

n = numpy.concatenate(([0],n))
b = numpy.array(bins[1:])
width = (b[1]-b[0])/2
b = b-width
b = numpy.concatenate(([0],b))
#plot(b, n, color='k')

#xscale('log')
#xlim(0,7)
ax = gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
xlabel(r'Displacement (mean:%.2f), $\Delta r$ ($\mu$m)' %(numpy.mean(data)), fontsize=lineFontSize)
ylabel(r'Probability density, $P$', fontsize=lineFontSize)
show()

