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
for name in glob.glob('IterateLog14.5.csv*'):
  data.extend((genfromtxt(name, delimiter=',').T)[0])
#print data
data = numpy.array(data)
n, bins, patches = hist(data, 48, color='r', normed=1, facecolor='green', alpha=0.75)

n = numpy.concatenate(([0],n))
b = numpy.array(bins[1:])
width = (b[1]-b[0])/2
b = b-width
b = numpy.concatenate(([0],b))
#plot(b, n, color='k')

#xscale('log')
ax = gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
xlabel(r'Residence Time (mean:%.2f), s' %(numpy.mean(data)), fontsize=lineFontSize)
ylabel(r'Probability density, $P$', fontsize=lineFontSize)
show()

