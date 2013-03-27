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
for name in glob.glob('2013.matsuoka.pcb.Fig4B.sim.csv*'):
  data.extend((genfromtxt(name, delimiter=',').T)[1])
#print data
data = numpy.array(data)*1e+6
n, bins, patches = hist(data, 42, color='r', histtype='step', normed=1, facecolor='green', alpha=0.75)

n = numpy.concatenate(([0],n))
b = numpy.array(bins[1:])
b = b-(b[1]-b[0])/2
b = numpy.concatenate(([0],b))
plot(b, n, color='b')

msd = 0.01
p = []
x = []
for i in range(len(bins)):
  x.append(bins[i])
  p.append(2.0*x[i]/msd*exp(-x[i]*x[i]/msd))

plot(x, p, color='c')

fileNames = ['2013.matsuoka.pcb.Fig4B.exp.csv']
legendTitles = ['PTEN (Matsuoka et al. 2013)', r'PTEN (sim D=0.072 um$^2$s$^{-1}$)']
speciesList = ['3State', '2State', '1State']
lines = ['-', '--', '-', '-']
colors = ['r', 'b', 'y', '#999999', 'm', 'c', 'y', 'k', 'k', 'k', 'k', 'k']

nFiles = len(fileNames)
nSpecies = len(speciesList)
data = []
for n in range(nFiles):
  data.append(genfromtxt(fileNames[n], delimiter=',').T)
  colSize = len(data[n])-1
#  for i in range(colSize):
#      plot(data[0][0], data[n][i+1], ls=lines[0], color=colors[nSpecies*n+i], label=legendTitles[n], linewidth=1, marker='o')

#grid(True)
#xscale('log')
#show()

ax = gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
xlabel(r'Displacement, $\Delta r$ ($\mu$m)', fontsize=lineFontSize)
ylabel(r'Probability density, $P$', fontsize=lineFontSize)
show()

