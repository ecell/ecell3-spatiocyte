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
for name in glob.glob('2013.matsuoka.pcb.Fig4B.sim.csv1??'):
  data.extend((genfromtxt(name, delimiter=',').T)[1])
#print data
data = numpy.array(data)*1e+6
n, bins, patches = hist(data, 54, histtype='step', facecolor='green', alpha=0.75)
b = numpy.array(bins[1:])
width = (b[1]-b[0])/2
plot(b-width, n)

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
  for i in range(colSize):
      plot(data[0][0], data[n][i+1]*100, ls=lines[0], color=colors[nSpecies*n+i], label=legendTitles[n], linewidth=1, marker='o')

grid(True)
xscale('log')
show()

