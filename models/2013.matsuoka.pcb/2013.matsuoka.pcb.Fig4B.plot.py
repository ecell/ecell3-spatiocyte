import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from matplotlib.ticker import MaxNLocator
import glob
#rc('text', usetex=True)


def P(r, D, t, e):
  return 2*r/(4*D*t+4*e*e)*exp(-r*r/(4*D*t+4*e*e))

def Prob(r, p, D, dt, e):
  aProb = []
  for i in range(len(r)):
      prob = 0
      for j in range(len(D)):
        prob = prob + p[j]*P(r[i], D[j]*1e+12, dt, e)
      aProb.append(prob)
  return aProb

labelFontSize = 18
tickFontSize = 18
legendFontSize = 18
lineFontSize = 18
data = []
for name in glob.glob('2013.matsuoka.pcb.Fig4B.sim.csv*'):
  data.extend((genfromtxt(name, delimiter=',').T)[1])
#print data
data = numpy.array(data)*1e+6
n, bins, patches = hist(data, 60, color='r', histtype='step', normed=1, facecolor='green', alpha=0.75)

n = numpy.concatenate(([0],n))
b = numpy.array(bins[1:])
width = (b[1]-b[0])/2
b = b-width
b = numpy.concatenate(([0],b))
plot(b, n, color='k')

# 1 state model
dt = 33e-3
D = [0.14773e-12]
p = [1]
e = 0
r = []
for i in range(1000):
  r.append(0.004+i*0.0008)
prob = Prob(r, p, D, dt, e) 
plot(r, prob, color='y')

# 2 state model:
D = [0.067e-12]
p = [0.884]
e = 0.0362
prob1 = numpy.array(Prob(r, p, D, dt, e)) 
plot(r, prob1, color='b')
D = [0.435e-12]
p = [0.116]
e = 0.0362
prob2 = numpy.array(Prob(r, p, D, dt, e)) 
plot(r, prob2, color='g')
plot(r, prob1+prob2, color='c')

# 3 state model:
D = [0.035e-12, 0.133e-12, 0.693e-12]
p = [0.479, 0.480, 0.041]
e = 0.036
prob = Prob(r, p, D, dt, e) 
plot(r, prob, color='m')

fileNames = ['2013.matsuoka.pcb.Fig4B.exp.csv']
legendTitles = ['PTEN (Matsuoka et al. 2013)', r'PTEN (sim D=0.072 um$^2$s$^{-1}$)']
speciesList = ['3State', '2State', '1State']
lines = ['-', '--', '-', '-']
colors = ['r', 'b', 'g', '#999999', 'm', 'c', 'y', 'k', 'k', 'k', 'k', 'k']

nFiles = len(fileNames)
nSpecies = len(speciesList)
data = []
for n in range(nFiles):
  data.append(genfromtxt(fileNames[n], delimiter=',').T)
  colSize = len(data[n])-1
  for i in range(colSize):
      plot(data[0][0], data[n][i+1], ls='.', color=colors[nSpecies*n+i], label=legendTitles[n], linewidth=1, marker='o')

#xscale('log')
ax = gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
xlabel(r'Displacement, $\Delta r$ ($\mu$m)', fontsize=lineFontSize)
ylabel(r'Probability density, $P$', fontsize=lineFontSize)
show()

