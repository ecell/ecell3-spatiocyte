import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from matplotlib.ticker import MaxNLocator
#rc('text', usetex=True)

labelFontSize = 18
tickFontSize = 18
legendFontSize = 18
lineFontSize = 18

fileNames = ['2013.matsuoka.pcb.Fig1G.exp.csv', '2013.matsuoka.pcb.Fig1G.sim.csv']
legendTitles = ['PTEN (Matsuoka et al. 2013)', r'PTEN (sim D=0.072 um$^2$s$^{-1}$)']
speciesList = ['PTENep', 'PTENsim']
lines = ['-', '--', '-', '-']
colors = ['g', 'b', 'y', '#999999', 'm', 'c', 'y', 'k', 'k', 'k', 'k', 'k']

nFiles = len(fileNames)
nSpecies = len(speciesList)
data = []
for n in range(nFiles):
  data.append(genfromtxt(fileNames[n], delimiter=',').T)
  colSize = len(data[n])-1
  for i in range(colSize):
    if(n):
      plot(data[0][0], data[n][i+1], ls=lines[0], color=colors[nSpecies*n+i], label=legendTitles[n], linewidth=1, marker='o')
    else:
      plot(data[0][0], data[n][i+1]*1e-6, ls=lines[0], color=colors[nSpecies*n+i], label=legendTitles[n], linewidth=1, marker='o')


ax = gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
leg = legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
xlim(2,7)
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')
xlabel('Time (s)', fontsize=lineFontSize)
ylabel('|Displacement| (m)', fontsize=lineFontSize)
show()

