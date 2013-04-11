import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from matplotlib.ticker import MaxNLocator

labelFontSize = 14
tickFontSize = 14
legendFontSize = 14
lineFontSize = 14

fileNames = ['actual100nM.csv', "estimate100nM.csv"]
legendTitles = ['MTKinesin (no diffusion)', 'MTKinesin (diffusion)', 'MTKinesinADP (diffusion)', 'Kinesin (diffusion)']
lines = ['-', '-', '-', '-']
colors = ['r', 'b', 'g', 'c', 'r', 'k', '#009955', '#ff9933', '#ff00ff', '#11dd00']

cnt = 0
for i in range(len(fileNames)):
  data = genfromtxt(fileNames[i], delimiter=',').T 
  colSize = len(data)-1
  for j in range(colSize):
    print cnt
    plot(data[0], data[j+1], ls=lines[j], color=colors[cnt], label=legendTitles[cnt], linewidth=1)
    cnt = cnt + 1

ax = gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
leg = legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')
xlabel('site')
ylabel('density')
show()

