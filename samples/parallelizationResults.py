import numpy
import csv
import math
from matplotlib import rc
from pylab import *

labelFontSize = 14
tickFontSize = 14
legendFontSize = 14
lineFontSize = 14

x = [1, 8, 64, 512, 4096]
y = [419.5, 65.85, 9.907, 1.657, 0.3118]

y2 = []
for val in x:
  y2.append(y[0]/val)

legendTitles = ['simulation', 'ideal']
lines = ['-', '--', '-', '-']
colors = ['#aa0077', 'k', 'g', 'r', 'c', 'k', '#009955', '#ff9933', '#ff00ff', '#11dd00']

plot(x, y, ls=lines[0], color=colors[0], label=legendTitles[0], linewidth=2)
plot(x, y2, ls=lines[1], color=colors[1], label=legendTitles[1], linewidth=2)

ax = gca()
leg = legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')
xlabel('# processors')
ylabel('time, t(s)')
ax.set_xscale('log')
ax.set_yscale('log')

show()


