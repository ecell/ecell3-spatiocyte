import numpy
import csv
import math
from matplotlib import rc
from pylab import *

labelFontSize = 14
tickFontSize = 14
legendFontSize = 14
lineFontSize = 14

x = [0, 10, 20, 30, 40, 50, 60, 70]
y100 = [100e-12, 83.2e-12, 65.7e-12, 47.7e-12, 18e-12, 0,  0, 0]
y10 = [10e-12, 8.32e-12, 6.58e-12, 4.79e-12, 1.8e-12, 0, 0, 0]

yA = []
for val in y100:
  yA.append(val/y100[0])

yB = []
for val in y10:
  yB.append(val/y10[0])
  

legendTitles = ['D=100e-12', 'D=10e-12']
lines = ['-', '--', '-', '-']
colors = ['#aa0077', 'b', 'g', 'r', 'c', 'k', '#009955', '#ff9933', '#ff00ff', '#11dd00']

plot(x, yA, ls=lines[0], color=colors[0], label=legendTitles[0], linewidth=2)
plot(x, yB, ls=lines[1], color=colors[1], label=legendTitles[1], linewidth=2)

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
xlabel('volume fraction %')
ylabel('normalized D')

show()


