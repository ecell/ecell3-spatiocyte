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

fileNames = ['2005.nishinari.prl.Fig3.left.csv']
legendTitles = ['State 2', 'State 1']
speciesList = ['E']
lines = ['-', '--', '-', '-']
colors = ['r', 'b', '#000000', 'black', 'c', 'k', '#009955', '#ff9933', '#ff00ff', '#11dd00']
data = genfromtxt(fileNames[0], delimiter=',', skiprows=1)
x = []
currTime = float(data[0][0])
for row in data:
  if(row[0] == currTime): 
    x.append(float(row[1]))
  else:
    break

bins = len(x)
timePoints = len(data)/bins
y1 = [0]*bins
y2 = [0]*bins
y3 = [0.562]*bins
y4 = [0.38]*bins
for i in range(timePoints):
  for j in range(bins):
    y1[j] = y1[j] + data[i*bins+j][2]
    y2[j] = y2[j] + data[i*bins+j][3]

for i in range(bins):
  y1[i] = y1[i]/timePoints/(700.0/bins)
  y2[i] = y2[i]/timePoints/(700.0/bins)

plot(x, y3, ls=lines[0], color=colors[2], linewidth=1)
plot(x, y4, ls=lines[0], color=colors[2], linewidth=1)
plot(x, y1, ls=lines[0], color=colors[0], label=legendTitles[0], linewidth=1)
plot(x, y2, ls=lines[0], color=colors[1], label=legendTitles[1], linewidth=1)

ax = gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
leg = legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
xlim(0,bins)
ylim(0,0.65)
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')
xlabel('site')
ylabel('density')
show()



