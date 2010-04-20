import numpy
import csv
import math
from matplotlib import rc
from pylab import *

labelFontSize = 14
tickFontSize = 14
legendFontSize = 14
lineFontSize = 14

fileNames = ['Fig.volume.rebind.csv.back']
legendTitles = ['0% fraction', 't^{-3/2}', '34% fraction']
speciesList = ['E']
lines = ['-', '--', '-', '-']
colors = ['#aa0077', 'b', 'g', 'r', 'c', 'k', '#009955', '#ff9933', '#ff00ff', '#11dd00']

data = csv.reader(open(fileNames[0]))
timePoints = []
input_data = []
for row in data:
  input_data.append(float(row[1]))
data = numpy.log10([value for value in input_data if value != 0])
n, bins = numpy.histogram(data, bins=20)
x = (10 ** bins[1: ] + 10 ** bins[: -1]) * 0.5
dx = (10 ** bins[1: ] - 10 ** bins[: -1])
y = n / float(len(input_data)) / dx
subplot(111)
plot(x, y, ls=lines[0], color=colors[0], label=legendTitles[0], linewidth=2)

for row in x:
  y = x**(-3.0/2.0)

plot(x, y, ls=lines[1], color=colors[1], label=legendTitles[1], linewidth=2)

fileNames = ['Fig.volume.rebind.34fraction.csv.back']

data = csv.reader(open(fileNames[0]))
timePoints = []
input_data = []
for row in data:
  input_data.append(float(row[1]))
data = numpy.log10([value for value in input_data if value != 0])
n, bins = numpy.histogram(data, bins=20)
x = (10 ** bins[1: ] + 10 ** bins[: -1]) * 0.5
dx = (10 ** bins[1: ] - 10 ** bins[: -1])
y = n / float(len(input_data)) / dx
subplot(111)
plot(x, y, ls=lines[2], color=colors[2], label=legendTitles[2], linewidth=2)


ax = gca()
ax.set_xscale('log')
ax.set_yscale('log')
leg = legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')
xlabel('t (s)')
ylabel('rebind frequency')

show()


