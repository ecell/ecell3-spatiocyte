import numpy
import csv
import math
from matplotlib import rc
from pylab import *

labelFontSize = 14
tickFontSize = 14
legendFontSize = 14
lineFontSize = 14

fileNames = ['Fig.volume.irreversible.csv.back']
legendTitles = ['Rebind Time', 't^{-3/2}', 'Spatiocyte (ka=3.77e-19, D=2.5e-12)', 'ODE (knet=3.43e-19) (koff=5)']
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

iterations = 88000.0

for row in x:
  y = y/iterations
plot(x, y, ls=lines[1], color=colors[1], label=legendTitles[1], linewidth=2)


ax = gca()
ax.set_xscale('log')
leg = legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')

show()


