import numpy as np
import csv
import glob
from pylab import *

labelFontSize = 18
tickFontSize = 18
legendFontSize = 18
lineFontSize = 18
lines = ['-', '-', '-', '-']
colors = ['r', 'g', 'b', 'k', 'c', 'y']

names = sorted(glob.glob("hislog.Translocating.Neurite?.csv"))
nFiles = len(names)
f = open(names[0], 'r')
legendTitles = f.readline().strip().split(",")
cols = len(legendTitles)
legendTitles = ['Neurite 1', 'Neurite 2', 'Neurite 3', 'Neurite 4']
data = np.genfromtxt(names[0], delimiter=',', skiprows=1)
rows = len(data)
bins = 0
logStart = float(data[0][0])
for row in data:
  if(row[0] == logStart): 
    bins = bins + 1
  else:
    break
startTime = 0
endTime = 10
stride = 3
timePoints = int((endTime-startTime)*60.0/stride)
startTimePoint = startTime*60/stride
binIndex = bins-1
includedBins = 1

for cnt in range(nFiles):
  data = np.genfromtxt(names[cnt], delimiter=',', skiprows=1)
  times = [0.0]*timePoints
  values = [0.0]*timePoints
  for i in range(timePoints):
    for j in range(includedBins):
      index = bins*(i+startTimePoint)*stride-j-1+bins
      times[i] = data[index][0]/3600.0-startTime
      for j in range(cols-2):
        values[i] = values[i]+data[index][j+2]
  plot(times, values, lw=2, ls=lines[0], color=colors[cnt], label=legendTitles[cnt])

ax = gca()
ax.grid(color='y', linestyle='--')
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

xlabel('Time (hours)', fontsize=lineFontSize)
ylabel('Number of molecules', fontsize=lineFontSize)
title('Kinesin Neurite Tip Localization')
show()
