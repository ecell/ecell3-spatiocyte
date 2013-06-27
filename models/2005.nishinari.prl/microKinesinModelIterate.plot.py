import numpy as np
import csv
import glob
from pylab import *

labelFontSize = 18
tickFontSize = 18
legendFontSize = 18
lineFontSize = 18

percent = [0, 20, 40, 60, 70, 80, 90, 100]
f = open('tmp.csv', 'wb') 
writer = csv.writer(f, delimiter='\t')
for i in percent:
  names = glob.glob(("t101csv/SingleNeurite.HisLog_p%d" %i)+"_t101_i*.csv")
  nFiles = len(names)
  f = open(names[0], 'r')
  legendTitles = f.readline().strip().split(",")
  binInterval = float(legendTitles[1])
  cols = len(legendTitles)
  data = np.genfromtxt(names[0], delimiter=',', skiprows=1)
  rows = len(data)
  x = []
  logStart = float(data[0][0])
  for row in data:
    if(row[0] == logStart): 
      x.append((float(row[1])+0.5))
    else:
      break
  bins = len(x)
  stride = 1
  ratioTimePoints = rows/bins/stride
  minus = [0.0]*ratioTimePoints
  plus = [0.0]*ratioTimePoints
  for cnt in range(nFiles):
    data = np.genfromtxt(names[cnt], delimiter=',', skiprows=1)
    for j in range(ratioTimePoints):
      minusIndex = j*stride*bins
      plusIndex = j*stride*bins+(bins-1)
      #print "minusIndex:", minusIndex
      #print "plusIndex:", plusIndex
      for k in range(5):
        minus[j] = minus[j] + float(data[minusIndex][k+2])
        plus[j] = plus[j] + float(data[plusIndex][k+2])
  ratios = [0.0]*(ratioTimePoints+1)
  ratios[0] = i
  y = [0.0]*ratioTimePoints
  x = [0.0]*ratioTimePoints
  for j in range(ratioTimePoints):
    ratios[j+1] = plus[j]/max(minus[j], 1.0)
    x[j] = j*stride
    y[j] = ratios[j+1]
  writer.writerow(ratios)
  plot(x, y, lw=2, label=str(i)+"% active MTs")

ax = gca()
ax.grid(color='k', linestyle='--')
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

xlabel('Time (s)', fontsize=lineFontSize)
ylabel('Kinesin Number Plus/Minus End Ratio', fontsize=lineFontSize)
title('Kinesin Plus End Localization by Active MTs')
show()
