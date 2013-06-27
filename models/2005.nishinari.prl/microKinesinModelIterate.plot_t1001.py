import numpy as np
import csv
import glob
from pylab import *

percent = [30, 40, 50, 60, 100]
f = open('tmp.csv', 'wb') 
writer = csv.writer(f, delimiter='\t')
for i in percent:
  names = glob.glob(("SingleNeurite.HisLog_p%d" %i)+"_t1001_i*.csv")
  names = names + glob.glob(("tcs5/SingleNeurite.HisLog_p%d" %i)+"_t1001_i*.csv")
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
  stride = 10
  ratioTimePoints = rows/bins/stride
  minus = [0.0]*ratioTimePoints
  plus = [0.0]*ratioTimePoints
  for cnt in range(nFiles):
    data = np.genfromtxt(names[cnt], delimiter=',', skiprows=1)
    print names[cnt]
    for j in range(ratioTimePoints):
      minusIndex = (j+1)*stride*bins
      plusIndex = (j+1)*stride*bins+(bins-1)
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
  plot(x, y)

show()
