import numpy as np
import csv

percent = [0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
f = open('tmp.csv', 'wb') 
writer = csv.writer(f, delimiter='\t')
for i in percent:
  name = ("SingleNeurite.HisLog%.2f" %i)+"_25s.csv"
  f = open(name, 'r')
  legendTitles = f.readline().strip().split(",")
  f.close()
  binInterval = float(legendTitles[1])
  data = np.genfromtxt(name, delimiter=',', skiprows=1)
  x = []
  logStart = float(data[0][0])
  for row in data:
    if(row[0] == logStart): 
      x.append((float(row[1])+0.5))
    else:
      break
  bins = len(x)
  stride = 25
  ratioTimePoints = len(data)/bins/stride
  ratios = [0]*(ratioTimePoints+2)
  ratios[0] = i*100
  for j in range(ratioTimePoints):
    ratios[j+1] = float(data[j*stride*bins+(bins-1)][4])/float(data[j*stride*bins][4])
  j = ratioTimePoints
  ratios[j+1] = float(data[j*stride*bins-1][4])/float(data[j*stride*bins-bins][4])
  writer.writerow(ratios)
