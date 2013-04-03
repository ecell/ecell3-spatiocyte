import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from matplotlib.ticker import MaxNLocator
import glob
#rc('text', usetex=True)


def P(r, D, t, e):
  return 2*r/(4*D*t+4*e*e)*exp(-r*r/(4*D*t+4*e*e))

def Prob(r, p, D, dt, e):
  aProb = []
  for i in range(len(r)):
      prob = 0
      for j in range(len(D)):
        prob = prob + p[j]*P(r[i], D[j]*1e+12, dt, e)
      aProb.append(prob)
  return aProb

labelFontSize = 18
tickFontSize = 18
legendFontSize = 18
lineFontSize = 18
data = []
filenames = glob.glob('2013.matsuoka.pcb.Fig4B.sim.anionic.DEO5.500000e+05.DT9.000000e-08.PRO1.700000e+00.csv*')
for name in filenames:
  tmp = genfromtxt(name, delimiter=',').T
  for i in range(len(tmp)-1):
    data.extend(tmp[i+1])
tmp = []
for i in data:
  if(float(i) != 0.0):
    tmp.append(i)
data = tmp
#print data
data = numpy.array(data)*1e+6
print len(data)
n, bins, patches = hist(data, 40, color='r', histtype='step', normed=1, facecolor='green', alpha=0.75)
clf()

n = numpy.concatenate(([0],n))
b = numpy.array(bins[1:])
width = (b[1]-b[0])/2
b = b-width
b = numpy.concatenate(([0],b))
plot(b, n, color='k', ls='--', label='simulation')

# 1 state model
dt = 33e-3
D = [0.14773e-12]
p = [1]
e = 0
r = []
step = b[len(b)-1]/1000.0
for i in range(1000):
  r.append(step*i)
prob = Prob(r, p, D, dt, e) 
plot(r, prob, color='g', label='1 state')

# 2 state model:
D = [0.067e-12]
p = [0.884]
#e = 0.0362
e = 0
prob1 = numpy.array(Prob(r, p, D, dt, e)) 
#plot(r, prob1, color='b')
D = [0.435e-12]
p = [0.116]
prob2 = numpy.array(Prob(r, p, D, dt, e)) 
#plot(r, prob2, color='g')
plot(r, prob1+prob2, color='b', label='2 states')
D = [0.435e-12]
p = [1]
prob3 = numpy.array(Prob(r, p, D, dt, e)) 
#plot(r, prob3, color='g', ls='--')

# 3 state model:
D = [0.035e-12, 0.133e-12, 0.693e-12]
p = [0.479, 0.480, 0.041]
e = 0
prob = Prob(r, p, D, dt, e) 
plot(r, prob, color='r', label='3 states')

xscale('log')
ax = gca()
ax.grid(color='#999999', linestyle='--')
leg = legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')

#ax.yaxis.set_major_locator(MaxNLocator(14))
xticks(fontsize=tickFontSize)
yticks(fontsize=tickFontSize)
xlabel(r'Displacement, $\Delta r$ ($\mu$m)', fontsize=lineFontSize)
ylabel(r'Probability density, $P$', fontsize=lineFontSize)
suptitle('Non-polarized cells', fontsize=lineFontSize)

show()


