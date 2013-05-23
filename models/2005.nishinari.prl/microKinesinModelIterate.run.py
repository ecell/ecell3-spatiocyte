#!/usr/bin/env python

import sys
import math
import os
import pickle
import random

#T = requested simulation duration
#V = volume in liters
#N = number of molecules
#R = voxel radius
#D = diffusion coefficient

def run(ActivePercent, Duration, Iteration):
  param = ("--parameters=\"{'ActivePercent':%d, 'Duration':%d, 'Iteration':%d}\"" %(ActivePercent, Duration, Iteration))
  os.system("ecell3-session " + param + " microKinesinModelIterate.py")

if __name__ == '__main__':
  ActivePercent = 0
  Duration = 101
  Iterations = 1000
  for i in range(Iterations):
    run(ActivePercent, Duration, i)
