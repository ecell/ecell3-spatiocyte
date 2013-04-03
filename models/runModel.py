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

def run(DT, DEO, PRO, FileStartCount, dt):
  param = ("--parameters=\"{'DT':%e, 'DEO':%e, 'PRO':%e, 'FileStartCount':%d, 'LogInterval':%e}\"" %(DT, DEO, PRO, FileStartCount, dt))
  os.system("ecell3-session " + param + " model.py &")

if __name__ == '__main__':
  dt = 33e-3
  params = [1e-7, 5e+5, 1.8, 4000], [5e-7, 3e+5, 3.7, 4000], [8e-8, 3.5e+5, 1.8, 4000], [9e-8, 5.5e+5, 1.7, 4000]]

  for p in params:
    run(p[0], p[1], p[2], p[3], dt)
