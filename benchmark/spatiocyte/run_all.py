#!/usr/bin/env python
import sys
import numpy
import math
import run_single

prefix = 'out_'

T = 10.
REPEAT = 1
R = 2.5e-9
D = 1e-12

def run_set(outfile, name, V_list, N_list, T_list):
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        outfile.write('# T=%g, N=%g, V=%g\n' % 
                      (T_list[i], N_list[i], V_list[i]))
        run_times = []
        est_times = []
        for c in range(REPEAT):
            run_time = run_single.run_single(T_list[i], V_list[i], N_list[i], R, D)
            est_time = run_time * (T / T_list[i])
            run_times.append(run_time)
            est_times.append(est_time)
        outfile.write('# run_times = %s\n' % str(run_times))
        outfile.write('%s,\n' % str(est_times))
        outfile.flush()
    outfile.write(']\n')

Vv = [2e-14, ] * 11
Nv = [100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000]
Tv = [max(1e-5, min(T, 1e0 / math.pow(N, 2.0 / 3.0))) for N in Nv]
Vc = [3.33e-15,1e-14, 3.33e-14,1e-13, 3.33e-13,1e-12, 3.33e-12,1e-11, 3.33e-11,1e-10, 3.33e-10,1e-9]#,3.33e-9]
Nc = [100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000,30000000]#,100000000]
Tc = [max(1e-3, min(T, 1e1 / math.pow(float(N), 1.0))) for N in Nc]
V300 = [40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19, 40e-20, 13.3e-20]#,40e-22]
N300 = [300, ] * 8
#T300 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]#, 1e-8]
T300 = [1e7 / math.pow(1.0/V, 2.0 / 3.0) for V in V300]
V3000 = [40e-13, 40e-14, 40e-15, 40e-16, 40e-17, 40e-18, 40e-19, 13.3e-19]#,40e-21]
N3000 = [3000, ] * 8
#T3000 = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-4, 1e-5, 1e-6]#, 1e-7]
T3000 = [2e6 / math.pow(1.0/V, 2.0 / 3.0) for V in V3000]

mode = "Vspatiocyte"
outfile = open(prefix+mode+'.py','w'); 
dataname = 'data_' + mode
run_set(outfile, dataname, Vv, Nv, Tv); outfile.write('\n\n')

