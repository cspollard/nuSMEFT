#!/usr/bin/env python3
import scipy.optimize
import yoda
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import scipy

def numpy1d(h):
  return \
    np.array \
    (  [ h.bin(xidx).val()
          for xidx in range(1, h.numBinsX()+1)
        ]
      
    )

def integrate(arr):
  integral = 0
  for val in arr:
    integral += val
  return integral

#import histos
scan_1 = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusn1stable_cdf/Rivet_merged.yoda")["/MyAnalysis/abs_u_pz_over_abs_d_pz_before"]
scan_2 = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable_cdf/Rivet_merged_weight_1.yoda")["/MyAnalysis/abs_u_pz_over_abs_d_pz_before"]
scan_3 = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable_cdf/Rivet_merged_weight_1.yoda")["/MyAnalysis/abs_u_pz_over_abs_d_pz_weighted_before"]

#normalise to number of events, i.e. divide sm histos by 5
hist_1=numpy1d(scan_1)
hist_2=numpy1d(scan_2) / 5
hist_3=numpy1d(scan_3) / 5
#preserve number of events in weighted graph vs unweighted
hist_3 = hist_3 * integrate(hist_2)/integrate(hist_3)
edges = scan_1.xEdges()
bins = np.linspace(-2.475,2.475,100)

# bsm fraction
i = 1000000000000
hist_sm_bsm = (hist_2 + i * hist_1) / (1 + i)

weights = hist_sm_bsm / hist_2

x = np.linspace(-2.5,2.5,100)

# def f(x,b):
#   return 0.5 * (x - b)**2 * (x + 0.5)**2 + 1.25

# a, pcov = scipy.optimize.curve_fit(f,bins,weights)

# for i in range(0,len(a)):
#   print("param " + str(i) + ": " + str(a[i]) + ", std dev " + str( np.sqrt( pcov[i,i] ) ))

def f(x):
  return  0.8 * (x - 0.75)**2 * (x + 1.2)**2 + 1.25

plt.stairs(weights, edges)
plt.plot(x,f(x))
plt.savefig("python_plots/weights_test.pdf")

# with open('weights/weightings_wplus_quark_pz.txt', 'w') as f_plus:
#   for x in weights:
#     f_plus.write(str(x) + "\n")

# def fct(x,frac):
#   a = 3.52
#   b = 5.92
#   c = 3.65
#   d = 0.585
#   e = 2.15
#   g = 1.25

#   weight = np.exp(-1 * (a * x + b)) - d*(1+np.tanh(e*x - c)) + g
#   return (frac * weight + 1)/(1+frac)

# x = np.linspace(-2.5,2.5,100)
# plt.stairs(weights,edges)
# plt.plot(x,fct(x,i))
# plt.savefig("python_plots/weights_test2.pdf")