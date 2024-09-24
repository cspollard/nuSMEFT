#!/usr/bin/env python3
import scipy.optimize
import yoda
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import scipy

#import histos
scan_sm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusve_lhcrun3/Rivet_merged_30_08.yoda")["/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
scan_bsm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1_lhcrun3/Rivet_merged_30_08.yoda")["/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
# scan_sm_weighted = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusve_lhcr3un/Rivet_merged.yoda")["/MyAnalysis/e_eta_weight2_before"]

scan_sm_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusve_lhcrun3/Rivet_merged_30_08.yoda")["/RAW/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
scan_bsm_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1_lhcrun3/Rivet_merged_30_08.yoda")["/RAW/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]

def numEntries1d(h):
  tot = 0
  for xidx in range(1,h.numBinsX()+1):
    tot+=h.bin(xidx).numEntries()
  return tot

def numpy1d(h):
  return \
    np.array \
    (  [ h.bin(xidx).val()
          for xidx in range(1, h.numBinsX()+1)
        ]
      
    )

def numEntries2d(h):
  tot = 0
  for xidx in range(1,h.numBinsX()+1):
    for yidx in range(1,h.numBinsY()+1):
      tot+=h.bin(xidx,yidx).numEntries()
  return tot
    

def numpy2d(h):
  return \
    np.array \
    ( [ [ h.bin(xidx, yidx).val()
          for yidx in range(1, h.numBinsY()+1)
        ]
        for xidx in range(1, h.numBinsX()+1)
      ]
    )

def integrate(arr):
  integral = 0
  for val in arr:
    integral += val
  return integral

#normalise to number of events, i.e. divide sm histos by no scripts combined
hist_sm= numpy1d(scan_sm)
hist_bsm= numpy1d(scan_bsm)
# hist_sm_weighted = numpy1d(scan_sm_weighted)

hist_sm *= (1/integrate(hist_sm)) * (scan_sm_raw.numEntries() / 25000000) * 200000 * 283.2
hist_bsm *= (1/integrate(hist_bsm)) * (scan_bsm_raw.numEntries() / 5000000) * 200000 * 75.6
# hist_sm_weighted *= (1/integrate(hist_sm_weighted)) * (scan_sm_raw.numEntries() / 40000000) * 8800 * 866.4

print(scan_sm_raw.numEntries())
# hist_3 *= integrate(hist_bsm)/integrate(hist_3)

hist_sm *= 1/integrate(hist_sm)
hist_bsm *= 1/integrate(hist_bsm)

edges = scan_sm.xEdges()

# bsm fraction
i = 0.48
hist_sm_bsm = (hist_sm + i * hist_bsm)

plt.stairs( hist_sm ,edges,label="sm")
plt.stairs( hist_bsm,edges,label="bsm")
# plt.stairs(hist_sm_weighted, edges,label = "sm, weighted",color = "xkcd:neon green",ls="--")
# plt.stairs(hist_sm_bsm, edges, label = "sm + 0.48*bsm",color = "red")
plt.xlabel(r"$2L_p -1$")
# plt.ylabel(r"count (8.8 fb$^{-1}$)")
plt.ylabel("count (normalised to unity)")
plt.title("$2L_p -1$ for $W^-$ events with \n $p_T^W$ > 50 GeV (reconstructed, after cuts)")
plt.legend()
plt.savefig("python_plots_run3/Lp_reco_cuts_50_minus.pdf",bbox_inches="tight")

# # # cross-sections
# # CDF:
# w- sm:  866.4 pb
# w+ sm:  3464 pb
# w- bsm: 866.4 pb
# w+ bsm: 202.0 pb
# # Atlas at 7TeV:
# w- sm:  1909 pb
# w+ sm:  3464 pb
# w- bsm: 615.6 pb
# w+ bsm: 610.1 pb
# # Atlas at 13.6TeV:
# w- sm:  283.2 pb
# w+ sm:  397.8 pb
# w- bsm: 75.6 pb
# w+ bsm: 90.3 pb