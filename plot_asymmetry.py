#!/usr/bin/env python3

import yoda
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

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

#read in data from yoda
scan_bsm_plus = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusn1stable_cdf/Rivet_merged.yoda")["/MyAnalysis/e_eta_before"]
scan_bsm_minus = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1stable_cdf/Rivet_merged.yoda")["/MyAnalysis/e_eta_before"]
scan_sm_plus = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable_cdf/Rivet_merged.yoda")["/MyAnalysis/e_eta_before"]
scan_sm_minus = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/MyAnalysis/e_eta_before"]
scan_sm_plus_weighted = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable_cdf/Rivet_merged.yoda")["/MyAnalysis/e_eta_weight2_before"]
scan_sm_minus_weighted = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/MyAnalysis/e_eta_weight2_before"]
edges = scan_sm_plus.xEdges()

scan_bsm_plus_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusn1stable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/e_eta_before"]
scan_bsm_minus_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1stable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/e_eta_before"]
scan_sm_plus_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/e_eta_before"]
scan_sm_minus_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/e_eta_before"]
scan_sm_plus_weighted_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/e_eta_weight2_before"]
scan_sm_minus_weighted_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/e_eta_weight2_before"]

#convert to numpy arrays and div by number of merged yodas
hist_sm_plus = numpy1d(scan_sm_plus)
hist_sm_minus = numpy1d(scan_sm_minus)
hist_sm_plus_weighted = numpy1d(scan_sm_plus_weighted)
hist_sm_minus_weighted = numpy1d(scan_sm_minus_weighted)
hist_bsm_plus = numpy1d(scan_bsm_plus)
hist_bsm_minus = numpy1d(scan_bsm_minus)

hist_sm_plus *= (1 / integrate(hist_sm_plus)) * (scan_sm_plus_raw.numEntries() / 40000000) * 866.4*8800
hist_sm_minus *= (1 / integrate(hist_sm_minus)) * (scan_sm_minus_raw.numEntries() / 40000000) * 866.4*8800
hist_bsm_plus *= (1 / integrate(hist_bsm_plus)) * (scan_bsm_plus_raw.numEntries() / 8000000) * 202*8800
hist_bsm_minus *= (1 / integrate(hist_bsm_minus)) * (scan_bsm_minus_raw.numEntries() / 8000000) * 202*8800
hist_sm_plus_weighted *= (1 / integrate(hist_sm_plus_weighted)) * (scan_sm_plus_weighted_raw.numEntries() / 40000000) * 866.4*8800
hist_sm_minus_weighted *= (1 / integrate(hist_sm_minus_weighted)) * (scan_sm_minus_weighted_raw.numEntries() / 40000000) * 866.4*8800

#do plus - minus
hist_bsm_diff = ( hist_bsm_plus - hist_bsm_minus ) / ( hist_bsm_plus + hist_bsm_minus )
hist_sm_diff = ( hist_sm_plus - hist_sm_minus ) / ( hist_sm_minus + hist_sm_plus )
hist_sm_diff_weighted = ( hist_sm_plus_weighted - hist_sm_minus_weighted ) / ( hist_sm_plus_weighted + hist_sm_minus_weighted )

#create histo of sm + f*bsm
f=0.48
hist_comb_plus = ( hist_sm_plus + f*hist_bsm_plus )
hist_comb_minus = ( hist_sm_minus + f*hist_bsm_minus )
hist_comb_diff = ( hist_comb_plus - hist_comb_minus ) / ( hist_comb_plus + hist_comb_minus ) 


plt.stairs(hist_sm_diff,edges, label='sm')
plt.stairs(hist_bsm_diff,edges, label='bsm')
# plt.stairs(hist_sm_diff_weighted,edges, label='sm, weighted')
# plt.stairs(hist_comb_diff,edges, label='sm + 0.48 * bsm',ls="--")
plt.xlabel("$\eta$ (Pseudorapidity)")
plt.ylabel("Lepton Asymmetry")
plt.title("Lepton Asymmetry CDF")
plt.legend()
plt.savefig("python_plots_cdf/asymmetry_plot.pdf")