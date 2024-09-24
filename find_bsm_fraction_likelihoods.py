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

def numpy2d(h):
  return \
    np.array \
    ( [ [ h.bin(xidx, yidx).val()
          for yidx in range(1, h.numBinsY()+1)
        ]
        for xidx in range(1, h.numBinsX()+1)
      ]
    )

def rawnumpy1d(h):
  return \
    np.array \
    (  [ h.bin(xidx).numEntries()
          for xidx in range(1, h.numBinsX()+1)
        ]
      
    )

def integrate(arr):
  integral = 0
  for val in arr:
    integral += val
  return integral

eps = 10e-6

# import histos
scan_sm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusve_lhcrun3/Rivet_merged_30_08.yoda")["/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
scan_bsm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1_lhcrun3/Rivet_merged_30_08.yoda")["/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
# scan_sm_weighted = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusve_lhcrun3/Rivet_merged.yoda")["/MyAnalysis/e_abscostheta_0_15_2_reco_weight2"]

scan_sm_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusve_lhcrun3/Rivet_merged_30_08.yoda")["/RAW/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
scan_bsm_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1_lhcrun3/Rivet_merged_30_08.yoda")["/RAW/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]

# read in to numpy arrays
hist_sm = numpy1d(scan_sm)
hist_bsm = numpy1d(scan_bsm)
# hist_sm_weighted = numpy1d(scan_sm_weighted)

# manually normalise histos to cross section
hist_sm *= (1 / integrate(hist_sm)) * (scan_sm_raw.numEntries() / 25000000) * 283.2 * 200000
hist_bsm *= (1 / integrate(hist_bsm)) * (scan_bsm_raw.numEntries() / 5000000) * 75.6 * 200000
# normalise weighted sm histogram to unweighted sm
# hist_sm_weighted *= integrate(hist_sm)/integrate(hist_sm_weighted)


fracs = np.arange(0,0.05,0.000001)

lls = np.array([],dtype=np.double)

for f in fracs:
  
  # read in data as sm + bsm
  data = hist_sm + f * hist_bsm

  # normalise template to data
  # temp = hist_sm_weighted * integrate(data)/integrate(hist_sm_weighted)
  temp = hist_sm * integrate(data)/integrate(hist_sm)

  ll = 0
  for j in range(0,len(hist_sm)):
      
    #signal
    n = data[j]
    #template
    m = temp[j] + eps
  
    ll += n * np.log(m) - m - (n + 1/2) * np.log(n + 1) + n

  lls = np.append(lls,ll)

# get delta log likelihood
print("max log-likelihood: " + str(lls.max()))
lls -= lls.max()

fbest = fracs[lls.argmax()]
print("max log-likelihood fraction: " + str(fbest))
try:
  print("1-sigma down: " + str(fracs[np.abs(lls[:lls.argmax()] + 1/2 ).argmin()]))
except:
  print("error: couldn't print 1-sigma lower limit")
  
print("1-sigma up: " + str(fracs[lls.argmax() + np.abs(lls[lls.argmax():] + 1/2 ).argmin()]))

try:
  print("2-sigma down: " + str( fracs[ np.abs( lls[:lls.argmax()] + 1 ).argmin() ] ))
except:
  print("error: couldn't print 2-sigma lower limit")

print("2-sigma up: " + str( fracs[lls.argmax() + np.abs( lls[lls.argmax():] + 1 ).argmin() ] ))

databest = hist_sm + fbest*hist_bsm
# tempbest = hist_sm_weighted * (integrate(databest) / integrate(hist_sm_weighted))
tempbest = hist_sm * (integrate(databest) / integrate(hist_sm))
plt.stairs(tempbest,scan_sm.xEdges(),label='sm, weighted')
plt.stairs(databest,scan_sm.xEdges(),label='sm + bsm',ls='--')
plt.legend()
plt.savefig("python_plots/test2.pdf")

plt.clf()

plt.plot(fracs,lls)
plt.savefig("python_plots/test.pdf")

# # # cross-sections
# # CDF:
# w- sm:  866.4 pb
# w+ sm:  3464 pb
# w- bsm: 866.4 pb
# w+ bsm: 202.0 pb
# # Atlas:
# w- sm:  1909 pb
# w+ sm:  3464 pb
# w- bsm: 615.6 pb
# w+ bsm: 610.1 pb
# # Atlas at 13.6TeV:
# w- sm:  283.2 pb
# w+ sm:  397.8 pb
# w- bsm: 75.6 pb
# w+ bsm: 90.3 pb