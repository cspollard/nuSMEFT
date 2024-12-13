import yoda
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

scan_sm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/MyAnalysis/hist_pT_miss_2d"]
scan_bsm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1stable_cdf/Rivet_merged.yoda")["/MyAnalysis/hist_pT_miss_2d"]
# scan_sm_weighted = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/MyAnalysis/hist_pT_miss_2d_weight2"]

scan_sm_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/hist_pT_miss_2d"]
scan_bsm_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1stable_cdf/Rivet_merged.yoda")["/RAW/MyAnalysis/hist_pT_miss_2d"]

def getNumEntries2d(h):
  # assumes each y histogram has the same number of entries, which should be true in this rivet analysis
  numEntries = 0
  for yidx in range(1, h.numBinsY()+1):
    numEntries += h.bin(1,yidx).numEntries()
  return numEntries

# neglects overflows
def numpy2d(h):
  return \
    np.array \
    ( [ [ h.bin(xidx, yidx).val()
          for yidx in range(1, h.numBinsY()+1)
        ]
        for xidx in range(1, h.numBinsX()+1)
      ]
    )

eps = 10e-6

def integrate(arr):
  integral = 0.0
  for val in arr:
    integral += val
  return integral

x = np.arange(0,1.365,0.005)                       #BSM signal fractions to test
w = (np.array(scan_sm.xEdges()) - 79.824 )*1000   #candidate W masses
y = np.array([])                                  #empty output array for most likely w masses
y_uerr = np.array([])
y_derr = np.array([])

hists_sm = numpy2d(scan_sm)
hists_bsm = numpy2d(scan_bsm)
# hists_sm_weighted = numpy2d(scan_sm_weighted)

# # normalise nominal bin (bin 74) for unweighted histos:
# TODO: this needs to be 4000000 , 8000000 -> nevents
# 8800: lumi
# 866.4 , 202 : xsec
# 74: nominal histogram
hists_sm = ( hists_sm / integrate(hists_sm[74,:]) ) * ( getNumEntries2d(scan_sm_raw) / 40000000 ) * 8800 * 866.4 
hists_bsm = ( hists_bsm / integrate(hists_bsm[74,:]) ) * ( getNumEntries2d(scan_bsm_raw) / 8000000 ) * 8800 * 202
# hists_sm_weighted = hists_sm_weighted * integrate(hists_sm[74,:]) / integrate(hists_sm_weighted[74,:])

#loop through different (sm) templates (labelled by i) with one (sm + bsm) data
for r in x:

  hists = hists_sm + r * hists_bsm

  data = hists[74,:] #data histogram, 74 chosen as corresponds to mw = 79.824 GeV
  
  lls = np.array([]) #log likelihoods

  # # loop over potential w mass bins
  for i in range(0,hists.shape[0]):

    temp = hists_sm[i,:] #template histogram
    temp = temp*np.sum(data)/np.sum(temp) #normalised template histogram

    # # calculate log likelihood
    ll = 0
    for j in range(0,hists.shape[1]):

      n = data[j]
      m = temp[j] + eps
      ll += n * np.log(m) - m - (n + 0.5) * np.log(n + 1) + n

    lls = np.append(lls, ll)

  # # get delta log likelihood
  lls -= np.max(lls)

  yi = w[np.argmax(lls)]
  # # add most likely mass to mass curve
  y = np.append(y, yi)

  # # find up and down errors, by finding m closest to delta log likelihood = 1/2, above and below argmax
  y_uerr = np.append(y_uerr, w[ np.argmax(lls) + np.abs(lls[np.argmax(lls):] + 1/2).argmin() ])

  try:
    y_derr = np.append(y_derr, w[ np.abs(lls[:np.argmax(lls)] + 1/2).argmin() ])
  except:
    y_derr = np.append(y_derr,y_derr[-1])
    print("tried to get argmin of empty sequence, skipping.")

print(y_uerr)
print(y_derr)

print(x[np.abs(y-76).argmin()])

plt.plot(x,y,label="central value")
plt.fill_between(x,y_derr,y_uerr,color="xkcd:powder blue",label="$1\sigma$")
plt.title("Likelihood fit to SM template using missing ${p}_T$ for $W^-$ (CDF)")
plt.xlabel("BSM fraction in signal")
plt.ylabel("$\Delta m_W$ (MeV)")
plt.legend(loc='upper left')
plt.savefig("likelihood_fits/minus_cdf_ptmiss.pdf")



# does `bsm` mean CNHe = 1?
# # # cross-sections
# # CDF at 2 TeV:
# w- sm:  866.4 pb
# w+ sm:  202.0 pb
# w- bsm: 866.4 pb
# w+ bsm: 202.0 pb
# INCLUSIVE CROSS SECTIONS
# # Atlas 7 TeV:
# w- sm:  1909 pb - validated
# w+ sm:  3464 pb - validated
# w- bsm: 615.6 pb - validated
# w+ bsm: 610.1 pb - validated
# I don't believe these can be correct. Are they flipped with the above?
# these are probably not used anywhere.
# # Atlas at 13.6TeV:
# w- sm:  283.2 pb
# w+ sm:  397.8 pb
# w- bsm: 75.6 pb
# w+ bsm: 90.3 pb

# # luminosities:
# CDF: 8800 pb-1
# Atlas: 4700 pb-1
# Atlas at 13.6TeV: 200000 pb-1
