import scipy.optimize
import yoda
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

scan_sm_p = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable/Rivet_merged.yoda")["/MyAnalysis/hist_e_pT_2d"]
scan_bsm_p = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusn1stable/Rivet_merged.yoda")["/MyAnalysis/hist_e_pT_2d"]

scan_sm_p_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusvestable/Rivet_merged.yoda")["/RAW/MyAnalysis/hist_e_pT_2d"]
scan_bsm_p_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeplusn1stable/Rivet_merged.yoda")["/RAW/MyAnalysis/hist_e_pT_2d"]

scan_sm_m = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable/Rivet_merged.yoda")["/MyAnalysis/hist_e_pT_2d"]
scan_bsm_m = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1stable/Rivet_merged.yoda")["/MyAnalysis/hist_e_pT_2d"]

scan_sm_m_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable/Rivet_merged.yoda")["/RAW/MyAnalysis/hist_e_pT_2d"]
scan_bsm_m_raw = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1stable/Rivet_merged.yoda")["/RAW/MyAnalysis/hist_e_pT_2d"]

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

x = np.arange(0,0.175,0.001)                       #BSM signal fractions to test
w = (np.array(scan_sm_p.xEdges()) - 79.824 )*1000   #candidate W masses
y_p = np.array([])                                  #empty output array for most likely w masses
y_m = np.array([])

hists_sm_p = numpy2d(scan_sm_p)
hists_bsm_p = numpy2d(scan_bsm_p)

hists_sm_m = numpy2d(scan_sm_m)
hists_bsm_m = numpy2d(scan_bsm_m)

# # normalise nominal bin (bin 74) for unweighted histos:
hists_sm_p = ( hists_sm_p / integrate(hists_sm_p[74,:]) ) * ( getNumEntries2d(scan_sm_p_raw) / 40000000 ) * 3464 * 4700 
hists_bsm_p = ( hists_bsm_p / integrate(hists_bsm_p[74,:]) ) * ( getNumEntries2d(scan_bsm_p_raw) / 8000000 ) * 610.1 * 4700

hists_sm_m = ( hists_sm_m / integrate(hists_sm_m[74,:]) ) * ( getNumEntries2d(scan_sm_m_raw) / 40000000 ) * 1909 * 4700
hists_bsm_m = ( hists_bsm_m / integrate(hists_bsm_m[74,:]) ) * ( getNumEntries2d(scan_bsm_m_raw) / 8000000 ) * 615.6 * 4700

#loop through different (sm) templates (labelled by i) with one (sm + bsm) data
for r in x:

  hists = hists_sm_p + r * hists_bsm_p

  data = hists[74,:] #data histogram, 74 chosen as corresponds to mw = 79.824 GeV
  
  lls = np.array([]) #log likelihoods

  # # loop over potential w mass bins
  for i in range(0,hists.shape[0]):

    temp = hists_sm_p[i,:] *np.sum(data)/np.sum(hists_sm_p[i,:]) #template histogram

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
  y_p = np.append(y_p, yi)

#loop through different (sm) templates (labelled by i) with one (sm + bsm) data
for r in x:

  hists = hists_sm_m + r * hists_bsm_m

  data = hists[74,:] #data histogram, 74 chosen as corresponds to mw = 79.824 GeV
  
  lls = np.array([]) #log likelihoods

  # # loop over potential w mass bins
  for i in range(0,hists.shape[0]):

    temp = hists_sm_m[i,:] *np.sum(data)/np.sum(hists_sm_m[i,:]) #template histogram

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
  y_m = np.append(y_m, yi)

# do plus - minus
y = y_p - y_m

# fit to linear
def linear(x,a,b):
  return a*x + b

[a,b], pcov = scipy.optimize.curve_fit(linear,x,y)
print("a = " + str(a) + " ± " + str(pcov[0,0]))
print("b = " + str(b) + " ± " + str(pcov[1,1]))

plt.plot(x,y,label = "data")
plt.plot(x,linear(x,a,b),label = "fit", ls="--")
plt.title("Likelihood fit to SM template using $e$ ${p}_T$ (Atlas)")
plt.xlabel("BSM fraction in signal")
plt.ylabel("$m_{W^+}-m_{W^-}$ (MeV)")
plt.legend()
plt.savefig("likelihood_fits/atlas_mw_plus_minus_mw_minus.pdf")




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

# # luminosities:
# CDF: 8800 pb-1
# Atlas: 4700 pb-1
# Atlas at 13.6TeV: 200000 pb-1