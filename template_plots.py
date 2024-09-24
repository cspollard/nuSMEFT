import yoda
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

#scan = yoda.read(argv[1])["/MyAnalysis/hist_e_pT_2d"]
scan_sm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusvestable_cdf/Rivet_merged.yoda")["/MyAnalysis/hist_e_pT"]
scan_bsm = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusn1stable_cdf/Rivet_merged.yoda")["/MyAnalysis/hist_e_pT"]

# neglects overflows
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

hists_sm = numpy1d(scan_sm)
#hists_bsm = numpy2d(scan_bsm) * integrate(numpy2d(scan_sm)[n,:])/integrate(numpy2d(scan_bsm)[n,:])
hists_bsm = numpy1d(scan_bsm)

print(integrate(hists_sm))
print(integrate(hists_bsm))

#x = np.linspace(0, hists_sm[0,-1], hists_sm[0,:].size)

nums = [0,0.2,0.4,0.6,0.8,100]

for i in nums:
  hist = hists_sm + i * hists_bsm
  hist = hist * 1/integrate(hist)
  plt.stairs(hist)

plt.savefig("template_plot.png")