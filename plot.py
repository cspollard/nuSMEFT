#!/usr/bin/env python3
from matplotlib.figure import Figure
import yoda
import numpy

FIGSIZE = (4, 4)

smfname = "yoda/Rivet_Run3_ve.yoda"
bsmfname = "yoda/Rivet_Run3_n1.yoda"

#import histos
scan_sm = yoda.read(smfname)["/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
scan_bsm = yoda.read(bsmfname)["/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
# scan_sm_weighted = yoda.read("/data/atlas/users/wadh6495/rivet/ppeminusve_lhcr3un/Rivet_merged.yoda")["/MyAnalysis/e_eta_weight2_before"]

scan_sm_raw = yoda.read(smfname)["/RAW/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]
scan_bsm_raw = yoda.read(bsmfname)["/RAW/MyAnalysis/Lp_histogram_50_reconstructed_cuts"]


def numEntries1d(h):
  tot = 0
  for xidx in range(1,h.numBinsX()+1):
    tot+=h.bin(xidx).numEntries()
  return tot

def numpy1d(h):
  values = \
    numpy.array \
    ( [ b.val() for b in h.bins() ]
    )
  return numpy.array(h.xEdges()) , values / h.annotation("ScaledBy")

def integrate(arr):
  integral = 0
  for val in arr:
    integral += val
  return integral

def double(xs, ys):
  return numpy.vstack([xs, ys]).reshape((-1,), order='F')

def h2curve(h):
  xs = double(h[0], h[0])
  ys = numpy.array([0] + list(double(h[1], h[1])) + [0])

  return xs, ys

#normalise to number of events, i.e. divide sm histos by no scripts combined
xs , sm = h2curve(numpy1d(scan_sm))
_ , bsm = h2curve(numpy1d(scan_bsm))
# hist_sm_weighted = numpy1d(scan_sm_weighted)

# Chris: I have no idea what these numbers are.
# hist_sm *= (1/integrate(hist_sm)) * (scan_sm_raw.numEntries() / 25000000) * 200000 * 283.2
# hist_bsm *= (1/integrate(hist_bsm)) * (scan_bsm_raw.numEntries() / 5000000) * 200000 * 75.6
# hist_sm_weighted *= (1/integrate(hist_sm_weighted)) * (scan_sm_raw.numEntries() / 40000000) * 8800 * 866.4

# hist_3 *= integrate(hist_bsm)/integrate(hist_3)

# hist_sm *= 1/integrate(hist_sm)
# hist_bsm *= 1/integrate(hist_bsm)

edges = scan_sm.xEdges()

fig = Figure(FIGSIZE)

plt = fig.add_subplot(3, 1, (1, 2))

plt.plot( xs, sm/integrate(sm), label="$W^- \\to \\ell^- \\nu_L$" , color="blue", lw=2)
plt.plot( xs, bsm/integrate(bsm), label="$W^- \\to \\ell^- \\nu_R$", color="red", lw=2, ls=":")

plt.set_xlim((-1, 1))
plt.set_xticks([])

plt.set_ylabel("$\\frac{1}{\\sigma} \\frac{\\mathrm{d}\\sigma}{\\mathrm{d}(2L_p - 1)}$")

plt.legend(title="$pp$, $\\sqrt{s}=13.6$ TeV")

_, ymax = plt.get_ylim()
plt.set_ylim(0, ymax)

ratioplt = fig.add_subplot(3, 1, 3)

ratioplt.plot((-2, 2), (1, 1), lw=1, color=("gray", 0.5))
ratioplt.plot(xs, sm*integrate(bsm)/bsm/integrate(sm), color="red", lw=2, ls=":")
ratioplt.set_xlim((-1, 1))
ratioplt.set_xlabel(r"$2L_p -1$")
ratioplt.set_ylabel("Ratio")

fig.tight_layout()

fig.savefig("Lp_reco_cuts.pdf",bbox_inches="tight")
