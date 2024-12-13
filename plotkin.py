import yoda
import numpy
from matplotlib.figure import Figure

lumiCDF = 8800
lumiATLAS = 4700

xsecs = \
  { "CDF" :
    { "SM+" : 866.4
    , "SM-" : 866.4
    , "BSM+" : 202
    , "BSM-" : 202
    }
  , "ATLAS" :
    { "SM+" : 3464
    , "SM-" : 1909
    , "BSM+" : 610.1
    , "BSM-" : 615.6
    }
  }

axisdict = \
  { "pTl" : r"$p_T^\ell$ [GeV]"
  , "mT" : r"$m_T$ [GeV]"
  }

def numpy2d(h):
  values = \
    numpy.array \
    ( [ [ h.bin(xidx, yidx).val()
          for yidx in range(1, h.numBinsY()+1)
        ]
        for xidx in range(1, h.numBinsX()+1)
      ]
    )

  return numpy.array(h.xEdges()) , values / h.annotation("ScaledBy")


def numpy1d(h):
  values = \
    numpy.array \
    ( [ h.bin(xidx).val() for xidx in range(1, h.numBinsX()+1) ]
    )

  return numpy.array(h.xEdges()) , values / h.annotation("ScaledBy")


def double(xs, ys):
  return numpy.vstack([xs, ys]).reshape((-1,), order='F')


def gethists(yodadict):
  sumW = yodadict["/RAW/_EVTCOUNT"].sumW()
  pTl = numpy1d(yodadict["/MyAnalysis/hist_e_pT_reconstructed"])
  mT = numpy1d(yodadict["/MyAnalysis/hist_mT_reconstructed"])
  pTl = (pTl[0] , pTl[1] / sumW)
  mT = (mT[0] , mT[1] / sumW)
  return { "pTl" : pTl , "mT" : mT }


def h2curve(h):
  xs = double(h[0], h[0])
  ys = numpy.array([0] + list(double(h[1], h[1])) + [0])

  return xs, ys


cdfsmplus = gethists(yoda.read("yoda/ppeplusvestable_cdf/Rivet_merged.yoda"))
cdfsmminus = gethists(yoda.read("yoda/ppeminusvestable_cdf/Rivet_merged.yoda"))
cdfbsmplus = gethists(yoda.read("yoda/ppeplusn1stable_cdf/Rivet_merged.yoda"))
cdfbsmminus = gethists(yoda.read("yoda/ppeminusn1stable_cdf/Rivet_merged.yoda"))

atlassmplus = gethists(yoda.read("yoda/ppeplusvestable/Rivet_merged.yoda"))
atlassmminus = gethists(yoda.read("yoda/ppeminusvestable/Rivet_merged.yoda"))
atlasbsmplus = gethists(yoda.read("yoda/ppeplusn1stable/Rivet_merged.yoda"))
atlasbsmminus = gethists(yoda.read("yoda/ppeminusn1stable/Rivet_merged.yoda"))


for k in ["pTl", "mT"]:
  edges , smplus = cdfsmplus[k]
  _ , smminus = cdfsmplus[k]
  _ , bsmplus = cdfbsmplus[k]
  _ , bsmminus = cdfbsmminus[k]

  smpred = \
      smplus * xsecs["CDF"]["SM+"] \
    + smminus * xsecs["CDF"]["SM-"]

  bsmpred = \
      bsmplus * xsecs["CDF"]["BSM+"] \
    + bsmminus * xsecs["CDF"]["BSM-"]

  smpred = smpred * lumiCDF
  bsmpred = bsmpred * lumiCDF

  print("CDF total SM prediction: %.2e events" % smpred.sum())
  print("CDF total BSM prediction: %.2e events" % bsmpred.sum())

  xs , sm = h2curve((edges, smpred))
  _ , bsm = h2curve((edges, bsmpred))

  fig = Figure((6, 6))
  plt = fig.add_subplot(3, 1, (1, 2))

  plt.plot(xs, sm, color="blue", label="SM")
  plt.plot(xs, bsm, color="orange", label="BSM")
  plt.set_xticks([])
  plt.set_ylabel("events / bin")

  ratioplt = fig.add_subplot(3, 1, 3)
  ratioplt.plot(xs, bsm/sm, color="orange")
  ratioplt.set_xlabel(axisdict[k])
  ratioplt.set_ylabel("ratio")

  plt.legend(title="CDF\n$W^+$ and $W^-$\nelectron channel")

  fig.tight_layout()
  fig.savefig("CDF-%s.png" % k)
