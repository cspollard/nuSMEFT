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
    ( [ b.val() for b in h.bins() ]
    )

  return numpy.array(h.xEdges()) , values / h.annotation("ScaledBy")


def double(xs, ys):
  return numpy.vstack([xs, ys]).reshape((-1,), order='F')


def gethists(yodadict):
  evtcount = yodadict["/RAW/_EVTCOUNT"]
  sumW = evtcount.sumW()

  hpTl = yodadict["/MyAnalysis/hist_e_pT_reconstructed"]
  hmT = yodadict["/MyAnalysis/hist_mT_reconstructed"]

  pTl = numpy1d(hpTl)
  mT = numpy1d(hmT)

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

  widths = edges[1:] - edges[:-1]

  smpred = \
      smplus * xsecs["CDF"]["SM+"] \
    + smminus * xsecs["CDF"]["SM-"]

  bsmpred = \
      bsmplus * xsecs["CDF"]["BSM+"] \
    + bsmminus * xsecs["CDF"]["BSM-"]

  smpred = smpred * lumiCDF
  bsmpred = bsmpred * lumiCDF

  print("CDF total SM prediction: %.2e events" % (smpred * widths).sum())
  print("CDF total BSM prediction: %.2e events" % (bsmpred * widths).sum())

  xs , sm = h2curve((edges, smpred))
  _ , bsm = h2curve((edges, bsmpred))


  fig = Figure((5, 5))
  plt = fig.add_subplot(3, 1, (1, 2))

  plt.plot(xs, sm, color="blue", label=r"$W \to \ell \nu$", lw=2)
  plt.plot(xs, bsm, color="orange", label=r"$W \to \ell N$ ", lw=2, ls=":")
  plt.set_xticks([])
  plt.set_ylabel("events / GeV")

  ratioplt = fig.add_subplot(3, 1, 3)
  ratioplt.plot(xs, bsm/sm, color="orange", lw=2, ls=":")
  ratioplt.set_xlabel(axisdict[k])
  ratioplt.set_ylabel("ratio")

  plt.legend(title="CDF electron channel")

  fig.tight_layout()
  fig.savefig("CDF-%s.png" % k)
  fig.savefig("CDF-%s.pdf" % k)
  fig.clf()


for k in ["pTl", "mT"]:
  edges , smplus = atlassmplus[k]
  _ , smminus = atlassmplus[k]
  _ , bsmplus = atlasbsmplus[k]
  _ , bsmminus = atlasbsmminus[k]

  widths = edges[1:] - edges[:-1]

  smpredplus = smplus * xsecs["ATLAS"]["SM+"] * lumiATLAS
  smpredminus = smminus * xsecs["ATLAS"]["SM-"] * lumiATLAS

  bsmpredplus = bsmplus * xsecs["ATLAS"]["BSM+"] * lumiATLAS
  bsmpredminus = bsmminus * xsecs["ATLAS"]["BSM-"] * lumiATLAS

  print("ATLAS total SM+ prediction: %.2e events" % (smpredplus * widths).sum())
  print("ATLAS total SM- prediction: %.2e events" % (smpredminus * widths).sum())
  print("ATLAS total BSM+ prediction: %.2e events" % (smpredplus * widths).sum())
  print("ATLAS total BSM- prediction: %.2e events" % (smpredminus * widths).sum())

  xs , smplus = h2curve((edges, smpredplus))
  _ , smminus = h2curve((edges, smpredminus))
  _ , bsmplus = h2curve((edges, bsmpredplus))
  _ , bsmminus = h2curve((edges, bsmpredminus))

  plt = fig.add_subplot(3, 1, (1, 2))

  plt.plot(xs, smplus, color="blue", label=r"$W^+ \to \ell^+ \nu_\ell$", lw=2)
  plt.plot(xs, smminus, color="orange", label=r"$W^- \to \ell^- \bar \nu_\ell$", lw=2)
  plt.plot(xs, bsmplus, color="blue", label=r"$W^+ \to \ell^+ N$", lw=2, ls=":")
  plt.plot(xs, bsmminus, color="orange", label=r"$W^- \to \ell^- \bar N$", lw=2, ls=":")
  plt.set_xticks([])
  plt.set_ylabel("events / GeV")

  ratioplt = fig.add_subplot(3, 1, 3)
  ratioplt.plot(xs, bsmplus/smplus, color="blue", lw=2, ls=":")
  ratioplt.plot(xs, bsmminus/smminus, color="orange", lw=2, ls=":")
  ratioplt.set_xlabel(axisdict[k])
  ratioplt.set_ylabel("ratio")

  plt.legend(title="ATLAS electron channel")

  fig.tight_layout()
  fig.savefig("ATLAS-%s.png" % k)
  fig.savefig("ATLAS-%s.pdf" % k)
  fig.clf()
