import yoda
import numpy
from matplotlib.figure import Figure


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

# lumi in pb
lumi = \
  { "CDF" : 8800
  , "ATLAS" : 4700
  , "CMS" : 16800
  }

# all in pb
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

  , "CMS" :
    { "SM+" : 5666
    , "SM-" : 3601
    , "BSM+" : 1053
    , "BSM-" : 1123
    }
  }

procs = \
  { "CDF" :
    { "SM+" : gethists(yoda.read("yoda/ppeplusvestable_cdf/Rivet_merged.yoda"))
    , "SM-" : gethists(yoda.read("yoda/ppeminusvestable_cdf/Rivet_merged.yoda"))
    , "BSM+" : gethists(yoda.read("yoda/ppeplusn1stable_cdf/Rivet_merged.yoda"))
    , "BSM-" : gethists(yoda.read("yoda/ppeminusn1stable_cdf/Rivet_merged.yoda"))
    }

  , "ATLAS" :
    { "SM+" : gethists(yoda.read("yoda/ppeplusvestable/Rivet_merged.yoda"))
    , "SM-" : gethists(yoda.read("yoda/ppeminusvestable/Rivet_merged.yoda"))
    , "BSM+" : gethists(yoda.read("yoda/ppeplusn1stable/Rivet_merged.yoda"))
    , "BSM-" : gethists(yoda.read("yoda/ppeminusn1stable/Rivet_merged.yoda"))
    }
  
  , "CMS" :
    { "SM+" : gethists(yoda.read("yoda/Rivet_SMW+_CMS.yoda"))
    , "SM-" : gethists(yoda.read("yoda/Rivet_SMW-_CMS.yoda"))
    , "BSM+" : gethists(yoda.read("yoda/Rivet_BSMW+_CMS.yoda"))
    , "BSM-" : gethists(yoda.read("yoda/Rivet_BSMW-_CMS.yoda"))
    }
  }

axisdict = \
  { "pTl" : r"$p_T^\ell$ [GeV]"
  , "mT" : r"$m_T$ [GeV]"
  }


fig = Figure((5, 5))
for experiment in ["CMS" , "CDF" , "ATLAS"]:
  for kin in ["pTl", "mT"]:
    edges , smplus = procs[experiment]["SM+"][kin]
    _ , smminus = procs[experiment]["SM-"][kin]
    _ , bsmplus = procs[experiment]["BSM+"][kin]
    _ , bsmminus = procs[experiment]["BSM-"][kin]

    widths = edges[1:] - edges[:-1]

    smpredplus = smplus * xsecs[experiment]["SM+"] * lumi[experiment]
    smpredminus = smminus * xsecs[experiment]["SM-"] * lumi[experiment]
    smpred = smpredplus + smpredminus

    bsmpredplus = bsmplus * xsecs[experiment]["BSM+"] * lumi[experiment]
    bsmpredminus = bsmminus * xsecs[experiment]["BSM-"] * lumi[experiment]
    bsmpred = bsmpredplus + bsmpredminus

    print("%s total SM+ prediction: %.2e events" % (experiment , (smpredplus * widths).sum()))
    print("%s total SM- prediction: %.2e events" % (experiment , (smpredminus * widths).sum()))
    print("%s total BSM+ prediction: %.2e events" % (experiment , (smpredplus * widths).sum()))
    print("%s total BSM- prediction: %.2e events" % (experiment , (smpredminus * widths).sum()))
    print("%s total SM prediction: %.2e events" % (experiment , (smpred * widths).sum()))
    print("%s total BSM prediction: %.2e events" % (experiment , (bsmpred * widths).sum()))

    xs , sm = h2curve((edges, smpred))
    _ , bsm = h2curve((edges, bsmpred))

    plt = fig.add_subplot(3, 1, (1, 2))

    plt.plot(xs, sm, color="blue", label=r"$W \to \ell \nu$", lw=2)
    plt.plot(xs, bsm, color="orange", label=r"$W \to \ell N$ ", lw=2, ls=":")
    plt.set_xticks([])
    plt.set_ylabel("events / GeV")

    ratioplt = fig.add_subplot(3, 1, 3)
    ratioplt.plot(xs, bsm/sm, color="orange", lw=2, ls=":")
    ratioplt.set_xlabel(axisdict[kin])
    ratioplt.set_ylabel("ratio")

    plt.legend(title="%s electron channel" % experiment)

    fig.tight_layout()
    fig.savefig("%s-%s-total.png" % (experiment , kin))
    fig.savefig("%s-%s-total.pdf" % (experiment , kin))
    fig.clf()


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
    ratioplt.set_xlabel(axisdict[kin])
    ratioplt.set_ylabel("ratio")

    plt.legend(title="%s electron channel" % experiment)

    fig.tight_layout()
    fig.savefig("%s-%s.png" % (experiment , kin))
    fig.savefig("%s-%s.pdf" % (experiment , kin))
    fig.clf()
