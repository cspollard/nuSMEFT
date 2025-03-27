import numpy
from matplotlib.figure import Figure

CRANGE = (-0.5, 2)

ATLASmTCs = [0, 934]
ATLASpTlCs = [0, 714]

CMSmTCs = [0, 529, -94.8]
CMSpTlCs = [0, 555, -132.1]

CDFmTCs = [0, 67.9, -16.6, 3.3]
CDFpTlCs = [0, 58.1, -14.3, 2.51]

# https://arxiv.org/abs/1701.07240
# assume ~all sensitivity is in pTl and there's no strong |\eta_\ell| dependence
cvATLASpTl = -29.2
uncertATLASpTl = 28

cvCDFmT = 80433.5 - 80357
uncertCDFmT = 9.4

cvCMSpTl = 57
uncertCMSpTl = 30

lim95 = 3.84
lim68 = 1.0

def poly(cs):
  def f(xs):
    return \
      sum \
      ( [ cs[i] * xs**i for i in range(len(cs)) ]
      , start=numpy.zeros_like(xs)
      )

  return f

curveATLASmT = poly(ATLASmTCs)
curveATLASpTl = poly(ATLASpTlCs)

curveCMSmT = poly(ATLASmTCs)
curveCMSpTl = poly(ATLASpTlCs)

curveCDFmT = poly(CDFmTCs)
curveCDFpTl = poly(CDFpTlCs)


def calcllh(calib, cs, cv, uncert):
  mus = cs * cs
  diffs = calib(mus) - cv
  return - diffs**2 / 2.0 / uncert**2


def dmax(llhs):
  return 2 * (llhs.max() - llhs)


def combine(llhs):
  return sum(llhs)


def lims(lim, xs, ys):
  allowed = xs[ys < lim]
  return numpy.min(allowed) , numpy.max(allowed)


def plotllhs(plt, cs, llhs, labels, colors, lss):
  plt.set_xlabel(r"$c$")
  plt.set_ylabel(r"$-2\Delta \log \mathcal{L}$")

  # add lines at \Delta log likelihood = 0.5, 1
  plt.plot([-100, 100], [lim68, lim68], lw=1, color=("gray", 0.5))
  plt.plot([-100, 100], [lim95 , lim95], lw=1, color=("gray", 0.5))


  for (llh , label , color , ls) in zip(llhs, labels, colors, lss):
    dllh = dmax(llh)

    print(label)
    print("95%% CL interval: [%0.3f, %0.3f]" % lims(lim95, cs, dllh) )
    print("68%% CL interval: [%0.3f, %0.3f]" % lims(lim68, cs, dllh) )
    print("central value: %0.3f" % cs[numpy.argmin(dllh)] )
    print()
    print()

    plt.plot(cs, dmax(llh), label=label, color=color, lw=2, ls=ls)

  plt.set_ylim(0, 4)

  return plt



fig = Figure((6, 4))
plt = fig.add_subplot()

cs = numpy.arange(CRANGE[0], CRANGE[1], 0.00001)

ATLASexp = calcllh(curveATLASpTl, cs, 0, uncertATLASpTl)
CMSexp = calcllh(curveCMSpTl, cs, 0, uncertCMSpTl)
LHCexp = combine([ATLASexp, CMSexp])

# expected
expllhs = \
  [ ATLASexp
  , CMSexp
  , LHCexp
  , calcllh(curveCDFmT, cs, 0, uncertCDFmT)
  ]


ATLASobs = calcllh(curveATLASpTl, cs, cvATLASpTl, uncertATLASpTl)
CMSobs = calcllh(curveCMSpTl, cs, cvCMSpTl, uncertCMSpTl)
LHCobs = combine([ATLASobs, CMSobs])

# observed
obsllhs = \
  [ ATLASobs
  , CMSobs
  , LHCobs
  , calcllh(curveCDFmT, cs, cvCDFmT, uncertCDFmT)
  ]

labels = \
  [ r"ATLAS: $m_W$ via $p_T^\ell$"
  , r"CMS: $m_W$ via $p_T^\ell$"
  , "LHC combination"
  , r"CDF: $m_W$ via $m_T^W$"
  ]

colors = ["blue", "orange", "black", "red"]

linestyles = ["--", "--", ":", "-"]

print("observed:")
print()
plt = \
  plotllhs \
  ( plt
  , cs
  , obsllhs
  , labels
  , colors
  , linestyles
  )

plt.set_xlim(CRANGE[0], CRANGE[1])

plt.legend(title="observed")

fig.tight_layout()
fig.savefig("llhs-observed.pdf")
fig.clf()

print("expected:")
print()

plt = fig.add_subplot()
plt = \
  plotllhs \
  ( plt
  , cs
  , expllhs
  , labels
  , colors
  , linestyles
  )

plt.set_xlim(CRANGE[0], CRANGE[1])

plt.legend(title="expected")

fig.tight_layout()
fig.savefig("llhs-expected.pdf")
