import numpy
from matplotlib.figure import Figure

MURANGE = (-0.5, 2)

def poly(cs):
  def f(xs):
    return \
      sum \
      ( [ cs[i] * xs**i for i in range(len(cs)) ]
      , start=numpy.zeros_like(xs)
      )

  return f

curveATLASmT = poly([0, 934])
curveATLASpTl = poly([0, 714])

# https://arxiv.org/abs/1701.07240
# assume ~all sensitivity is in pTl and there's no strong |\eta_\ell| dependence
cvATLASpTl = -29.2
uncertATLASpTl = 28

curveCDFmT = poly([0, 67.9, -16.6, 3.3])
curveCDFpTl = poly([0, 58.1, -14.3, 2.51])

cvCDFmT = 80433.5 - 80357
uncertCDFmT = 9.4

cvCMSpTl = 57
uncertCMSpTl = 30


def calcllh(calib, mus, cv, uncert):
  diffs = calib(mus) - cv
  return - diffs**2 / 2.0 / uncert**2

def dmax(llhs):
  return llhs.max() - llhs


def plotllhs(plt, mus, llhs, labels, colors, lss):
  plt.set_xlabel(r"$\mu$")
  plt.set_ylabel(r"$\Delta \log \mathcal{L}$")

  # add lines at \Delta log likelihood = 0.5, 1
  plt.plot([-100, 100], [0.5, 0.5], lw=1, color=("gray", 0.5))
  plt.plot([-100, 100], [1, 1], lw=1, color=("gray", 0.5))

  for (llh , label , color , ls) in zip(llhs, labels, colors, lss):
    plt.plot(mus, dmax(llh), label=label, color=color, lw=2, ls=ls)

  plt.set_ylim(0, 2)

  return plt

def combine(llhs):
  return sum(llhs)


fig = Figure((6, 4))
plt = fig.add_subplot()

mus = numpy.arange(MURANGE[0], MURANGE[1], 0.00001)

ATLASexp = calcllh(curveATLASpTl, mus, 0, uncertATLASpTl)
CMSexp = calcllh(curveATLASpTl, mus, 0, uncertCMSpTl)

# expected
expllhs = \
  [ ATLASexp
  , CMSexp
  , combine([ATLASexp, CMSexp])
  , calcllh(curveCDFmT, mus, 0, uncertCDFmT)
  ]


ATLASobs = calcllh(curveATLASpTl, mus, cvATLASpTl, uncertATLASpTl)
CMSobs = calcllh(curveATLASpTl, mus, cvCMSpTl, uncertCMSpTl)

# observed
obsllhs = \
  [ ATLASobs
  , CMSobs
  , combine([ATLASobs, CMSobs])
  , calcllh(curveCDFmT, mus, cvCDFmT, uncertCDFmT)
  ]

labels = \
  [ r"ATLAS: $m_W$ via $p_T^\ell$"
  , r"CMS: $m_W$ via $p_T^\ell$"
  , "LHC combination"
  , r"CDF: $m_W$ via $m_T^W$"
  ]

colors = ["blue", "orange", "black", "red"]

linestyles = ["--", "--", ":", "-"]

plt = \
  plotllhs \
  ( plt
  , mus
  , obsllhs
  , labels
  , colors
  , linestyles
  )

plt.set_xlim(MURANGE[0], MURANGE[1])

plt.legend(title="observed")

fig.tight_layout()
fig.savefig("llhs-observed.pdf")
fig.clf()


plt = fig.add_subplot()
plt = \
  plotllhs \
  ( plt
  , mus
  , expllhs
  , labels
  , colors
  , linestyles
  )

plt.set_xlim(MURANGE[0], MURANGE[1])

plt.legend(title="expected")

fig.tight_layout()
fig.savefig("llhs-expected.pdf")
