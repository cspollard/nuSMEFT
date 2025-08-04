# import yoda
import numpy
from matplotlib.figure import Figure

# TODO
# all of these are ATLAS results so far?
# mean xsec
# https://arxiv.org/abs/1612.03016
mu_sigma_p = 2947
atlas_theory_p = 2850
mu_sigma_m = 1964
atlas_theory_m = 1918

# error on xsecs
alpha_sigma_p = numpy.sqrt(1**2 + 15**2 + 53**2 + 82**2)
alpha_sigma_m = numpy.sqrt(1**2 + 11**2 + 35**2 + 57**2)

# mW+ - mW-
mu_mw = -29
# error on mw+ - mW-
alpha_mw = 28

# maybe Sam was normalizing + and - to the ATLAS predictions
# but doesn't this neglect the acceptance effects in the BSM scenario?
# xsec for mw as a function of bsm fraction
sigma_sm_p = 1767.8
sigma_bsm_p = 262.6
sigma_sm_m = 1103
sigma_bsm_m = 216.8

# mu = a * mW from a linear fit
a = 714


def loglikelihood_sigma_plus(mus):
    return \
        - ( (atlas_theory_p / sigma_sm_p) * (sigma_sm_p + mus * sigma_bsm_p) \
            - atlas_theory_p \
          )**2 \
        / (2*alpha_sigma_p**2)


def loglikelihood_sigma_minus(mus):
    return \
        - ( (atlas_theory_m / sigma_sm_m) * (sigma_sm_m + mus * sigma_bsm_m)
            - atlas_theory_m
          )**2 \
        / (2*alpha_sigma_m**2)


def loglikelihood_mass(mus):
    return - (a*mus)**2 / (2*alpha_mw**2)


# bsm fractions
cs = numpy.arange(0, 1, 0.001)
c2s = cs * cs

# add and remove function calls inside numpy.append to change which curves are
# being fitted for.

def withmax(ll):
    return ll.max() - ll

llsigplus = loglikelihood_sigma_plus(c2s)
llsigminus = loglikelihood_sigma_minus(c2s)
llmass = loglikelihood_mass(c2s)

lls = llsigplus + llsigminus + llmass

print("max log-likelihood: " + str(lls.max()))
llsigplus = 2 * withmax(llsigplus)
llsigminus = 2 * withmax(llsigminus)
llmass = 2 * withmax(llmass)
lls = 2 * withmax(lls)

# fbest = cs[lls.argmin()]
# print("max log-likelihood fraction: " + str(fbest))
# print("1-sigma down: " + str(cs[numpy.abs(lls[:lls.argmax()] + 1/2 ).argmin()]))
# print("1-sigma up: " + str(cs[lls.argmax() + numpy.abs(lls[lls.argmax():] + 1/2 ).argmin()]))
# print("2-sigma down: " + str( cs[ numpy.abs( lls[:lls.argmax()] + 1 ).argmin() ] ))
# print("2-sigma up: " + str( cs[lls.argmax() + numpy.abs( lls[lls.argmax():] + 1 ).argmin() ] ))

fig = Figure((6, 4))
plt = fig.add_subplot()

plt.set_xlabel(r"$|c_{HNe}|$")
plt.set_ylabel(r"$-2 \Delta \log L$")

# add lines at \Delta log likelihood = 0.5, 1
plt.plot([-100, 100], [1, 1], lw=1, color=("gray", 0.5))
plt.plot([-100, 100], [3.84, 3.84], lw=1, color=("gray", 0.5))

plt.plot(cs, lls, label="combined", color="black", lw=2)
plt.plot(cs, llmass, label="$m_W$", color="orange", lw=2, ls=":")
plt.plot(cs, llsigplus, label=r"$\sigma_{fid}(W^+)$", color="blue", lw=2, ls="--")
plt.plot(cs, llsigminus, label=r"$\sigma_{fid}(W^-)$", color="red", lw=2, ls="--")

plt.set_ylim(0, 4)
plt.set_xlim(cs[0], cs[-1])

plt.legend(title="ATLAS")

fig.tight_layout()
fig.savefig("ATLASllh.pdf")
