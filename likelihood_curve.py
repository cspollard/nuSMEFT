# import yoda
import numpy
from matplotlib.figure import Figure

# TODO
# all of these are ATLAS results so far?
# mean xsec
mu_sigma_p = 2747
atlas_theory_p = 2850
mu_sigma_m = 1964
atlas_theory_m = 1918

# error on xsec
alpha_sigma_p = numpy.sqrt(1**2 + 15**2 + 53**2 + 82**2)
alpha_sigma_m = numpy.sqrt(1**2 + 11**2 + 35**2 + 57**2)

# from where?
# mean mw
mu_mw = -29
# error on mw
alpha_mw = 28

# maybe Sam was normalizing + and - to the ATLAS predictions
# but doesn't this neglect the acceptance effects in the BSM scenario?
# xsec for mw as a function of bsm fraction
sigma_sm_p = 1767.8
sigma_bsm_p = 262.6
sigma_sm_m = 1103
sigma_bsm_m = 216.8

# mu = a * mW from a linear fit
a = 750


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
mus = numpy.arange(-0.1,0.1,0.00001)

# add and remove function calls inside numpy.append to change which curves are
# being fitted for.

def withmax(ll):
    return ll.max() - ll

llsigplus = loglikelihood_sigma_plus(mus) 
llsigminus = loglikelihood_sigma_minus(mus) 
llmass = loglikelihood_mass(mus)

lls = llsigplus + llsigminus + llmass

print("max log-likelihood: " + str(lls.max()))
llsigplus = withmax(llsigplus)
llsigminus = withmax(llsigminus)
llmass = withmax(llmass)
lls = withmax(lls)

# fbest = mus[lls.argmin()]
# print("max log-likelihood fraction: " + str(fbest))
# print("1-sigma down: " + str(mus[numpy.abs(lls[:lls.argmax()] + 1/2 ).argmin()]))
# print("1-sigma up: " + str(mus[lls.argmax() + numpy.abs(lls[lls.argmax():] + 1/2 ).argmin()]))
# print("2-sigma down: " + str( mus[ numpy.abs( lls[:lls.argmax()] + 1 ).argmin() ] ))
# print("2-sigma up: " + str( mus[lls.argmax() + numpy.abs( lls[lls.argmax():] + 1 ).argmin() ] ))

fig = Figure((6, 4))
plt = fig.add_subplot()

plt.set_xlabel(r"$\mu$")
plt.set_ylabel(r"$\Delta \log \mathcal{L}$")

# add lines at \Delta log likelihood = 0.5, 1
plt.plot([-100, 100], [0.5, 0.5], lw=1, color=("gray", 0.5))
plt.plot([-100, 100], [1, 1], lw=1, color=("gray", 0.5))

plt.plot(mus, lls, label="combined", color="black", lw=2)
plt.plot(mus, llmass, label="$m_W$", color="orange", lw=2, ls=":")
plt.plot(mus, llsigplus, label=r"$\sigma_{fid}(W^+)$", color="blue", lw=2, ls="--")
plt.plot(mus, llsigminus, label=r"$\sigma_{fid}(W^-)$", color="red", lw=2, ls="--")

plt.set_ylim(0, 2)
plt.set_xlim(-0.1, 0.1)

plt.legend(title="ATLAS")

fig.tight_layout()
fig.savefig("ATLASllh.pdf")
