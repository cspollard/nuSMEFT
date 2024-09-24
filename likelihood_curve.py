import yoda
import numpy as np
import matplotlib.pyplot as plt

# mean xsec
mu_sigma_p = 2747
atlas_theory_p = 2850
mu_sigma_m = 1964
atlas_theory_m = 1918

# error on xsec
alpha_sigma_p = np.sqrt(1**2 + 15**2 + 53**2 + 82**2)
alpha_sigma_m = np.sqrt(1**2 + 11**2 + 35**2 + 57**2)

# mean mw
mu_mw = -29
# error on mw
alpha_mw = 28

# TODO: these look incorrect? need to ask Sam
# these must be fiducial cross sections.
# also why are these symmetric?
# maybe Sam was normalizing + and - to the ATLAS predictions
# but doesn't this neglect the acceptance effects in the BSM scenario?
# xsec for mw as a function of bsm fraction
sigma_sm_p = 1103
sigma_bsm_p = 216.7
sigma_sm_m = 1103
sigma_bsm_m = 216.7

# mu = a * mW from a linear fit
a = 750


# why not ratio of \sigma+ / \sigma-?
# Chris is a bit confused by this to be honest...


def loglikelihood_sigma_plus(mus):
    return \
        - ( (atlas_theory_p/sigma_sm_p) \
           * (sigma_sm_p + mus * sigma_bsm_p) \
           - atlas_theory_p \
            )**2 \
        / (2*alpha_sigma_p**2)


def loglikelihood_sigma_minus(mus):
    return \
        - ( (atlas_theory_m/sigma_sm_m) \
            * (sigma_sm_m + mus * sigma_bsm_m) - atlas_theory_m)**2 \
        / (2*alpha_sigma_m**2)

def loglikelihood_mass(mus):
    return - (a*mus)**2 / (2*alpha_mw**2)
    
# likelihoods
lls = np.array([])

# bsm fractions
mus = np.arange(-0.1,0.1,0.00001)

# add and remove function calls inside np.append to change which curves are being fitted for.
lls = loglikelihood_sigma_plus(mus) + loglikelihood_sigma_minus(mus) + loglikelihood_mass(mus)

print("max log-likelihood: " + str(lls.max()))
lls -= lls.max()

fbest = mus[lls.argmax()]
print("max log-likelihood fraction: " + str(fbest))
print("1-sigma down: " + str(mus[np.abs(lls[:lls.argmax()] + 1/2 ).argmin()]))
print("1-sigma up: " + str(mus[lls.argmax() + np.abs(lls[lls.argmax():] + 1/2 ).argmin()]))
print("2-sigma down: " + str( mus[ np.abs( lls[:lls.argmax()] + 1 ).argmin() ] ))
print("2-sigma up: " + str( mus[lls.argmax() + np.abs( lls[lls.argmax():] + 1 ).argmin() ] ))

plt.xlabel("Signal Strength")
plt.ylabel(r"$\Delta \ln \mathcal{L}$")
plt.title("Joint Cross Section And Mass Likelihood")

plt.plot(mus,lls)
plt.savefig("likelihood_fits/atlas_joint_fit_sigma_mw.pdf")