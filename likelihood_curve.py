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

# xsec for mw as a function of bsm fraction
sigma_sm_p = 1103
sigma_bsm_p = 216.7
sigma_sm_m = 1103
sigma_bsm_m = 216.7

# parameters for mw as a function of bsm fraction
a = 750
b = 0.92

# only defined for f>0
def loglikelihood(f):
    ll = - ( (sigma_sm_p + f*sigma_bsm_p)/sigma_sm_p - atlas_theory_p/atlas_theory_p)**2/(2*alpha_sigma_p**2) - ( (sigma_sm_m + f*sigma_bsm_m)/sigma_sm_m - atlas_theory_m/atlas_theory_m)**2/(2*alpha_sigma_m**2) - (a*f + b - mu_mw)**2/(2*alpha_mw**2)

    return ll

def loglikelihood_sigma_plus(f):
    ll = - ( (atlas_theory_p/sigma_sm_p)*(sigma_sm_p + f*sigma_bsm_p) - atlas_theory_p)**2/(2*alpha_sigma_p**2)

    return ll

def loglikelihood_sigma_minus(f):
    ll = - ( (atlas_theory_m/sigma_sm_m)*(sigma_sm_m + f*sigma_bsm_m) - atlas_theory_m)**2/(2*alpha_sigma_m**2)

    return ll

def loglikelihood_mass(f):
    ll = - (a*f + b)**2/(2*alpha_mw**2)

    return ll
    
# likelihoods
lls = np.array([])

# bsm fractions
fracs = np.arange(-0.5,0.5,0.00001)

# add and remove function calls inside np.append to change which curves are being fitted for.
# e.g. lls = np.append(lls, loglikelihood_sigma_plus(x) + loglikelihood_sigma_minus(x)) fits for cross section only.
for x in fracs:
    lls = np.append(lls, loglikelihood_sigma_plus(x) + loglikelihood_sigma_minus(x) + loglikelihood_mass(x))

print("max log-likelihood: " + str(lls.max()))
lls -= lls.max()

fbest = fracs[lls.argmax()]
print("max log-likelihood fraction: " + str(fbest))
print("1-sigma down: " + str(fracs[np.abs(lls[:lls.argmax()] + 1/2 ).argmin()]))
print("1-sigma up: " + str(fracs[lls.argmax() + np.abs(lls[lls.argmax():] + 1/2 ).argmin()]))
print("2-sigma down: " + str( fracs[ np.abs( lls[:lls.argmax()] + 1 ).argmin() ] ))
print("2-sigma up: " + str( fracs[lls.argmax() + np.abs( lls[lls.argmax():] + 1 ).argmin() ] ))

plt.xlabel("Signal Strength")
plt.ylabel(r"$\Delta \ln \mathcal{L}$")
plt.title("Joint Cross Section And Mass Likelihood")

plt.plot(fracs,lls)
plt.savefig("likelihood_fits/atlas_joint_fit_sigma_mw.pdf")