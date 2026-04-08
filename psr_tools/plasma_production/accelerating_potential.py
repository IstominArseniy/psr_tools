import numpy as np
from ..interface import RadioPulsar

def RS_potential(PSR:RadioPulsar, x, phi=0):
    return PSR.hRS(x, phi)**2 * 2 * np.pi * PSR.rhoGJ()

def vac_potential(PSR:RadioPulsar, x, phi=0):
    return np.pi * PSR.rhoGJ() * PSR.R0**2 * np.cos(PSR.chi) * (1-x**2)

def simple_potential(PSR:RadioPulsar, x, phi=0):
    return np.minimum(RS_potential(PSR, x, phi), vac_potential(PSR, x, phi))

def psi_to_gamma(psi):
    return 5.9e-4 * psi
