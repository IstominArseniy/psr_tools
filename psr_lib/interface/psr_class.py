import numpy as np
from ..utils import constants
from scipy import integrate
from ..utils import utils_funcs as ufunc

"""
Class with basic information about RadioPulsar
"""

class RadioPulsar():
    def __init__(self, PSR_NAME, P, B12, chi_deg, Rkm=12, M_solar=1.4, Ir=100):
        # ----------------Set basic parameters------------------------
        self.PSR_NAME = PSR_NAME
        self.P = P # in s
        self.B_surf12 = B12
        self.B_surf = B12 * 1e12
        self.chi_deg = chi_deg
        self.Rkm = Rkm
        self.M_solar = M_solar # in solar masses
        self.Ir = Ir # in Solar masses * km^2
        # ----------------Set derived parameters------------------------
        self.chi = chi_deg / 180 * np.pi
        # B_type options: ATNF, BGI, MHD, USER
        #-************************
        # if B_model_type == 'ATNF':
        #     self.B0_12 = P**0.5 * P_DOT**0.5
        # elif B_model_type == 'BGI':
        #     self.B0_12 = min((P**0.5 * P_DOT**0.5)/np.cos(self.CHI)**2, 10)
        # elif B_model_type == 'MHD':
        #     self.B0_12 = P**0.5 * P_DOT**0.5 / (1 + np.sin(self.CHI)**2)
        # elif B_model_type == 'USER':
        #     self.B0_12 = USER_B
        # else:
        #     print('INCORRECT EVOLUTION MODEL')
        #*********************************
        self.Lambda = 20 # free pass length Lambda parameter
        self.R = self.Rkm * 1e5
        #----------gravitational corretions----
        self.epsGR = 3 * self.M_solar / self.R_km # gravitational corrections parametre
        self.KRc = 1 - 1/2 * self.epsGR        # General Reltivistic curvature radius (R_curv) correction
        self.KB = 1 + 3/4 * self.epsGR        # General Reltivistic magnetic field corrction
        self.Kpsi = (1 - 0.24*(self.Ir/100) / (self.R_km/12)**3) / (1 - self.epsGR)      # General Reltivistic rho_GJ correction
        self.sR0 = 1 - 3/8 * self.epsGR        # General Reltivistic polar cap radius correction
        #---------------------------------------
        self.Omega = 2 * np.pi / self.P # pulsar angular velocity
        self.R0 = self.sR0 * 1.25 * np.sqrt(1 + 0.2*(np.sin(self.chi))**2) * self.R * (self.Omega * self.R / constants.c)**0.5   # polar cap radius in cm
        self.OmegaB = constants.qe * self.B_surf / constants.me / constants.c # synchrotron frequency on the polar cap surface in s^-1
        self.RLC = constants.c / self.Omega # light cylinder radius in cm


    def read_from_file(self, filename): # probabably make as an additional init
        pass

    def write_to_file(self, filename):
        pass

    def Rc(self, r, units='cm'):
        """
        curvature radius
        """
        EPS = self.R**2 / self.R0 / self.RLC
        rc_cm =  self.KRc * 4/3 * self.R**2 / self.R0 / (r + EPS)
        if units == 'cm':
            return rc_cm
        elif units == 'R0':
            return rc_cm / self.R0
        elif units == 'R':
            return rc_cm / self.R
        else:
            print('Incorrect unit') # TODO make proper Raise Exeption

    def hRS(self, r, phi, units='cm'):
        """
        Ruderman-Sutherlend height
        units can be cm and R0 and R
        """
        cs=np.abs(np.cos(self.chi-1.5 * np.sin(phi) * np.sin(self.chi) * r * self.R0/ self.R ))
        h_cm = 1.1e2 * cs**(-3/7) * self.Rc(r, units='cm')**(2/7) * self.P**(3/7) / self.B_surf12**(4/7)
        if units == 'cm':
            return h_cm
        elif units == 'R0':
            return h_cm / self.R0
        elif units == 'R':
            return h_cm / self.R
        else:
            print('Incorrect unit') # TODO make proper Raise Exeption
            

class ObservedRadioPulsar(RadioPulsar):
    def __init__(self, PSR_NAME, P, Pdot, B12, chi_deg, beta_deg, freq=600, Rkm=12, M_solar=1.4, Ir=100):
        super().__init__(PSR_NAME, P, B12, chi_deg, Rkm, M_solar, Ir)
        self.Pdot = Pdot # in 1e-15
        self.beta_deg = beta_deg # impact angle in degrees
        self.beta = self.beta_def * np.pi / 180
        self.freq = freq # observation frequency

    def set_model_B(self, model):
        if model == 'ATNF':
            self.B_surf12 = self.P**0.5 * self.Pdot**0.5
        elif model == 'BGI':
            self.B_surf12 = min((self.P**0.5 * self.Pdot**0.5)/np.cos(self.chi)**2, 10)
        elif model == 'MHD':
            self.B_surf12 = self.P**0.5 * self.Pdot**0.5 / (1 + np.sin(self.chi)**2)
        else:
            print('INCORRECT EVOLUTION MODEL') # TODO Proper raise exception


class PSRutils:
    def __init__(self, PSR):
        self.PSR = PSR
    
    def B12(self, r):
        """
        magnetic field model
        """
        return self.PSR.B_surf12 * (1/r)**3

    def psi_inf(self, r, x):
        return r/self.PSR.Rc(r, x, units='R')

    def w(self, re, xe, r, E_ph):
        """
        pair production probability
        E_ph should be normalized to m_e c^2
        return: probability density (in w dl, dl expressed in R star)
        """
        if (r - re) / re < 1e-15:
            return 0
        psi = self.psi_inf(re, xe) * (1 - re/r)
        return 0.23 * 1 / 137 / (constants.e_lambda_bar/self.PSR.R) * self.B12(r) / (constants.Bcr12) * psi * np.exp(-8/3 / E_ph * constants.Bcr12 / self.B12(r) / psi)

    def E_ph_min_exact(self, re, xe):
        """
        binary search computation of the minimal photon energy required to pair produce
        """
        E1 = 10
        E2 = 1e8
        while((E2 - E1) / E2 > 1e-3):
            E = (E1 + E2)/2
            tau_inf = integrate.quad(lambda h: self.w(re, xe, h, E), re, 100)[0]
            if tau_inf > 1:
                E2 = E
            else:
                E1 = E
        return E 
    
    def E_c(self, r, x, gamma_e):
        """
        curvature radiation characteristic energy
        """
        return 3 / 2 * (constants.e_lambda_bar / self.PSR.R) / self.PSR.Rc(r, x) * gamma_e**3

    def n_curv(self, r, x, E_ph, gamma_e):
        """
        curvature radiation spectrum
        dNph = n_curv * dE_ph * dh
        """
        if callable(gamma_e): # if gamma_e is a function of height (gamma = gamma(h))
            return np.sqrt(3) / 2 / np.pi * constants.alpha_e / self.PSR.Rc(r, x) * gamma_e(r) * ufunc.F_curavture(E_ph/self.E_c(r, x, gamma_e(r))) / E_ph
        else:
            return np.sqrt(3) / 2 / np.pi * constants.alpha_e / self.PSR.Rc(r, x) * gamma_e * ufunc.F_curavture(E_ph/self.E_c(r, x, gamma_e)) / E_ph
        