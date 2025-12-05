import numpy as np
from ..utils import constants
from scipy import integrate
from ..utils import utils_funcs as ufunc

import tomlkit

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
        self.Lambda = 20 # free pass length Lambda parameter
        self.R = self.Rkm * 1e5
        #----------gravitational corretions----
        self.epsGR = 3 * self.M_solar / self.Rkm # gravitational corrections parametre
        self.KRc = 1 - 1/2 * self.epsGR        # General Reltivistic curvature radius (R_curv) correction
        self.KB = 1 + 3/4 * self.epsGR        # General Reltivistic magnetic field corrction
        self.Kpsi = (1 - 0.24*(self.Ir/100) / (self.Rkm/12)**3) / (1 - self.epsGR)      # General Reltivistic rho_GJ correction
        self.sR0 = 1 - 3/8 * self.epsGR        # General Reltivistic polar cap radius correction
        #---------------------------------------
        self.Omega = 2 * np.pi / self.P # pulsar angular velocity
        # self.R0 = self.sR0 * 1.25 * np.sqrt(1 + 0.2*(np.sin(self.chi))**2) * self.R * (self.Omega * self.R / constants.c)**0.5   # polar cap radius in cm
        self.R0 = self.R * (self.Omega * self.R / constants.c)**0.5   # polar cap radius in cm (simplified version)

        self.OmegaB = constants.qe * self.B_surf / constants.me / constants.c # synchrotron frequency on the polar cap surface in s^-1
        self.RLC = constants.c / self.Omega # light cylinder radius in cm
        self.val = 0.825

    @classmethod
    def from_file(cls, filename): 
        with open(filename, mode="rt", encoding="utf-8") as fb:
            contents = tomlkit.load(fb).unwrap()
            Name, P, B12, chi = contents['general']['Name'], contents['general']['P'], contents['general']['B12'], contents['general']['chi']
            Rkm, M_solar, Ir = contents['mechanics']['Rs'], contents['mechanics']['M'], contents['mechanics']['Ir']
        psr = cls(Name, P, B12, chi, Rkm, M_solar, Ir)
        return psr

    def write_to_file(self, filename=None):
        if filename is None:
            filename = self.PSR_NAME 
        toml_doc = tomlkit.document()
        general = tomlkit.table()
        general.add("Name", self.PSR_NAME)
        general["Name"].comment("Pulsar Name")
        general.add("P", self.P)
        general["P"].comment("Rotation period in s")
        general.add("B12", self.B_surf12)
        general["B12"].comment("Magnetic filed on the polar cap, nomralized to 1e12 G")
        general.add("chi", self.chi_deg)
        general["chi"].comment("inclination angle in degrees")
        mechanics = tomlkit.table()
        mechanics.add("Rs", self.Rkm)
        mechanics["Rs"].comment("Neutron star radius in km")
        mechanics.add("M", self.M_solar)
        mechanics["M"].comment("Neutron star mass in solar masses")
        mechanics.add("Ir", self.Ir)
        mechanics["Ir"].comment("Neutron star moment of inertia in Solar masses * km^2")
        toml_doc.add("general", general)
        toml_doc.add("mechanics", mechanics)
        with open(filename+".toml", "w") as fp:
            tomlkit.dump(toml_doc, fp)

    def get_half_opening_angle(self, Rem, units='rad'):
        angle = 3/2* np.sqrt(Rem) * self.R0/self.R
        if units == 'rad':
            return angle
        elif units == 'deg':
            return angle * 180 / np.pi
        else:
            raise ValueError("Incorrect units name.")



    def Rc(self, x, r, units='cm'):
        """
        x - position on polar cap (normalized to R0)
        r - distance from the star center (normalized to R)
        units can be cm, R0 and R
        return: curvature radius
        """
        EPS = self.R**2 / self.R0 / self.RLC
        EPS=0
        Rc_cm =  self.KRc * 4/3 * self.R**2 / self.R0 / (x + EPS) * np.sqrt(r)
        if units == 'cm':
            return Rc_cm
        elif units == 'R0':
            return Rc_cm / self.R0
        elif units == 'R':
            return Rc_cm / self.R
        else:
            raise ValueError("Incorrect units name.")

    def hRS(self, r, phi, units='cm'):
        """
        Ruderman-Sutherlend height
        units can be cm, R0 and R
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
    
    def rhoGJ(self, x=0, phi=0):
        return self.Omega *  self.B_surf * np.cos(self.chi) / 2 /np.pi / constants.c
    
    def w(self, x_em, r_em, r, Eph):
        """
        pair production probability
        Eph should be normalized to m_e c^2
        return: probability density (in w dl, dl expressed in R star)
        """
        if (r - r_em) / r_em < 1e-15:
            return 0
        psi = self.psi_inf(x_em, r_em) * (1 - r_em/r)
        return 0.23 * 1 / 137 / (constants.e_lambda_bar/self.R) * self.B12(r) / (constants.Bcr12) * psi * np.exp(-8/3 / Eph * constants.Bcr12 / self.B12(r) / psi)


    def B12(self, r):
        """
        magnetic field model
        """
        return self.B_surf12 / r**3
    
    def l_gamma(self, x, h, Eph, units='cm', approximation='coarse'):
        r = 1 + h * self.R0 / self.R
        if approximation == 'coarse':
            return 8 / 3 / self.Lambda * constants.Bcr12 / self.B12(r) / Eph  * self.Rc(x, r, units=units)
        elif approximation == 'fine':
            pass
        elif approximation == 'exact':
            pass
        else:
            raise ValueError("Incorrect approximation type.")
        
    def x_absorption(self, x, h, Eph, units='cm', approximation='coarse'):
        return x * (1 - 3/8 * self.l_gamma(x, h, Eph, units=units, approximation=approximation)**2) 
    
    def r_absorption(self, x, h, Eph, units='cm', approximation='coarse'):
        ra_cm = h * self.R0 + self.R + self.l_gamma(x, h, Eph, units='cm', approximation=approximation)
        if units == 'cm':
            return ra_cm
        elif units == 'R0':
            return ra_cm / self.R0
        elif units == 'R':
            return ra_cm / self.R
        else:
            raise ValueError("Incorrect units name.")
    
    def psi_inf(self, x, r):
        return r / self.Rc(x, r, units='R')
        
    def Eph_min_exact(self, x, r):
        """
        binary search computation of the minimal photon energy required to pair produce
        """
        E1 = 10
        E2 = 1e8
        while((E2 - E1) / E2 > 1e-3):
            E = (E1 + E2)/2
            tau_inf = integrate.quad(lambda h: self.w(x, r, h, E), r, 100)[0] # upper limmit was 100 - need to be tested
            if tau_inf > 1:
                E2 = E
            else:
                E1 = E
        return E 

    def Ec(self, x, r, gamma_e):
        """
        curvature radiation characteristic energy
        """
        return 3 / 2 * (constants.e_lambda_bar / self.R) / self.Rc(x, r, units='R') * gamma_e**3
    

    def curvature_emission_denisty(self, x, r, Eph, gamma_e):
        """
        dNph = n_curv * dEph * dh
        """
        if callable(gamma_e):
            return np.sqrt(3) / 2 / np.pi * constants.alpha_e / self.Rc(x, r, units='R') * gamma_e(r) * ufunc.F_curavture(Eph/self.Ec(x, r, gamma_e(r))) / Eph
        else:
            return np.sqrt(3) / 2 / np.pi * constants.alpha_e / self.Rc(x, r, units='R') * gamma_e * ufunc.F_curavture(Eph/self.Ec(x, r, gamma_e)) / Eph
        
    def curvature_multipliticy(self, x, gamma_e, tol=1e-2):
        return integrate.dblquad(lambda Eph, r: self.curvature_emission_denisty(x, r, Eph, gamma_e), 1, 100, lambda r: self.Eph_min_exact(x, r), np.inf, epsrel=1e-2)[0]

class ObservedRadioPulsar(RadioPulsar):
    def __init__(self, PSR_NAME, P, Pdot, B12, chi_deg, beta_deg, freq=600, Rkm=12, M_solar=1.4, Ir=100):
        super().__init__(PSR_NAME, P, B12, chi_deg, Rkm, M_solar, Ir)
        self.Pdot = Pdot # in 1e-15
        self.beta_deg = beta_deg # impact angle in degrees
        self.beta = self.beta_deg * np.pi / 180
        self.freq = freq # observation frequency in Mhz

    def set_beta(self, new_beta_deg):
        self.beta_deg = new_beta_deg
        self.beta = self.beta_deg*np.pi/180

    def set_chi(self, new_chi_deg):
        self.chi_deg = new_chi_deg
        self.chi = self.chi_deg*np.pi/180

    def set_model_B(self, model):
        if model == 'ATNF':
            self.B_surf12 = self.P**0.5 * self.Pdot**0.5
        elif model == 'BGI':
            self.B_surf12 = min((self.P**0.5 * self.Pdot**0.5)/np.cos(self.chi)**2, 10)
        elif model == 'MHD':
            self.B_surf12 = self.P**0.5 * self.Pdot**0.5 / (1 + np.sin(self.chi)**2)
        else:
            print('INCORRECT EVOLUTION MODEL') # TODO Proper raise exception

    @classmethod
    def from_file(cls, filename): 
        with open(filename, mode="rt", encoding="utf-8") as fb:
            contents = tomlkit.load(fb).unwrap()
            Name, P, Pdot, B12, chi, beta, freq = contents['general']['Name'], contents['general']['P'], contents['general']['Pdot'], \
            contents['general']['B12'], contents['general']['chi'], contents['general']['beta'], contents['general']['frequency']
            Rkm, M_solar, Ir = contents['mechanics']['Rs'], contents['mechanics']['M'], contents['mechanics']['Ir']
        psr = cls(Name, P, Pdot, B12, chi, beta, freq, Rkm, M_solar, Ir)
        return psr

    def write_to_file(self, filename=None): 
        if filename is None:
            filename = self.PSR_NAME 
        toml_doc = tomlkit.document()
        general = tomlkit.table()
        general.add("Name", self.PSR_NAME)
        general["Name"].comment("Pulsar Name")
        general.add("P", self.P)
        general["P"].comment("Rotation period in s")
        general.add("Pdot", self.Pdot)
        general["Pdot"].comment("Rotation period derivative in 1e-15")
        general.add("B12", self.B_surf12)
        general["B12"].comment("Magnetic filed on the polar cap, nomralized to 1e12 G")
        general.add("chi", self.chi_deg)
        general["chi"].comment("inclination angle in degrees")
        general.add("beta", self.beta_deg)
        general["beta"].comment("impact angle in degrees")
        general.add("frequency", self.freq)
        general["frequency"].comment("observational frequency in Mhz")
        mechanics = tomlkit.table()
        mechanics.add("Rs", self.Rkm)
        mechanics["Rs"].comment("Neutron star radius in km")
        mechanics.add("M", self.M_solar)
        mechanics["M"].comment("Neutron star mass in solar masses")
        mechanics.add("Ir", self.Ir)
        mechanics["Ir"].comment("Neutron star moment of inertia in Solar masses * km^2")
        toml_doc.add("general", general)
        toml_doc.add("mechanics", mechanics)
        with open(filename+".toml", "w") as fp:
            tomlkit.dump(toml_doc, fp)


