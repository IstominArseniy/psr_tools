import numpy as np
import itertools

from multiprocessing import Pool
from psr_tools.radiative_trasfer import loki_python_binding
from psr_tools.interface import PulsarProfile


class FixedHeightModel:
    def __init__(self, multiplicity, gamma, Rem, fr=1, fphi=1):
        self.fr = fr
        self.fphi = fphi
        self.multiplicity = multiplicity
        self.gamma = gamma
        self.Rem = Rem


class ProfileCalculator:
    def __init__(self, PSR, model):
        self.model = model
        self.PSR = PSR
        self.psr_dict = {"B12":PSR.B_surf12, "Period":PSR.P, "chi_deg":PSR.chi_deg, "beta_deg":PSR.beta_deg, "freqGHz":PSR.freq*1e-3, "Rs":PSR.R}
        self.model_dict = {"fr":model.fr, "fphi":model.fphi, "Rem":model.Rem, "lambda":model.multiplicity, "gamma0":model.gamma}
        self.profile_calculator = loki_python_binding.ProfileCalculator(self.psr_dict, self.model_dict)

    @staticmethod
    def _worker(calculator, phi, mode):
        return calculator.find_ILVPA(phi, mode)
    
    def find_ILVPA(self, phi, mode, with_absorption=True):
        return self.profile_calculator.find_ILVPA(phi, mode, with_absorption)
    
    def find_I(self, phi, with_absorption=True):
        return self.profile_calculator.find_I(phi, 0, with_absorption)

    def calculate_profile(self, phi_start, phi_end, phi_step, mode, with_absroption=True):
        phi_arr = np.arange(phi_start, phi_end, phi_step)
        mode_arr = np.ones(phi_arr.shape[0], dtype=np.int8) * mode
        res = []
        # with Pool(10) as pool:
            # res = pool.starmap(ProfileCalculator._worker, zip(itertools.repeat(self.profile_calculator), phi_arr, mode_arr))
        res = list(map(lambda phi: self.profile_calculator.find_ILVPA(phi, mode, with_absroption), phi_arr))
        Is = []
        Ls = []
        Vs = []
        PAs = []
        for ILVPA in res:
            Is.append(ILVPA['I'])
            Ls.append(ILVPA['L'])
            Vs.append(ILVPA['V'])
            PAs.append(ILVPA['PA'])
        profile = PulsarProfile.from_ILVPA(Is, Ls, Vs, PAs, phi_step, normalize=False)
        return profile
    
    def get_r_esc(self):
        return 1e3*(self.model.multiplicity/1e4)**0.4 * (self.model.gamma/1e2)**(-1.2)*self.PSR.B_surf12**0.4*(self.PSR.freq*1e-3)**(-0.4)*self.PSR.P**(-0.2)
                
    def calculate_intencity_profile(self, phi_start, phi_end, phi_step, with_absorption=True):
        phi_arr = np.arange(phi_start, phi_end, phi_step)
        res = list(map(lambda phi: self.find_I(phi, with_absorption), phi_arr))
        Is = []
        for I in res:
            Is.append(I)
        return np.array(Is)
    