import numpy as np
from ..utils import constants

class Bmodel:
    def __init__(self, chi_deg, P):
        self.chi = chi_deg * np.pi/180
        self.P = P
        self.Omega = 2*np.pi/P
        self.Omega_vec = np.array((0, 0, self.Omega), dtype=np.float32)
        self.Rs = 1.2e6
        self.RLC = constants.c/self.Omega / self.Rs
        self.fr = 1
        self.fphi = 1
    
    def B_vec(self, vR, m)->np.array:
        R = np.linalg.norm(vR)
        n = vR/R
        Bdipole = 3 * m.dot(n) * n - m
        Bwind = np.zeros(3)
        Bwind[0] = R/self.RLC * self.fr * n[0] - (R/self.RLC)**2 * self.fphi * self.fr * (-n[1])
        Bwind[1] = R/self.RLC * self.fr * n[1] - (R/self.RLC)**2 * self.fphi * self.fr * n[0]
        Bwind[2] = R/self.RLC * self.fr * n[2]
        return Bdipole + Bwind
    
    def b_vec(self, vR, m)->np.array:
        return self.B_vec(vR, m) / np.linalg.norm(self.B_vec(vR, m))
    
    def BetaR(self, vR)->np.array:
        return self.Rs/constants.c*np.cross(self.Omega_vec, vR)
    
    
class Ray:
    def __init__(self, Bmodel, observer_vec, phase_deg, emission_height):
        self.observer_vec = observer_vec
        self.phase = phase_deg * np.pi/180
        self.emission_height = emission_height
        self.Bmodel = Bmodel
        self.theta_em, self.phi_em = self.find_emission_point()

    def _DX(function, x, y):
        dx = 1e-5
        return (function(x+dx, y) - function(x-dx, y))/2/dx
    
    def _DY(function, x, y):
        dy = 1e-5
        return (function(x, y+dy) - function(x, y-dy))/2/dy

    def find_emission_point(self):
        theta_em = self.Bmodel.chi
        phi_em = self.phase
        magnetic_moment_vec = self.vM(0)
        def func(theta, phi):
            vPoint = np.zeros(3)
            vPoint[0] = self.emission_height * np.sin(theta) * np.cos(phi)
            vPoint[1] = self.emission_height * np.sin(theta) * np.sin(phi)
            vPoint[2] = self.emission_height * np.cos(theta)
            return np.cross(self.Bmodel.B_vec(vPoint, magnetic_moment_vec), self.observer_vec)
        for i in range(5):
            f1x = Ray._DX(lambda phi, theta: func(phi, theta)[0], theta_em, phi_em)
            f2x = Ray._DX(lambda phi, theta: func(phi, theta)[1], theta_em, phi_em)
            f1y = Ray._DY(lambda phi, theta: func(phi, theta)[0], theta_em, phi_em)
            f2y = Ray._DY(lambda phi, theta: func(phi, theta)[1], theta_em, phi_em)
            f1 = func(theta_em, phi_em)[0]
            f2 = func(theta_em, phi_em)[1]
            dX = (f1y * f2 - f1 * f2y) / (f1x * f2y - f1y * f2x)
            dY = (f1x * f2 - f1 * f2x) / (f1y * f2x - f2y * f1x)
            theta_em += dX
            phi_em += dY
        return (theta_em, phi_em)

    def vR(self, l)->np.array:
        n0 = np.zeros(3)
        n0[0] = np.sin(self.theta_em) * np.cos(self.phi_em) 
        n0[1] = np.sin(self.theta_em) * np.sin(self.phi_em)
        n0[2] = np.cos(self.theta_em)
        return self.emission_height * n0 + l * self.observer_vec

    def vM(self, l)->np.array:
        magnetic_moment_vec = np.zeros(3)
        magnetic_moment_vec[0] = np.sin(self.Bmodel.chi) * np.cos(self.phase + l / self.Bmodel.RLC)
        magnetic_moment_vec[1] = np.sin(self.Bmodel.chi) * np.sin(self.phase + l / self.Bmodel.RLC)
        magnetic_moment_vec[2] = np.cos(self.Bmodel.chi)
        return magnetic_moment_vec

    def pis_m(self, l):
        pass

    def b_vec(self, l):
        return self.Bmodel.b_vec(self.vR(l), self.vM(l))

    def BetaR(self, l):
        return self.Bmodel.BetaR(self.vR(l))


    def q_perp_vec(self, l)->np.array:
        return -np.cross(self.observer_vec, np.cross(self.b_vec(l),self.BetaR(l)))

    def beta_eff(self, l):
        if l<1e-1:
            l=1e-1
        ex = self.Bmodel.Omega_vec - self.observer_vec.dot(self.Bmodel.Omega_vec) * self.observer_vec
        ex /= np.linalg.norm(ex)
        ey = np.cross(self.observer_vec, ex)
        PA_vec = self.b_vec(l) + self.q_perp_vec(l)
        PAx = PA_vec.dot(ex)
        PAy = PA_vec.dot(ey)
        return np.atan2(PAy, PAx)
    
    

class PAmodel:
    def __init__(self, chi_deg, beta_deg, emission_height, P):
        self.chi = chi_deg * np.pi / 180
        self.beta = beta_deg * np.pi/180
        self.dzeta = self.chi - self.beta
        self.emission_height = emission_height
        self.P = P    
        self.Bmodel = Bmodel(chi_deg, P)
        self.observer_vec = np.array((np.sin(self.dzeta), 0, np.cos(self.dzeta)))

    def get_info(self):
        print(f"Half-openning angle = {3/2* np.sqrt(self.emission_height) / np.sqrt(self.Bmodel.RLC) * 180 / np.pi}", f"x = {self.beta / (3/2* np.sqrt(self.emission_height) / np.sqrt(self.Bmodel.RLC))}")

    def get_PA(self, phase, r_esc):
        if len(np.array(r_esc).shape) == 0:
            ray = Ray(self.Bmodel, self.observer_vec, phase, self.emission_height)
            return ray.beta_eff(r_esc) * 180 /np.pi
        else:
            r_esc_arr = np.array(r_esc)
            ray = Ray(self.Bmodel, self.observer_vec, phase, self.emission_height)
            beta_eff_arr = np.zeros_like(r_esc_arr)
            for ind, r in enumerate(r_esc_arr):
                beta_eff_arr[ind] = ray.beta_eff(r) * 180 /np.pi
            return beta_eff_arr

            
    

