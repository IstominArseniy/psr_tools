import numpy as np
from matplotlib import pyplot as plt
import scipy.interpolate
import scipy.signal
import scipy.interpolate

from psr_tools import data_processing as processing



class PulsarProfile:
    def __init__(self, Ncounts):
        """
        initialize blank profile
        """
        self.Ncounts = Ncounts # number of points
        self.I = np.zeros(Ncounts)
        self.Q = np.zeros(Ncounts)
        self.U = np.zeros(Ncounts)
        self.V = np.zeros(Ncounts)
        self.L = np.zeros(Ncounts)
        self.PA = np.zeros(Ncounts)

    @classmethod
    def from_IQUV(cls, I, Q, U, V, d_phi=None, normalize=True):
        """
        initialize profile class using Stocks parameters in order: I, Q, U, V
        d_phi can be provided if phases do not span 360 degrees. In this case profile I, Q, U, V are padded with zeros.
        """
        I = np.array(I)
        Q = np.array(Q)
        U = np.array(U)
        V = np.array(V)
        if d_phi is not None:
            I = PulsarProfile._pad_arr_with_zeros(I, d_phi)
            Q = PulsarProfile._pad_arr_with_zeros(Q, d_phi)
            U = PulsarProfile._pad_arr_with_zeros(U, d_phi)
            V = PulsarProfile._pad_arr_with_zeros(V, d_phi)

        if (I.shape[0] != Q.shape[0] or I.shape[0] != U.shape[0] or I.shape[0] != V.shape[0]):
            raise ValueError("I, Q, U, V arrays must have the same length.")
        Ncounts = I.shape[0]
        profile = cls(Ncounts)
        profile.I = I
        profile.Q = Q
        profile.U = U
        profile.V = V
        if normalize:
            Imax = np.max(profile.I)
            profile.I /= Imax
            profile.Q /= Imax
            profile.U /= Imax
            profile.V /= Imax
        profile.L = np.sqrt(Q**2 + U**2)
        profile.L -= processing.noise_mean(profile.L)
        profile.PA = processing.shift_angle(0.5 * np.arctan2(Q, U)) * 180 /np.pi
        return profile

    @classmethod
    def from_ILVPA(cls, I, L, V, PA, d_phi=None, normalize=True):
        """
        initialize profile class using I, L, V, PA
        d_phi can be provided if phases do not span 360 degrees. In this case profile I, L, V, are padded with zeros and PA with Nans
        """
        I = np.array(I)
        L = np.array(L)
        V = np.array(V)
        PA = np.array(PA)
        if d_phi is not None:
            I = PulsarProfile._pad_arr_with_zeros(I, d_phi)
            L = PulsarProfile._pad_arr_with_zeros(L, d_phi)
            V = PulsarProfile._pad_arr_with_zeros(V, d_phi)
            PA = PulsarProfile._pad_arr_with_zeros(PA, d_phi)
        if (I.shape[0] != L.shape[0] or I.shape[0] != V.shape[0] or I.shape[0] != PA.shape[0]):
            raise ValueError("I, L, V, PA arrays must have the same length.")
        Ncounts = I.shape[0]
        profile = cls(Ncounts)
        profile.I = I
        profile.L = L
        profile.V = V
        profile.PA = PA
        if normalize:
            Imax = np.max(profile.I)
            profile.I /= Imax
            profile.L /= Imax
            profile.V /= Imax
        profile.Q = profile.I * np.cos(profile.PA * np.pi / 180)
        profile.U = profile.I * np.sin(profile.PA * np.pi / 180)
        return profile
    
    @classmethod
    def from_LOKI_dat(cls, filename, normalize=True):
        """
        tmp method to read profile from .dat files of LOKI-pulsar-propagation program
        """
        phis = []
        Is = []
        Ls = []
        Vs = []
        PAs = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                phi, I, V, PA = map(float, line.split(' '))
                phis.append(phi)
                Is.append(I)
                Vs.append(V)
                PAs.append(PA)    
                Ls.append(np.sqrt(I**2 - V**2))
        idx_sort = np.argsort(phis)
        phis = np.array(phis)[idx_sort]
        Is = np.array(Is)[idx_sort]
        Ls = np.array(Ls)[idx_sort]
        Vs = np.array(Vs)[idx_sort]
        PAs = np.array(PAs)[idx_sort]
        
        #------------- padding with zeros on sides ----------
        d_phi = phis[1] - phis[0]
        # Ncounts_target = np.round(360 / d_phi) + 1
        # Ncounts_current = Is.shape[0]
        # pad = np.zeros(int((Ncounts_target - Ncounts_current)/2))
        # Is = np.concatenate((pad, Is, pad))
        # Ls = np.concatenate((pad, Ls, pad))
        # Vs = np.concatenate((pad, Vs, pad))
        # PAs = np.concatenate((pad, PAs, pad))
        Is = PulsarProfile._pad_arr_with_zeros(Is, d_phi)
        Ls = PulsarProfile._pad_arr_with_zeros(Ls, d_phi)
        Vs = PulsarProfile._pad_arr_with_zeros(Vs, d_phi)
        PAs = PulsarProfile._pad_arr_with_zeros(PAs, d_phi)
        #----------------------------------------------------
        if (Is.shape[0] != Ls.shape[0] or Is.shape[0] != Vs.shape[0] or Is.shape[0] != PAs.shape[0]):
            raise ValueError("I, L, V, PA arrays must have the same length.")
        Ncounts = Is.shape[0]
        profile = cls(Ncounts)
        profile.I = Is
        profile.L = Ls
        profile.V = Vs
        profile.PA = PAs
        if normalize:
            Imax = np.max(profile.I)
            profile.I /= Imax
            profile.L /= Imax
            profile.V /= Imax
        profile.Q = profile.I * np.cos(profile.PA * np.pi / 180)
        profile.U = profile.I * np.sin(profile.PA * np.pi / 180)
        return profile
    
    @staticmethod
    def _pad_arr_with_zeros(arr, d_phi):
        Ncounts_target = np.round(360 / d_phi) + 1
        Ncounts_current = arr.shape[0]
        pad = np.zeros(int((Ncounts_target - Ncounts_current)/2))
        arr = np.concatenate((pad, arr, pad))
        return arr
    
    @staticmethod
    def _pad_arr_with_nans(arr, d_phi):
        Ncounts_target = np.round(360 / d_phi) + 1
        Ncounts_current = arr.shape[0]
        pad = np.zeros(int((Ncounts_target - Ncounts_current)/2), dtype=np.float32)
        pad.fill(np.nan)
        arr = np.concatenate((pad, arr, pad))
        return arr
    
    def get_Wa(self, a):
        """
        params: a - level (from 0 to 100 %)
        returns: Wa - width on level a in degrees
        """
        height_a = np.max(self.I) * (a/100)
        # oversampling--------------------------------------------------
        phase = np.linspace(0, 1, self.Ncounts)
        phase_x100 = np.linspace(0, 1, 100 * self.Ncounts)
        Is = scipy.interpolate.interp1d(phase, self.get_smoothed_profile())(phase_x100) # TODO test smoothed profile (as it was before)
        #---------------------------------------------------------------
        try:
            left_ind = np.where(np.isclose(Is, height_a, rtol=5e-2, atol=0))[0][0] # leftmost point on level a
            right_ind = np.where(np.isclose(Is, height_a, rtol=5e-2, atol=0))[0][-1] #rightmost point on level a
        except:
            return np.nan
        return 360 * (right_ind - left_ind) / 100 / self.Ncounts
    
    def get_W10(self):
        return self.get_Wa(10)
    
    def get_W50(self):
        return self.get_Wa(50)
    
    def get_level_bounds(self, level):
        """
        params: level - from 0 to 100 %
        returns: leftmost and rightmost indicies correspoding to level
        """
        height = np.max(self.I) * (level/100)
        # oversampling------------------------------------------
        phase = np.linspace(0, 1, self.Ncounts)
        phase_x100 = np.linspace(0, 1, 100 * self.Ncounts)
        Is = scipy.interpolate.interp1d(phase, self.I)(phase_x100) # TODO test smoothed profile (as it was before)
        #--------------------------------------------------------
        inds = np.where(np.isclose(Is, height, rtol=5e-2))[0]

        if inds.shape[0] == 0:
            left_ind = 0
            right_ind = self.Ncounts - 1
        else:
            left_ind = inds[0] // 100
            right_ind = inds[-1] // 100

        if right_ind - left_ind < 3: # provide at least three points interval # TODO make more accurate bounds in this case
            left_ind = max(0, left_ind - 3)
            right_ind = min(self.Ncounts-1, right_ind + 3)

        return (left_ind, right_ind)

    def find_emission_mode(self, boarder_value=0.4):
        left_ind, right_ind = self.get_level_bounds(10)
        Vs = self.V[left_ind:right_ind]
        PAs = processing.shift_angle(self.PA[left_ind:right_ind]*np.pi/180)*180/np.pi
        Imax = np.max(self.I)
        Is = self.I[left_ind:right_ind] / Imax
        Ls = self.L[left_ind:right_ind] / Imax
        noise = processing.noise_estimation(self.I) / Imax

        quality_L_mask = ((Ls / (np.abs(Is) + 0.01)) > 0.1) & (Is > 4 * noise)
        N = Vs.shape[0]
        xs = np.arange(0, N, 1) # arbitrary "even steps" array
        if len(xs[quality_L_mask]) <= 3:
            return '?' # Not enough quality points 
        # PA qubic spline interpolation ------------------------------------------------
        spl = scipy.interpolate.splrep(xs[quality_L_mask], PAs[quality_L_mask], s=N*5**2)
        PA_func = scipy.interpolate.BSpline(*spl)
        # ------------------------------------------------------------------------------
        ders = PA_func(xs, 1) 
        count = 0
        for i in range(N):
            count += ((Vs[i]>0) * 2 - 1) * ((ders[i]>0) * 2 - 1)
        if count / N > boarder_value:
            return 'X'
        elif count / N < -boarder_value:
            return 'O'
        else:
            return '?'
        
    def get_smoothed_profile(self): # TODO Write general routine to smooth arrays based on noise estimation
        """
        smooth intencity profile using smooth qubic spilnes
        """
        noise = processing.noise_estimation(self.I)
        phase = np.linspace(0, 1, self.Ncounts)
        spl = scipy.interpolate.splrep(phase, self.I, s=1.4 * self.Ncounts*noise**2)
        I_func = scipy.interpolate.BSpline(*spl)
        return I_func(phase)

    def find_peaks(self):
        """
        find peaks in intensity profile
        """
        noise = processing.noise_estimation(self.I)
        smoothed_profile = self.get_smoothed_profile()
        peaks, info = scipy.signal.find_peaks(smoothed_profile, prominence = 4 * noise, height=max(4 * noise, np.max(smoothed_profile) * 0.01))
        return peaks

    def find_profile_type(self):
        peaks_arr = self.find_peaks()
        if peaks_arr.shape[0] == 1:
            return 'single'
        elif peaks_arr.shape[0] >= 2:
            dists = scipy.spatial.distance.cdist(peaks_arr.reshape((-1, 1)), peaks_arr.reshape((-1, 1)), lambda u, v: min((u-v)%self.Ncounts, (v-u)%self.Ncounts))
            maxdist = np.max(dists)
            if maxdist > 0.9 * self.Ncounts / 2:
                return 'orthogonal'
            else:
                if peaks_arr.shape[0] == 2:
                    return 'double'
                elif peaks_arr.shape[0] == 3:
                    return 'tripple'
                else:
                    return 'complex'
        else:
            print("Error occured during profile type identifiction. Nan has been returned.")
            return np.nan
    
    def plot_profile(self, plot_polarisation=True, plot_fit=False, zoom=False):
        fig, axs = plt.subplots(2, height_ratios=[1, 4])
        noise = processing.noise_estimation(self.I)
        left_ind = 0 # leftmost index to plot
        right_ind = self.Ncounts - 1 # rightmost index to plot
        if zoom == True: # find W10 and add 0.25 of W10 length to left and to right
            left_ind, right_ind = self.get_level_bounds(10)
            left_ind, right_ind = \
                int(max(left_ind - 0.25 * (right_ind - left_ind), 0)), int(min(right_ind + 0.25 * (right_ind - left_ind), self.Ncounts - 1))
            
        phase_arr = np.linspace(left_ind / self.Ncounts, (right_ind + 1)/self.Ncounts, right_ind - left_ind + 1)
        Is = self.I[left_ind : right_ind + 1]
        Ls = self.L[left_ind : right_ind + 1]
        Vs = self.V[left_ind : right_ind + 1]
        PAs = self.PA[left_ind : right_ind + 1]
        quality_L_array = ((Ls / (np.abs(Is) + 0.01 * np.max(Is))) > 0.1) & (Is > 4 * noise)
        axs[0].set_xlim(left_ind / self.Ncounts, (right_ind + 1)/self.Ncounts)
        axs[0].scatter(phase_arr[quality_L_array], processing.shift_angle(PAs[quality_L_array] * np.pi/180)*180/np.pi, c='black', s=3)
        axs[1].plot(phase_arr, Is, c='black', label='I', linewidth=1)
        if plot_polarisation == True:
            axs[1].plot(phase_arr, Vs, c='blue', label='V')
            axs[1].plot(phase_arr, Ls, c='red', label='L')
        if plot_fit == True:
            axs[1].plot(phase_arr, self.get_smoothed_profile()[left_ind : right_ind + 1], c='yellow')
        fig.legend()
        fig.show()
        return fig, axs
        

