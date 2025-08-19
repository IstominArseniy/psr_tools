import numpy as np
from matplotlib import pyplot as plt
import scipy.interpolate

from psr_lib import data_processing as processing



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
    def from_IQUV(cls, I, Q, U, V):
        """
        initialize profile class using Stocks parameters in order: I, Q, U, V
        """
        I = np.array(I)
        Q = np.array(Q)
        U = np.array(U)
        V = np.array(V)
        if (I.shape[0] != Q.shape[0] or I.shape[0] != U.shape[0] or I.shape[0] != V.shape[0]):
            raise ValueError("I, Q, U, V arrays must have the same length.")
        Ncounts = I.shape[0]
        profile = cls(Ncounts)
        profile.I = I
        profile.Q = Q
        profile.U = U
        profile.V = V
        profile.L = np.sqrt(Q**2 + U**2)
        profile.L -= processing.noise_mean(profile.L)
        profile.PA = processing.shift_angle(0.5 * np.arctan2(Q, U)) * 180 /np.pi
        return profile

    @classmethod
    def from_ILVPA(cls, I, L, V, PA):
        """
        initialize profile class using I, L, V, PA
        """
        I = np.array(I)
        L = np.array(L)
        V = np.array(V)
        PA = np.array(PA)
        if (I.shape[0] != L.shape[0] or I.shape[0] != V.shape[0] or I.shape[0] != PA.shape[0]):
            raise ValueError("I, L, V, PA arrays must have the same length.")
        Ncounts = I.shape[0]
        profile = cls(Ncounts)
        profile.I = I
        profile.L = L
        profile.V = V
        profile.PA = PA
        profile.Q = profile.I * np.cos(profile.PA)
        profile.V = profile.I * np.sin(profile.PA)
        return profile
    
    def get_Wa(self, a):
        """
        params: a - level (from 0 to 100 %)
        returns: Wa - width on level a in degrees
        """
        height_a = np.max(self.I) * (a/100)
        # oversampling--------------------------------------------------
        phase = np.linspace(0, 1, self.Ncounts)
        phase_x100 = np.linspace(0, 1, 100 * self.Ncounts)
        Is = scipy.interpolate.interp1d(phase, self.I)(phase_x100) # TODO test smoothed profile (as it was before)
        #---------------------------------------------------------------
        try:
            left_ind = np.where(np.isclose(Is, height_a, rtol=5e-2))[0][0] # leftmost point on level a
            right_ind = np.where(np.isclose(Is, height_a, rtol=5e-2))[0][-1] #rightmost point on level a
        except:
            return np.nan
        return 360 * (right_ind - left_ind) / 100 / self.Ncounts
    
    def get_W10(self):
        return self.get_Wa(10)
    
    def get_W50(self):
        return self.get_Wa(50)
    


        


    




# class PulsarProfile:
#     def load_profile_data(self):
#         """
#         loads profile data
#         data[0] = I, data[1] = Q, data[2] = U, data[3] = V, data[4] = L, data[5] = PA
#         returns: data 
#         """
#         if self.source == 'FAST':
#             data = DataLoaders.load_FAST_data(self.file_name)
#         elif self.source == 'MeerKAT':
#             data = DataLoaders.load_MeerKAT_mean_data(self.file_name)
#         elif self.source == 'EPN':
#             data = DataLoaders.load_EPN_data(self.file_name)
#         elif self.source == 'MeerKAT_one_channel':
#             data, freq = DataLoaders.load_MeerKAT_channel_data(self.file_name, self.channel)
#             self.freq = freq
#         else:
#             print('Unknown source!!!') #TODO raise error
#         return data

#     def __init__(self, file_name, source, channel=0):
#         self.file_name = file_name
#         self.source = source
#         self.channel=channel
#         self.freq = 0 # Hz
#         data = self.load_profile_data()
#         self.I = data[0]
#         self.Q = data[1]
#         self.U = data[2]
#         self.V = data[3]
#         self.L = data[4]
#         self.PA = data[5]
#         self.Ncounts = self.I.shape[0]
#         self.peaks_arr = self.find_peaks()
#         self.profile_type = self.find_profile_type()

#     def get_smoothed_profile(self):
#         noise = AuxFunctions.noise_estimation(self.I)
#         phase = np.linspace(0, 1, self.Ncounts)
#         spl = scipy.interpolate.splrep(phase, self.I, s=1.4 * self.Ncounts*noise**2)
#         I_func = scipy.interpolate.BSpline(*spl)
#         return I_func(phase)

#     def plot_profile(self, plot_pol=True, plot_fit=False, zoom=False):
#         fig, axs = plt.subplots(2, height_ratios=[1, 4])
#         fig.suptitle(self.file_name)
#         noise = AuxFunctions.noise_estimation(self.I)
#         left_ind = 0
#         right_ind = self.Ncounts - 1
#         if zoom == True:
#             left_ind, right_ind = self.get_level_bounds(10)
#             left_ind, right_ind = \
#                 int(max(left_ind - 0.25 * (right_ind - left_ind), 0)), int(min(right_ind + 0.25 * (right_ind - left_ind), self.Ncounts - 1))
#         phase = np.linspace(left_ind / self.Ncounts, (right_ind + 1)/self.Ncounts, right_ind - left_ind + 1)

#         good_L_array = ((self.L[left_ind : right_ind + 1] / (np.abs(self.I[left_ind : right_ind + 1]) + 0.01)) > 0.1) & (self.I[left_ind : right_ind + 1] > 4 * noise)
#         axs[0].set_xlim(left_ind / self.Ncounts, (right_ind + 1)/self.Ncounts)
#         axs[0].scatter(phase[good_L_array], self.PA[left_ind : right_ind + 1][good_L_array], c='black', s=3)
#         axs[1].plot(phase, self.I[left_ind : right_ind + 1], c='black', label='I', linewidth=1)
#         if plot_pol == True:
#             axs[1].plot(phase, self.V[left_ind : right_ind + 1], c='blue', label='V')
#             axs[1].plot(phase, self.L[left_ind : right_ind + 1], c='red', label='L')
#         if plot_fit == True:
#             axs[1].plot(phase, self.get_smoothed_profile()[left_ind : right_ind + 1], c='yellow')
#         fig.legend()
#         fig.show()
#         return fig, axs
        
#     def find_peaks(self):
#         noise = AuxFunctions.noise_estimation(self.I)
#         smoothed_profile = self.get_smoothed_profile()
#         # peaks, info = scipy.signal.find_peaks(smoothed_profile, prominence=np.max(smoothed_profile * 0.05), height=2 * noise) #before here was max(smoothed_profile * 0.05) and height=4*noise
#         peaks, info = scipy.signal.find_peaks(smoothed_profile, prominence = 4 * noise, height=max(4 * noise, np.max(smoothed_profile) * 0.01))
#         return peaks

#     def find_profile_type(self):
#         if self.peaks_arr.shape[0] == 1:
#             return 'single'
#         elif self.peaks_arr.shape[0] >= 2:
#             dists = scipy.spatial.distance.cdist(self.peaks_arr.reshape((-1, 1)), self.peaks_arr.reshape((-1, 1)), lambda u, v: min((u-v)%self.Ncounts, (v-u)%self.Ncounts))
#             maxdist = np.max(dists)
#             if maxdist > 0.9 * self.Ncounts / 2:
#                 return 'orthogonal'
#             else:
#                 if self.peaks_arr.shape[0] == 2:
#                     return 'double'
#                 elif self.peaks_arr.shape[0] == 3:
#                     return 'tripple'
#                 else:
#                     return 'complex'
#         else:
#             return 'error'
        
#     def get_level_bounds(self, level):
#         """
#         params: level - from 0 to 100
#         returns: leftmost and right most indicies correspoding to level
#         """
#         height = np.max(self.I) * (level/100)
#         # oversampling-------
#         phase = np.linspace(0, 1, self.Ncounts)
#         phase_x10 = np.linspace(0, 1, 10 * self.Ncounts)
#         Is = scipy.interpolate.interp1d(phase, self.get_smoothed_profile())(phase_x10)
#         #--------------------
#         inds = np.where(np.isclose(Is, height, rtol=5e-2))[0]
#         if inds.shape[0] == 0:
#             left_ind = 0
#             right_ind = self.Ncounts - 1
#         else:
#             left_ind = inds[0] // 10
#             right_ind = inds[-1] // 10
#         if right_ind - left_ind < 3:
#             left_ind = max(0, left_ind - 3)
#             right_ind = min(self.Ncounts-1, right_ind + 3)

#         return (left_ind, right_ind)

#     def get_Wa(self, a):
#         """
#         params: a - level (from 0 to 100)
#         returns: W_a - width on level a
#         """
#         height_a = np.max(self.I) * (a/100)
#         # oversampling-------
#         phase = np.linspace(0, 1, self.Ncounts)
#         phase_x10 = np.linspace(0, 1, 10 * self.Ncounts)
#         Is = scipy.interpolate.interp1d(phase, self.get_smoothed_profile())(phase_x10)
#         #--------------------
#         try:
#             left_ind = np.where(np.isclose(Is, height_a, rtol=1e-1))[0][0]
#             right_ind = np.where(np.isclose(Is, height_a, rtol=1e-1))[0][-1]
#         except:
#             return np.nan
#         # print(left_ind/self.Ncounts, right_ind/self.Ncounts, self.I[left_ind] / np.max(self.I), self.I[right_ind] / np.max(self.I))
#         return 360 * (right_ind - left_ind) / 10 / self.Ncounts

#         # iil = np.where(np.isclose(pulsar_data[0][:peaks[0]], min(pulsar_data[0][:peaks[0]], key=lambda x:abs(x-height50))))[0][0]
#         # iir = np.where(np.isclose(pulsar_data[0][peaks[0]:], min(pulsar_data[0][peaks[0]:], key=lambda x:abs(x-height50))))[0][0] + peaks[0]
#         # if iir <= peaks[1]:
#         #     height50 = np.max(pulsar_data[0][peaks[1]])/2
#         #     iir = np.where(np.isclose(pulsar_data[0][peaks[1]:], min(pulsar_data[0][peaks[1]:], key=lambda x:abs(x-height50))))[0][0] + peaks[1]

#     def get_mode(self):
#         left_ind, right_ind = self.get_level_bounds(10)
#         # print(left_ind, right_ind)
#         Varray = self.V[left_ind:right_ind]
#         PAarray = self.PA[left_ind:right_ind]
#         noise = AuxFunctions.noise_estimation(self.I)
#         good_L_array = ((self.L[left_ind:right_ind] / (np.abs(self.I[left_ind:right_ind]) + 0.01)) > 0.1) & (self.I[left_ind:right_ind] > 4 * noise)
#         N = Varray.shape[0]
#         xs = np.arange(0, N, 1)
#         if len(xs[good_L_array]) <= 3:
#             return '?'
#         spl = scipy.interpolate.splrep(xs[good_L_array], PAarray[good_L_array], s=N*5**2)
#         PA_func = scipy.interpolate.BSpline(*spl)
#         ders = PA_func(xs, 1) 
#         count = 0
#         for i in range(N):
#             count += ((Varray[i]>0) * 2 - 1) * ((ders[i]>0) * 2 - 1)
#         count = count / N
#         # print(count)
#         if count > 0.4:
#             return 'X'
#         elif count < -0.4:
#             return 'O'
#         else:
#             return '?'