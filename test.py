import numpy as np
from matplotlib import pyplot as plt

from psr_lib.test_package import test_module
from psr_lib import interface

phase = np.linspace(0, 1, 32)
I = phase * (1 - phase)
Q = I / 3
U = I / 3
V = I / 3
Profile = interface.PulsarProfile.from_IQUV(I, Q, U, V)
print(Profile.get_W10())
FAST_Profile = interface.read_FAST_data('Data/FAST_Profiles/J0139p3336ILVPA210110.dat')
print(FAST_Profile.get_W50())
MeerKAT_Profile = interface.read_MeerKAT_mean_data('Data/MeerKAT_Profiles/ar_files/J0536-7543_2019-10-18-23:49:36_zap.8chTS.fluxcal.ar')
print(MeerKAT_Profile.get_W50())
MeerKAT_channel_Profile = interface.read_MeerKAT_channel_data('Data/MeerKAT_Profiles/ar_files/J0555-7056_2020-06-26-05:32:56_zap.8chTS.fluxcal.ar', 2)[0]
print(MeerKAT_channel_Profile.get_W50())
print(MeerKAT_Profile.find_emission_mode(), MeerKAT_Profile.find_profile_type())
MeerKAT_Profile.plot_profile(plot_polarisation=True, zoom=True)

psr = interface.ObservedRadioPulsar('test_pulsar', 1, 1, 1, 45, 2, 600)
psr.write_to_file('test_pulsar')

new_psr = interface.ObservedRadioPulsar.from_file('test_pulsar.toml')
print(new_psr.freq)