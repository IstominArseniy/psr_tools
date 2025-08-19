from psr_lib.test_package import test_module
# from psr_lib.interface import profile_class
from psr_lib import interface
I = [0, 1, 0]
Q = [0, 0.5, 0]
U = [0, 0.5, 0]
V = [0, 0.1, 0]
Profile = interface.PulsarProfile.from_IQUV(I, Q, U, V)
print(Profile.get_W10())
FAST_Profile = interface.read_FAST_data('Data/FAST_Profiles/J0139p3336ILVPA210110.dat')
print(FAST_Profile.get_W50())
MeerKAT_Profile = interface.read_MeerKAT_mean_data('Data/MeerKAT_Profiles/ar_files/J0555-7056_2020-06-26-05:32:56_zap.8chTS.fluxcal.ar')
print(MeerKAT_Profile.get_W50())
MeerKAT_channel_Profile = interface.read_MeerKAT_channel_data('Data/MeerKAT_Profiles/ar_files/J0555-7056_2020-06-26-05:32:56_zap.8chTS.fluxcal.ar', 2)[0]
print(MeerKAT_channel_Profile.get_W50())



