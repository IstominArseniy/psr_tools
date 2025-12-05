import numpy as np
from astropy.io import fits
from psr_tools import data_processing as processing
from psr_tools import interface

from matplotlib import pyplot as plt

def read_FAST_data(file_name):
    Is = np.array(0)
    Ls = np.array(0)
    Vs = np.array(0)
    PAs = np.array(0)
    with open(file_name, 'r') as file:
        lines = file.readlines()
        Is = np.zeros(len(lines) - 1)
        Ls = np.zeros(len(lines) - 1)
        Vs = np.zeros(len(lines) - 1)
        PAs = np.zeros(len(lines) - 1)
        for i in range(1, len(lines)):
            ILVPA = list(map(float, lines[i].split(' ')))
            Is[i-1] = ILVPA[1]
            Ls[i-1] = ILVPA[2]
            Vs[i-1] = ILVPA[3]
            PAs[i-1] = ILVPA[4]
    return interface.PulsarProfile.from_ILVPA(Is, Ls, Vs, PAs)

def read_MeerKAT_mean_data(file_name):
    hdul = fits.open(file_name)
    DM = hdul[4].header['DM']
    RM = hdul[4].header['RM']
    DC = DM / 2.41e-16
    Jname = hdul[0].header['SRC_NAME']
    P = hdul[4].header['PERIOD']
    freqs = hdul[4].data[0]['DAT_FREQ'] * 1e6
    N = hdul[4].data[0]['DATA'][0][0].shape[0]
    Is = np.zeros(N)
    Qs = np.zeros(N)
    Us = np.zeros(N)
    Vs = np.zeros(N)
    add_shift = (np.argmax(hdul[4].data[0]['DATA'][0][0]) - N//4) 
    for channel in range(8):
        for stocks in range(4):
            hdul[4].data[0]['DATA'][stocks][channel] = np.roll(hdul[4].data[0]['DATA'][stocks][channel], int(DC * (1/freqs[0]**2 - 1/freqs[channel]**2) * N / P) - add_shift)
    for channel in range(8):
        weight = hdul[4].data[0]['DAT_WTS'][channel] / np.sum(hdul[4].data[0]['DAT_WTS'])
        #---------I, Q, U, V--------------
        Is += (hdul[4].data[0]['DATA'][0][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][0][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*0+channel] * weight
        Q = (hdul[4].data[0]['DATA'][1][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][1][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*1+channel] 
        U = (hdul[4].data[0]['DATA'][2][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][2][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*2+channel]
        beta = RM * ((299792458 / freqs[0])**2 - (299792458 / freqs[channel])**2) 
        Qs += weight * (np.cos(2 * beta) * Q - np.sin(2 * beta) * U)
        Us += weight * (np.sin(2 * beta) * Q + np.cos(2 * beta) * U)
        Vs += (hdul[4].data[0]['DATA'][3][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][3][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*3+channel] * weight
    return interface.PulsarProfile.from_IQUV(Is, Qs, Us, Vs)
    

def read_MeerKAT_channel_data(file_name, channel):
    hdul = fits.open(file_name)
    DM = hdul[4].header['DM']
    DC = DM / 2.41e-16
    Jname = hdul[0].header['SRC_NAME']
    P = hdul[4].header['PERIOD']
    freqs = hdul[4].data[0]['DAT_FREQ'] * 1e6
    N = hdul[4].data[0]['DATA'][0][0].shape[0]
    Is = np.zeros(N)
    Qs = np.zeros(N)
    Us = np.zeros(N)
    Vs = np.zeros(N)
    add_shift = (np.argmax(hdul[4].data[0]['DATA'][0][0]) - N//4) 
    # TODO remove RM
    for stocks in range(4):
        hdul[4].data[0]['DATA'][stocks][channel] = np.roll(hdul[4].data[0]['DATA'][stocks][channel], int(DC * (1/freqs[0]**2 - 1/freqs[channel]**2) * N / P) - add_shift)
    Is = (hdul[4].data[0]['DATA'][0][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][0][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*0+channel] * (hdul[4].data[0]['DAT_WTS'][channel] / np.sum(hdul[4].data[0]['DAT_WTS']))
    Qs = (hdul[4].data[0]['DATA'][1][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][1][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*1+channel] * (hdul[4].data[0]['DAT_WTS'][channel] / np.sum(hdul[4].data[0]['DAT_WTS']))
    Us = (hdul[4].data[0]['DATA'][2][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][2][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*2+channel] * (hdul[4].data[0]['DAT_WTS'][channel] / np.sum(hdul[4].data[0]['DAT_WTS']))
    Vs = (hdul[4].data[0]['DATA'][3][channel].astype('int64') - processing.noise_mean(hdul[4].data[0]['DATA'][3][channel].astype('int64'))) * hdul[4].data[0]['DAT_SCL'][8*3+channel] * (hdul[4].data[0]['DAT_WTS'][channel] / np.sum(hdul[4].data[0]['DAT_WTS']))
    return interface.PulsarProfile.from_IQUV(Is, Qs, Us, Vs), freqs[channel]

def load_EPN_data(file_name):
    pass


