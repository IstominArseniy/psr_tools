import numpy as np
import scipy.stats


def noise_estimation(data):
    N = data.shape[0]
    n8 = N // 8
    sigmas = []
    for i in range(8):
        sigmas.append(np.var(data[i * n8: (i+1) * n8]))
    return np.sqrt(np.min(sigmas))
    
def noise_mean(data):
    N = data.shape[0]
    n8 = N // 8
    std = np.max(data) - np.min(data)
    mean = 0
    for i in range(8): # find mean of the minimum std piece
        if np.sqrt(np.var(data[i * n8: (i+1) * n8])) < std:
            std = np.sqrt(np.var(data[i * n8: (i+1) * n8]))
            mean = np.mean(data[i * n8: (i+1) * n8])
    return mean



def calculate_width(chi, beta, rho):
    dzeta = beta + chi
    if np.sqrt(np.sin((rho + beta)/2) * np.sin((rho - beta)/2) / np.sin(chi) / np.sin(dzeta)) >= 1:
        return 2 * np.pi
    return 4 * np.arcsin(np.sqrt(np.sin((rho + beta)/2) * np.sin((rho - beta)/2) \
    / np.sin(chi) / np.sin(dzeta)))

def find_power_law_ind(xs, ys):
    log_xs = np.log(xs)
    log_ys = np.log(ys)
    res = scipy.stats.linregress(log_xs, log_ys)
    return res.slope, res.rvalue

def shift_angle(angle):
    angle = np.asarray(angle)
    scalar_in = (angle.ndim==0)
    angle_vals = np.array(angle, copy=False, ndmin=1, dtype=float)
    for i in range(angle_vals.shape[0]):
        while angle_vals[i] > np.pi / 2:
            angle_vals[i] -= np.pi
        while angle_vals[i] < -np.pi / 2:
            angle_vals[i] += np.pi
    if scalar_in:
        return angle_vals[0]
    else:
        return angle_vals