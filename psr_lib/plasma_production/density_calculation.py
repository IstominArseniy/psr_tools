from ..interface import RadioPulsar
from tqdm.autonotebook import tqdm
import numpy as np
import warnings
from ..interface import DensityProfile1D

def find_1D_density_profile(psr:RadioPulsar, gamma_arr, Npoints=10) -> DensityProfile1D:
    #TODO multiprocessing, better gamma_arr
    x_min = 0.05
    n_arr = np.zeros(Npoints)
    x_arr = np.linspace(x_min, 1, Npoints)
    with warnings.catch_warnings(record=True, ) as warning_list:
        for i in tqdm(range(Npoints)):
            n_arr[i] = psr.curvature_multipliticy(x_arr[i], gamma_arr[i])
    if len(warning_list) != 0:
        print(f"{len(warning_list)} warnings were catched!")
    return DensityProfile1D(x_arr, n_arr)


def find_2D_density_profile(psr):
    pass