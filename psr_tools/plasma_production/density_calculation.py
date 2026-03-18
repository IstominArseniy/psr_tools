from ..interface import RadioPulsar
from tqdm.autonotebook import tqdm
import numpy as np
import warnings
from ..interface import DensityProfile1D, DensityProfile2D
# from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map
from functools import partial
from itertools import product

def _find_multiplicity_1D(psr:RadioPulsar, x, gamma_func:callable):
    return psr.curvature_multipliticy(x, gamma_func(x))

def _find_multiplicity_2D(psr:RadioPulsar, x_phi_pair, gamma_func:callable):
    return psr.curvature_multipliticy(x_phi_pair[0], gamma_func(x_phi_pair[0], x_phi_pair[1]))

def find_density_profile_1D(psr:RadioPulsar, gamma_func:callable, Npoints=10) -> DensityProfile1D:
    x_min = 0.05 # TODO do smth with x_min
    x_arr = np.linspace(x_min, 1-x_min, Npoints)
    n_arr = np.zeros(Npoints)
    with warnings.catch_warnings(record=True) as warning_list:
        # worker = partial(_find_multiplicity_1D, psr=psr, gamma_func=gamma_func)
        # n_arr = process_map(worker, x_arr, max_workers=max_workers)
        for ind, x in enumerate(x_arr):
            n_arr[ind] = _find_multiplicity_1D(psr, x, gamma_func)
    if len(warning_list) != 0:
        print(f"{len(warning_list)} warnings were catched!")
    return DensityProfile1D(x_arr, n_arr)

# REDO 2D

# def find_density_profile_2D(psr:RadioPulsar, gamma_func:callable, Npoints_x=10, Npoints_phi=10, max_workers=12) -> DensityProfile2D:
#     x_min = 0.05
#     x_arr = np.linspace(x_min, 1, Npoints_x)
#     phi_arr = np.linspace(0, 2 * np.pi, Npoints_phi, endpoint=False)
#     n_arr = np.zeros((Npoints_x, Npoints_phi))
#     with warnings.catch_warnings(record=True) as warning_list:
#         worker = partial(_find_multiplicity_2D, psr=psr, gamma_func=gamma_func)
#         n_arr = np.array(process_map(worker, list(product(x_arr, phi_arr)), max_workers=max_workers)).reshape(Npoints_x, Npoints_phi)
#     if len(warning_list) != 0:
#         print(f"{len(warning_list)} warnings were catched!")
#     return DensityProfile2D(x_arr, phi_arr, n_arr)