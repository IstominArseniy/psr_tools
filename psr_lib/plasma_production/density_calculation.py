from ..interface import RadioPulsar
from tqdm.autonotebook import tqdm
import numpy as np
import warnings
from ..interface import DensityProfile1D
# from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map
from functools import partial

def _find_multiplicity_1D(x, psr:RadioPulsar, gamma_func:callable):
    return psr.curvature_multipliticy(x, gamma_func(x))

def find_density_profile_1D(psr:RadioPulsar, gamma_func:callable, Npoints=10, max_workers=12) -> DensityProfile1D:
    #TODO multiprocessing, better gamma_arr
    x_min = 0.05
    n_arr = np.zeros(Npoints)
    x_arr = np.linspace(x_min, 1, Npoints)
    with warnings.catch_warnings(record=True, ) as warning_list:
        worker = partial(_find_multiplicity_1D, psr=psr, gamma_func=gamma_func)
        n_arr = process_map(worker, x_arr, max_workers=max_workers)
    if len(warning_list) != 0:
        print(f"{len(warning_list)} warnings were catched!")
    return DensityProfile1D(x_arr, n_arr)


def find_density_profile_2D(psr):
    pass