import numpy as np
import scipy.interpolate
from matplotlib import pyplot as plt

class DensityProfile1D:
    def __init__(self, x_arr, n_arr):
        x_arr = np.array(x_arr)
        n_arr = np.array(n_arr)
        self.Npoints = n_arr.shape[0]
        self.x_arr = x_arr
        self.n_arr = n_arr
        self.interpolant = self._create_interpolation()
        self.multiplicity = self._find_multiplicity()

    def show(self):
        fig, ax = plt.subplots()
        xs = np.linspace(self.x_arr[0], self.x_arr[-1], 100)
        ax.plot(xs, self.get_n(xs))
        ax.scatter(self.x_arr, self.n_arr)
        fig.show()
        return fig, ax

    def _create_interpolation(self):
        return scipy.interpolate.make_interp_spline(self.x_arr, self.n_arr, k=1)

    def get_n(self, x, normalized=False):
        x = np.array(x)
        if np.any(x<0) or np.any(x>1):
            raise ValueError("x must be from 0 to 1.")
        if not normalized:
            return self.interpolant(x)
        else:
            return self.interpolant(x) / self.multiplicity
    
    def write_to_file(self, filename):
        pass

    def _find_multiplicity(self):
        scipy.integrate.trapezoid(self.n_arr, self.x_arr)

class DenistyProfile2D:
    def __init__(self, n_arr):
        pass

    def show(self):
        pass

    def get_1D_slice(self, phi):
        pass

    def _create_interpolation(self):
        pass
