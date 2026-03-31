import numpy as np
import scipy.interpolate
from matplotlib import pyplot as plt
import json
from ..utils import plotters

class DensityProfile1D:
    def __init__(self, x_arr, n_arr):
        x_arr = np.array(x_arr)
        n_arr = np.array(n_arr)
        self.Npoints = n_arr.shape[0]
        self.x_arr = x_arr
        self.n_arr = n_arr
        self.interpolant = self._create_interpolation()
        self.multiplicity = self._find_multiplicity()

    def show(self, file_name=None):
        fig, ax = plt.subplots()
        xs = np.linspace(self.x_arr[0], self.x_arr[-1], 100)
        ax.plot(xs, self.get_n(xs))
        ax.scatter(self.x_arr, self.n_arr)
        if file_name is not None:
            fig.savefig(file_name, dpi=400, bbox_inches='tight')        
        fig.show()
        return fig, ax

    def _create_interpolation(self):
        return scipy.interpolate.make_interp_spline(self.x_arr, self.n_arr, k=1)
        # return scipy.interpolate.RegularGridInterpolator((self.x_arr,), self.n_arr, method='linear', bounds_error=False, fill_value=0)
    

    def get_n(self, x, normalized=False):
        x = np.array(x)
        if np.any(x<0) or np.any(x>1):
            raise ValueError("x must be from 0 to 1.")
        if not normalized:
            return self.interpolant(x)
        else:
            return self.interpolant(x) / self.multiplicity
    
    def write_to_file(self, filename, normalized=False):
        tmp_dict = dict()
        tmp_dict['x_arr'] = self.x_arr.tolist()
        tmp_dict['n_arr'] = self.n_arr.tolist()
        with open (filename + '.json', 'w') as f:
            json.dump(tmp_dict, f)

    def _find_multiplicity(self):
        return scipy.integrate.trapezoid(self.n_arr, self.x_arr)

class DensityProfile2D:
    def __init__(self, x_arr, phi_arr, n_arr):
        x_arr = np.array(x_arr)
        phi_arr = np.array(phi_arr)
        n_arr = np.array(n_arr)
        self.Npoints_x = n_arr.shape[0]
        self.Npoints_phi = n_arr.shape[1]
        self.x_arr = x_arr
        self.phi_arr = phi_arr
        self.n_arr = n_arr
        self.interpolant = self._create_interpolation()
        self.multiplicity = self._find_multiplicity()

    def show(self, file_name=None):
        plotters.polar_plot(self.n_arr, self.x_arr, self.phi_arr, file_name=file_name)

    def get_n(self, x, phi):
        return self.interpolant((x, phi))

    def get_1D_slice(self, phi, Npoints=None):
        if Npoints is None:
            Npoints=self.Npoints_x
        xs = np.linspace(self.x_arr[0], self.x_arr[-1], Npoints)
        ns = self.get_n(xs, phi)
        return DensityProfile1D(xs, ns)

    def write_to_file(self, filename, normalized=False):
        tmp_dict = dict()
        tmp_dict['x_arr'] = self.x_arr.tolist()
        tmp_dict['phi_arr'] = self.phi_arr.tolist()
        tmp_dict['n_arr'] = self.n_arr.tolist()
        with open (filename + '.json', 'w') as f:
            json.dump(tmp_dict, f)

    def _create_interpolation(self): #TODO cyclic interpolation
        return scipy.interpolate.RegularGridInterpolator((self.x_arr, self.phi_arr), self.n_arr, method='linear', bounds_error=False, fill_value=0)
    
    def _find_multiplicity(self):
        return scipy.integrate.trapezoid([scipy.integrate.trapezoid(n_arr_x_cut * self.x_arr, self.x_arr) for n_arr_x_cut in self.n_arr], self.phi_arr)