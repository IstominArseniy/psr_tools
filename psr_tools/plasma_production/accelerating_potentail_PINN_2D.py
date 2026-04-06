import os
from pathlib import Path
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import scipy.special
import scipy.integrate
import scipy.optimize
import functools
import scipy.interpolate
import warnings
from typing import Callable, Iterable
import torch
import torch.nn as nn
from torch.autograd import Variable
# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
device = "cpu"
import psutil
import torch.nn.functional as F
from ..interface import RadioPulsar
from ..utils import constants


class PsiPINN2D(nn.Module):
    def __init__(self):
        super(PsiPINN2D, self).__init__()
        self.hidden_layer1 = nn.Linear(2,20)
        self.hidden_layer2 = nn.Linear(20,20)
        self.hidden_layer3 = nn.Linear(20,20)
        self.hidden_layer4 = nn.Linear(20,20)
        self.output_layer = nn.Linear(20,1)

    def forward(self, r, z):
        inputs = torch.cat([r,z],axis=1) # combined two arrays of 1 columns each to one array of 2 columns
        layer1_out = torch.tanh(self.hidden_layer1(inputs))
        layer2_out = torch.tanh(self.hidden_layer2(layer1_out))
        layer3_out = torch.tanh(self.hidden_layer3(layer2_out))
        layer4_out = torch.tanh(self.hidden_layer4(layer3_out))
        output = torch.tanh(z) * self.output_layer(layer4_out)
        return output
    
class PoissonSolver:
    def __init__(self, PSR:RadioPulsar, h0:callable, r_eps=0.0, w_bc = 0.8):
        self.PSR = PSR
        self.h0 = h0 # initial gap height profile
        self.net = PsiPINN2D() # neural network class to solve Poisson equation
        self.r_eps = r_eps # shift from 0 and 1 boarders
        self.w_bc = w_bc # boundary condition loss function weight
        self.net = self.net.to(device)
        self.mse_cost_function = torch.nn.MSELoss() # Mean squared error
        self.optimizer = torch.optim.Adam(self.net.parameters(), lr=1e-3) # NN optimizer
        # boarder coordinates
        self.r_bc1 = np.ones((500, 1)) * self.r_eps 
        self.r_bc2 = np.linspace(self.r_eps, 1 - self.r_eps, 500).reshape((500, 1))
        self.r_bc3 = np.ones((500, 1)) * (1 - self.r_eps)
        self.r_bc4 = np.linspace(self.r_eps, 1 - self.r_eps, 500).reshape((500, 1))
        self.z_bc1 = np.linspace(0, self.h0(self.r_eps), 500).reshape((500, 1))
        self.z_bc2 = self.h0(self.r_bc2)
        self.z_bc3 = np.linspace(0, self.h0(1 - self.r_eps), 500).reshape((500, 1))
        self.z_bc4 = np.zeros((500, 1))
        # boarder conditions
        self.dpsi_bc1 = np.zeros((500, 1))
        self.dpsi_bc2 = np.zeros((500, 1))
        self.psi_bc3 = np.zeros((500, 1))
        self.psi_bc4 = np.zeros((500, 1))
        self.losses = [] # total losses list
        self.pde_losses = [] # PDE losses list
        self.bc_losses = [] # boundary conditions losses list

    def PDEres(self, r, z):
        """
        PDE residual as a function of (r, z)
        r - transvers distance from magnetic axis in polar cap radius
        z - height in polar cap radius
        return: pde residual
        """
        psi = self.net(r, z)
        psi_r = torch.autograd.grad(psi.sum(), r, create_graph=True)[0]
        psi_z = torch.autograd.grad(psi.sum(), z, create_graph=True)[0]
        psi_rr = torch.autograd.grad(psi_r.sum(), r, create_graph=True)[0]
        psi_zz = torch.autograd.grad(psi_z.sum(), z, create_graph=True)[0]
        pde = psi_rr + psi_zz + 1 / r * psi_r + 4 * np.cos(self.PSR.chi) * self.PSR.Kpsi
        return pde

    def dpsi_dr(self, r, z):
        """
        derivative of the potential with respect to r
        """
        psi = self.net(r, z)
        psi_r = torch.autograd.grad(psi.sum(), r, create_graph=True)[0]
        return psi_r

    def dpsi_dz(self, r, z):
        """
        derivative of the potential with respect to z
        """
        psi = self.net(r, z)
        psi_z = torch.autograd.grad(psi.sum(), z, create_graph=True)[0]
        return psi_z

    def batch_generator(self, N):
        """
        generates a batch of N randomly uniformly distributed points within the computation domain 
        """
        delta = 0.01
        rs = np.random.uniform(low=delta, high=1.0 - delta, size=(N, 1))
        zs = np.zeros((N, 1))
        for i in range(N):
            zs[i] = np.random.uniform(low=0.0, high=self.h0(rs[i]))
        return rs, zs

    def optimize_psi(self, n_batches, epochs_per_batch, batch_size):
        """
        main function to optimize solution to the Poisson equation
        n_bathces - number of bathces (each bathc correspond to a new set of learngin potins)
        epochs_per_batch - number of epochs in one batch optimization
        batch_size = size of batch
        """
        self.net.train()
        rs = np.zeros((batch_size, 1))
        zs = np.zeros((batch_size, 1))
        all_zeros = np.zeros((batch_size, 1))
        for epoch in range(n_batches * epochs_per_batch): # main optimization loop
            self.optimizer.zero_grad() # to make the gradients zero
            # Boundary condition losses
            pt_r_bc1 = Variable(torch.from_numpy(self.r_bc1).float(), requires_grad=True).to(device)
            pt_r_bc2 = Variable(torch.from_numpy(self.r_bc2).float(), requires_grad=True).to(device)
            pt_r_bc3 = Variable(torch.from_numpy(self.r_bc3).float(), requires_grad=False).to(device)
            pt_r_bc4 = Variable(torch.from_numpy(self.r_bc4).float(), requires_grad=False).to(device)

            pt_z_bc1 = Variable(torch.from_numpy(self.z_bc1).float(), requires_grad=True).to(device)
            pt_z_bc2 = Variable(torch.from_numpy(self.z_bc2).float(), requires_grad=True).to(device)
            pt_z_bc3 = Variable(torch.from_numpy(self.z_bc3).float(), requires_grad=False).to(device)
            pt_z_bc4 = Variable(torch.from_numpy(self.z_bc4).float(), requires_grad=False).to(device)

            pt_dpsi_bc1 = Variable(torch.from_numpy(self.dpsi_bc1).float(), requires_grad=False).to(device)
            pt_dpsi_bc2 = Variable(torch.from_numpy(self.dpsi_bc2).float(), requires_grad=False).to(device)
            pt_psi_bc3 = Variable(torch.from_numpy(self.psi_bc3).float(), requires_grad=False).to(device)
            pt_psi_bc4 = Variable(torch.from_numpy(self.psi_bc4).float(), requires_grad=False).to(device)

            net_bc_out1 = self.dpsi_dr(pt_r_bc1, pt_z_bc1)
            net_bc_out2 = self.dpsi_dz(pt_r_bc2, pt_z_bc2)
            net_bc_out3 = self.net(pt_r_bc3, pt_z_bc3)
            net_bc_out4 = self.net(pt_r_bc4, pt_z_bc4)

            mse_bc_1 = self.mse_cost_function(net_bc_out1, pt_dpsi_bc1)
            mse_bc_2 = self.mse_cost_function(net_bc_out2, pt_dpsi_bc2)
            mse_bc_3 = self.mse_cost_function(net_bc_out3, pt_psi_bc3)
            mse_bc_4 = self.mse_cost_function(net_bc_out4, pt_psi_bc4)


            # PDE losses
            if epoch % epochs_per_batch == 0: # make new batch of points
                rs, zs = self.batch_generator(batch_size)
                print(psutil.virtual_memory().percent)


            pt_rs = Variable(torch.from_numpy(rs).float(), requires_grad=True).to(device) # TODO optimize data transfer to device
            pt_zs = Variable(torch.from_numpy(zs).float(), requires_grad=True).to(device) # TODO optimize data transfer to device
            pt_all_zeros = Variable(torch.from_numpy(all_zeros).float(), requires_grad=False).to(device)

            pde_residual = self.PDEres(pt_rs, pt_zs)
            mse_pde = self.mse_cost_function(pde_residual, pt_all_zeros)

            # Combining the loss functions
            loss = self.w_bc * (mse_bc_1 + mse_bc_2 + mse_bc_3 + mse_bc_4)/4.0 + (1 - self.w_bc) * mse_pde
            with torch.autograd.no_grad(): # to track progress
                if epoch % epochs_per_batch == 0:
                    print(epoch, loss.data)

            self.bc_losses.append(((mse_bc_1 + mse_bc_2 + mse_bc_3 + mse_bc_4)/4.0).data.item())
            self.pde_losses.append(mse_pde.data.item())
            self.losses.append(loss)

            loss.backward() # This is for computing gradients using backward propagation
            self.optimizer.step()

    def get_psi(self, r, z, normalize=True, cut_at_h0=True):
        """
        get potential from the device
        r - array of r coordinates
        z - array of z coordinates
        """
        self.net.eval()
        with torch.autograd.no_grad():
            pt_r = Variable(torch.from_numpy(r).float(), requires_grad=False).to(device)
            pt_z = Variable(torch.from_numpy(z).float(), requires_grad=False).to(device)
            border = Variable(torch.from_numpy(self.h0(r)).float(), requires_grad=False).to(device)
            pt_psi = self.net(pt_r, pt_z)
            pt_psi_border = self.net(pt_r, border)
            psi=pt_psi.data.cpu().numpy()
            psi_border = pt_psi_border.data.cpu().numpy()
            psi = np.where(z < self.h0(r), psi, psi_border)
            return psi

    def get_dpsi_dh(self, r, h, normalize=True, cut_at_h0=True):
        pass

    def set_h0(self, h_func):
        self.h0 = h_func
        self.z_bc1 = np.linspace(0, self.h0(self.r_eps), 500).reshape((500, 1))
        self.z_bc2 = self.h0(self.r_bc2)
        self.z_bc3 = np.linspace(0, self.h0(1 - self.r_eps), 500).reshape((500, 1))


    def show_psi(self, save=False):
        rs = np.arange(0, 1, 0.01)
        zs = np.arange(0, 1, 0.01)
        fig, ax = plt.subplots()
        grid_rs, grid_zs = np.meshgrid(rs, zs)
        rs = np.ravel(grid_rs).reshape(-1,1)
        zs = np.ravel(grid_zs).reshape(-1,1)
        psi = self.get_psi(rs, zs)
        grid_psi = psi.reshape(grid_rs.shape)
        ax.plot(np.arange(0, 1, 0.01), self.h0(np.arange(0, 1, 0.01)), c='r')
        im = ax.imshow(grid_psi, origin='lower', extent=(0, 1, 0, 1))
        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.tick_params(labelsize=14)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.set_xlabel(r'$r_m$', fontsize=18)
        ax.set_ylabel(r'$z$', fontsize=18)
        #plt.savefig('psi.png', dpi=500)
        if save:
            plt.savefig(f'{self.PSR.name}_psi.png', bbox_inches='tight', dpi=400)

    def save_model(self, folder=''):
        if folder=='':
            path = ''
        else:
            path = folder + '/'
        torch.save(self.net.state_dict(), path +  f'{self.PSR.name}_model.pt')
        rs = np.linspace(0, 1, 100)
        hs = self.h0(rs)
        np.savetxt(path + f'{self.PSR.name}_h0.txt', hs)

    def load_psi(self, model_file='model.pt', h0_file='h0.txt'):
        self.net.load_state_dict(torch.load(model_file))
        hs = np.loadtxt(h0_file)
        rs = np.linspace(0, 1, 100)
        self.h0 = scipy.interpolate.interp1d(rs, hs, bounds_error=False, fill_value='extrapolate')


    def save_psi(self, r_grid=np.arange(0, 1, 0.001),  z_grid=np.arange(0, 5, 0.001)):
        """
        Save accelerating potentail on grid
        """
        rs = r_grid
        zs = z_grid
        grid_rs, grid_zs = np.meshgrid(rs, zs)
        rs = np.ravel(grid_rs).reshape(-1,1)
        zs = np.ravel(grid_zs).reshape(-1,1)
        psi = self.get_psi(rs, zs)
        grid_psi = psi.reshape(grid_rs.shape)
        np.savetxt(f'{self.PSR.name}_psis.txt', grid_psi)

    def show_residual(self, save=False, scale='linear'):
        rs = np.arange(0, 1, 0.01)
        zs = np.arange(0, 1, 0.01)
        fig, ax = plt.subplots()
        grid_rs, grid_zs = np.meshgrid(rs, zs)
        rs = np.ravel(grid_rs).reshape(-1,1)
        zs = np.ravel(grid_zs).reshape(-1,1)
        pt_rs = Variable(torch.from_numpy(rs).float(), requires_grad=True).to(device)
        pt_zs = Variable(torch.from_numpy(zs).float(), requires_grad=True).to(device)
        pt_res = self.PDEres(pt_rs, pt_zs)
        res = pt_res.data.cpu().numpy()
        res /= (4 * np.cos(self.PSR.chi) * self.Kpsi)
        res = np.abs(res)
        res = np.where(zs < self.h0(rs), res, np.nan)
        grid_res = res.reshape(grid_rs.shape)
        if scale == 'log':
            im = ax.imshow(grid_res, origin='lower', extent=(0, 1, 0, 1), norm=matplotlib.colors.LogNorm())
        else:
            im = ax.imshow(grid_res, origin='lower', extent=(0, 1, 0, 1))
        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.tick_params(labelsize=14)
        ax.set_xlabel(r'$r_m$', fontsize=18)
        ax.set_ylabel(r'$z$', fontsize=18)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        if save:
            plt.savefig('f{self.PSR.name}_psi_res.png', bbox_inches='tight', pdi=400)


class PotentialCalculator:
    def __init__(self, PSR:RadioPulsar, N_iter=20, h0:callable=lambda r: 1 + 0*r):
        self.PSR = PSR
        self.N_iter = N_iter
        self.h_func = h0
        self.psi_solver = PoissonSolver(PSR, h0) 

    def find_accelerating_potential(self):
        rs = np.linspace(0.0, 0.99, 100) # r grid for interpolation
        zs = np.linspace(0.0, 5, 500) # z grid for interpolation
        h_mins_old = np.zeros(rs.shape[0]) # old gap height
        h_mins_new = self.h_func(rs) # new gap height
        for _ in range(self.N_iter): # main loop for iterative procedure            
            h_mins_old = h_mins_new.copy()
            self.psi_solver.set_h0(self.h_func) # setting new emission height
            self.psi_solver.optimize_psi(5, 500, 5000) # find approximate Poisson equation solution
            grid_rs, grid_zs = np.meshgrid(rs, zs)
            # acceleration potentail interpolation from PINN solution
            psi_interp = PotentialCalculator.scipy_interp_wrapper(scipy.interpolate.RegularGridInterpolator((rs, zs), self.psi_solver.get_psi(np.ravel(grid_rs).reshape(-1,1), np.ravel(grid_zs).reshape(-1,1)).reshape(grid_rs.shape).T, bounds_error=False, fill_value=None)) 
            h_mins_new = self.calc_h_mins(rs, psi_interp) # new gap height
            h_mins_new = 0.2 * h_mins_new + 0.8 * h_mins_old # mixing with old gap height to increase convergence
            h_func = scipy.interpolate.interp1d(rs, h_mins_new, bounds_error=False, fill_value='extrapolate') # gap height interpolation
            plt.plot(np.linspace(0, 1, 100), h_func(np.linspace(0, 1, 100)))
            print(f"{_}-th iteration done")
            print(f"Current height at x=0.5 is {h_func(0.5)}")
            eps = np.sum(np.abs(h_mins_new - h_mins_old)) /  np.sum(np.abs(h_mins_new)) # rellative change of the emission height
            print(eps)
            if eps < 0.02:
                break
        psi_interp_final = PotentialCalculator.scipy_interp_wrapper(scipy.interpolate.RegularGridInterpolator((rs, zs), self.psi_solver.get_psi(np.ravel(grid_rs).reshape(-1,1), np.ravel(grid_zs).reshape(-1,1)).reshape(grid_rs.shape).T, bounds_error=False, fill_value=None)) 
        h_func_final = scipy.interpolate.interp1d(rs, h_mins_new, bounds_error=False, fill_value='extrapolate') # gap height interpolation
        return psi_interp_final, h_func_final
    
    @staticmethod
    def scipy_interp_wrapper(interp):
        """
        convert scipy 2D interpolation object into function with two separate arguments
        """
        def new_func(r, z):
            return interp((r, z))
        return new_func
    
    def calc_h_mins(self, rs, psi_interpolation):
        hs = np.zeros(shape=(rs.shape[0]))
        for i in range(rs.shape[0]):
            hs[i] = self.calc_h_min(rs[i], psi_interpolation)
        return hs
    
    def calc_h_min(self, r, psi_func):
        """
        returns gap height in polar cap radii  
        """
        r_max = 10 # maximum possible gap height (upper boundary of the minimum search)
        res = scipy.optimize.minimize_scalar(lambda x: self.l_tot_curv(r, x, psi_func), bounds=(0.01, r_max), method='bounded')
        if not res.success:
            print(f'minimisation has not converged at {r}')
        if self.l_tot_curv(r, res.x, psi_func) > 10.0:
                return r_max
        return self.l_tot_curv(r, res.x, psi_func)
    
    def l_tot_curv(self, r, l_rad, psi):
        return l_rad + self.l_gamma(r, l_rad, psi) 
    
    def l_gamma(self, r, h, psi):
        """
        gamma quantum free path length in polar cap radii
        """
        return 16 / 9 / self.PSR.Lambda * self.PSR.Rc(r, 1, units='cm') / self.PSR.R0 * constants.Bcr / self.PSR.B_surf12 / 1e12 * self.PSR.Rc(r, 1, units='cm') / constants.e_lambda_bar / self.gamma_e(r, h, psi)**3
    
    def gamma_e(self, r, h, psi_func):
        """
        primary particel gamma factor
        """
        return self.PSR.OmegaB * self.PSR.Omega * self.PSR.R0**2 / 2 / constants.c**2 * np.abs(psi_func(r, h))
        

    

