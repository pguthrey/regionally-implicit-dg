import matplotlib.pyplot as plt
import numpy as np
import time
import math
import dg_error as err
import dg_projection as proj
import dg_ode as ode

import scipy
from scipy import fftpack
from scipy import interpolate

import matplotlib
import seaborn as sns

import ipywidgets as widgets
from ipywidgets import interactive

def solve_once(Gamma,kappa,guess):
    r_bins = 128
    #--------------------------------------------------------------------------------------------
    # run parameters
    num_iterations = 1000
    epsilon = 0.95
    alpha = 1.0
    
    # r binning
    R_max = 30.0
    del_r = R_max/r_bins
    # error is minimized when using bin centers
    r_array = np.linspace(del_r/2,R_max-del_r/2,r_bins)
    energy_iter = np.zeros(num_iterations)
    pressure_iter = np.zeros(num_iterations)
    
    # k binning
    del_k = np.pi/((r_bins+1)*del_r)
    K_max = del_k*r_bins
    k_array = np.linspace(del_k/2,K_max-del_k/2,r_bins)
    
    # precompute factors
    fact_r_2_k = 2*np.pi*del_r
    fact_k_2_r = del_k/(4.*np.pi**2)
    dimless_dens = 3./(4*np.pi)
    
    #--------------------------------------------------------------------------------------------
    
    # f(k) = FT[f(r)]
    def FT_r_2_k(input_array):
        from_dst = fact_r_2_k*fftpack.dst(r_array*input_array)
        return from_dst/k_array
    
    # f(r) = IFT[f(k)]
    def FT_k_2_r(input_array):
        from_idst = fact_k_2_r*fftpack.idst(k_array*input_array)
        return from_idst/r_array
    
    #--------------------------------------------------------------------------------------------
    #Gamma = params.Gamma
    #kappa = params.kappa
    #--------------------------------------------------------------------------------------------
    
    # Returns a tuple of arrays with [u(r),u_s(r)]; it is the short-ranged u_s(r) that is used in the HNC solver.
    # Note that u(r) itself is not needed for HNC, but is used for computing the energy, which is used to test convergence.
    # This function assumes global variables Gamma, kappa and alpha.
    def u(r_in):
        return [Gamma*np.exp(-kappa*r_in)/r_in,Gamma*np.exp(-(alpha + kappa)*r_in)/r_in]
    
    # when du/dr is known, precompute it here; only needed for computing the pressure, not for the HNC solution itself
    def deriv_u_r(r_in):
        return -Gamma*np.exp(-kappa*r_in)*(1 + kappa*r_in)/r_in**2
    
    # The initial condition is given through c(k), defined here. This should be consistent with the other definitions
    # in this cell.
    def initial_c_k(k_in):
        return -4*np.pi*Gamma/(k_array**2 + kappa**2)
    
    # Long range part of the potential in Fourier space. This could be computed with the Fourier routines, but
    # I chose to precopmute things analytically when possible.
    def u_long_l(k_in):
        return 4*np.pi*Gamma*(alpha**2 + 2*alpha*kappa)/((k_array**2 + kappa**2)*(k_array**2+(alpha + kappa)**2))
    
    #--------------------------------------------------------------------------------------------
    u_r, u_s_r = u(r_array)
    d_u_d_r = deriv_u_r(r_array)
    c_k = initial_c_k(k_array)
    u_l_k = u_long_l(k_array)
    c_s_k = c_k + u_l_k
    S_k = 1+dimless_dens*c_k/(1-dimless_dens*c_k)
    # start the iteration loop
    for iteration in np.arange(num_iterations):
        gamma_s_k = (dimless_dens*c_s_k*c_k - u_l_k)/(1-dimless_dens*c_k)
        gamma_s_r = FT_k_2_r(gamma_s_k)
        g_r = np.exp(gamma_s_r - u_s_r)
        energy_iter[iteration] = 1.5*del_r*np.sum(r_array**2*u_r*(g_r-1))
        pressure_iter[iteration] = -0.5*del_r*np.sum(r_array**3*d_u_d_r*(g_r-1))
        new_c_s_r = g_r - 1 - gamma_s_r
        c_s_k = FT_r_2_k(new_c_s_r)
        c_k = c_s_k - u_l_k
    
    g_data = np.array([r_array,g_r]).T  # g(r)
    c_data = np.array([k_array,c_k]).T  # c(k)
    S_data = np.array([k_array,S_k]).T  # S(k)

    #Let C(|x-y|) = c(|x-y|)|x-y|  
    # C(x)        = c(|x|)|x|
    # dFdn        = int n(x)C(|x-y|) dy 
    # F[dFdn](k)  = F[n](k)F[C](k)    
    # dFdn(x)     = F^-1[ F[n](k)F[C](k) ](x)
    c_x = k_array         # x=k values
    c_f = FT_k_2_r(c_k)   # c real space
    c_F = FT_r_2_k(c_f)   # c in fourier space
    c_fr = c_f*c_x        # c*r in real space
    c_Fr = FT_r_2_k(c_fr) # c*r in fourier space

    return c_x,c_f,c_F,c_fr,c_Fr,energy_iter[-1],pressure_iter[-1]

