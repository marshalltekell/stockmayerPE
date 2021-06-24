#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 09:39:38 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package

# Import user-defined scripts
#from cython_functions_C1 import MUAC
import cfunctions;

# Import parameters specific to simulation
from parameters import M, N, I, n_t, mu_mag;

#%% Load in data from run

data_mu = np.load('run_data_mu_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
mu_int = data_mu['a'];
mu_int = mu_int.flatten();

#%% Calculate the dipole autocorrelation function

muac = np.zeros((n_t), dtype=np.double);
ct = np.arange(1,n_t+1,1, dtype=np.intc);
ct = ct[::-1];
ct[0] = n_t-1;
ct = ct.flatten();

print('Calculation the dipole auto-correlation function...')
cfunctions.MUAC(muac, mu_int, ct, N*M, n_t);
muac = muac/(mu_mag**2);
print('Done with the calculation.')

#%% Save data in compressed file

np.savez_compressed('MUAC_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=muac);