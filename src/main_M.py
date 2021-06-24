#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 09:39:38 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package
from numpy import linalg as LA

# Import user-defined scripts
from analysis_functions import musum
from parameters import N, M, I, n_t, mu_mag;

#%% Load in data from run

data_mu = np.load('run_data_mu_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
mu_int = data_mu['a'];

M_int = np.zeros((n_t,3));
M_mag = np.zeros((n_t));

for i in range(0,n_t):
    
    for j in range(0,3):
        
        M_int[i,j] = N*np.mean(mu_int[:,3*i+j]);
    
    M_mag[i] = LA.norm(M_int[i,:])

#%% Calculate the overall dipole moment autocorrelation function

M_acf = np.zeros((n_t,1));

print('Calculation the overall dipole moment auto-correlation function...')
M_acf = musum(M_int, N*M, n_t);
M_acf = M_acf/(M_acf[0,0]);
print('Done with the calculation.')

#%% Save data in compressed file

np.savez_compressed('M_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=M_acf[:,0]);