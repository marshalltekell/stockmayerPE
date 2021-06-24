#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 09:39:38 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package

# Import user-defined scripts
from analysis_functions import musum_log
from functions import update_progress

from parameters import N, M, I, n_tlog, mu_mag, numloop;

#%% Load in data from run

data_mu_log = np.load('mu_ELEC_LOG.npz');  # a=mu_int
mu_full = data_mu_log['a'];

#%% Calculate the overall dipole for the simulation box

# Initialize overall array
M_full = np.zeros((n_tlog,3*numloop))

print('Calculating the overall dipole moments...')

for i in range(0, numloop):
    
    # Initialize arrays
    M_int = np.zeros((n_tlog,3));
    M_mag = np.zeros((n_tlog));
    
    # Get data for this loop
    mu_int = mu_full[i*N*M:(i+1)*N*M,:];
    
    for j in range(0,n_tlog):
        
        for k in range(0,3):
            
            M_int[j,k] = N*np.mean(mu_int[:,3*j+k]);
    
    # Add data to overall array
    M_full[:,i*3:(i+1)*3] = M_int

print('Done with the calculation.')
    
#%% Calculate the overall dipole moment autocorrelation function

M_acf_full = np.zeros((n_tlog,numloop));
M_acf_t = np.zeros((n_tlog));

print('Calculation the overall dipole moment auto-correlation function...')

for i in range(0,numloop):
    
    # Get data for this loop
    M_int = M_full[:,i*3:(i+1)*3];

    # Initialize the array
    M_acf = np.zeros((n_tlog));
    
    # Perform the calculation
    M_acf = musum_log(M_int, N, n_tlog);
    
    # Normalize to first point
    M_acf = M_acf/(M_acf[0]);
    
    # Add data to overall array
    M_acf_full[:,i] = M_acf;
    
    # Add data to average array
    M_acf_t += M_acf;
    
    update_progress(i/numloop);
    
    #M_acf = M_acf/(M_acf[0,0]);
    
# Find average
M_acf_t = M_acf_t/numloop;
    
print('Done with the calculation.')


#%% Save data in compressed file

np.savez_compressed('M_log_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=M_acf_t);