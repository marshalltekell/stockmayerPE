#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 09:39:38 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package

# Import user-defined scripts
#from cython_functions_C1 import MUAC_log
import cfunctions;

# Import parameters for the simulation
from parameters import N, M, I, mu_mag, numloop, n_tlog


#%% Load in data from run

data_mu_log = np.load('mu_ELEC_LOG.npz');  # a=mu_int;
mu_full = data_mu_log['a'];

#%% Calculate the dipole autocorrelation function

print('Calculation the dipole auto-correlation function...')

# Initialize arrays
muac_full = np.zeros((n_tlog,numloop))
muac_t = np.zeros((n_tlog))

for i in range(0, numloop):
    
    # Initialize arrays
    muac = np.zeros((n_tlog));
 
    # Get data for current loop
    mu_int = mu_full[i*N*M:(i+1)*N*M,:];
    mu_int = mu_int.flatten();
    
    # Perform calculation
    cfunctions.MUAC_log(muac, mu_int, N*M, n_tlog);
    muac = muac/(mu_mag**2);
    
    # Assign current run to overall
    muac_full[:,i] = muac;
    
    # Add run to total
    muac_t += muac;
    
# Find the average
muac_t = muac_t/numloop;
    
print('Done with the calculation.')

#%% Save data in compressed file

np.savez_compressed('MUAC_log_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=muac_t);