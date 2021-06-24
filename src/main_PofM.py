#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:06:48 2020

@author: marshalltekell

Find various pair distribution functions
"""
# Equation of state of the Lennard-Jones fluid

# Import packages
import numpy as np #python array manipulation package
#from classes import RDF
#from functions import pofn
#from pofn import pofn
from cfunctions import pofM
from analysis_functions import update_progress
#from analysis_functions import fcn_slow

# Import necessary parameters from master list specific to this simulation
from parameters import M, N, I, n_t, mu_mag;
from param_NVT import L, V;

#%% Additional paramters

lo = 2.0;
hi = 4.0;

#%% Load in data from Terremoto zip

data_r = np.load('run_data_r_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int

# Assign position data
r_int = data_r['a'];
rc_int = data_r['b'];
ra_int = data_r['c'];

# Concatenate data for functions
rcb_int = np.concatenate((rc_int, r_int));


#%% Perform calculation

# Get array to hold p[i] for all of the timesteps
p_t = np.zeros((I, n_t));

print('Calculating coordination statistics.')

# Loop over the number of times
for i in range(0, n_t):
        
    # Update progress
    update_progress(i/n_t);
    
    # Get positions for this time step
    r = rcb_int[:,i*3:(i+1)*3];
    r = r.flatten();
    
    # Get new array
    p = np.zeros((I),dtype=np.intc);
    
    # Get instantaneous coordination number
    pofM(p, r, I, N*M, lo, hi, L, M);
    
    # Assign to total array
    p_t[:,i] = p;
    
print('Done calculating coordination statistics.')
    
#%% Flatten p_t array
p_t = p_t.flatten();

# Sort integers from low to high
p_t = np.sort(p_t);

# Get min and max
p_min = np.int(np.min(p_t));
p_max = np.int(np.max(p_t));

# Count number of entries for each
cn_run = np.arange(p_min, p_max+1, 1, dtype=np.int);
cn_ct = np.zeros((cn_run.shape[0], 1));

# Loop over the entire array
for i in range(0, p_t.shape[0]):
    
    # Loop over the integer options
    for j in range(0, cn_run.shape[0]):
        
        if (p_t[i] == cn_run[j]):
            
            cn_ct[j,0] += 1;
            
# Normalize
p_end = cn_ct/(n_t*I);

#%% Save result

np.savez('PofM_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=p_end, b=cn_run)
    
    

