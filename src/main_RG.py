#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:06:48 2020

@author: marshalltekell

Calculate radius of gyration and end-to-end distances
"""
# Equation of state of the Lennard-Jones fluid

# Import packages
import numpy as np #python array manipulation package

# Import user-defined functions from separate scripts
from analysis_functions import calc_RG
from analysis_functions import calc_RG2
from analysis_functions import update_progress

# Import parameters for simulation
from parameters import n_t, N, M, I, mu_mag;
from param_NVT import L;

#%% Load in data

data_r = np.load('run_data_r_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int
r_int = data_r['a'];
   
#%% Compute the radius of gyration

# Initialize arrays
rg = np.zeros((M,7*n_t));
rg_av = np.zeros((n_t,7));
av_rg = np.zeros((7));
std_rg = np.zeros((7));

print('Calculating radius of gyration and end-to-end distances.')
update_progress(0/n_t);

for i in range(0, n_t):
    
    # Get data for this time sample
    r = r_int[:,i*3:(i+1)*3];
    
    # If the chain crosses the edge of the box, add or subtract one
    count = calc_RG2(N, r, L, M);
    count = count*L;
    
    # Unfold position data
    r += count;
    
    # Get rg, rg4, r2, r4, and gyration tensor eigenvalues
    rg[:,i*7:(i+1)*7], w, s3 = calc_RG(N, M, L, r);
    
    # Calculate the averages
    for j in range(0, 7):
        
        rg_av[i,j] = np.mean(rg[:,7*i+j])
        
    # Normalize gyration tensor values
    rg_av[i,5] = rg_av[i,5]/rg_av[i,4];
    rg_av[i,6] = rg_av[i,6]/rg_av[i,4];
    rg_av[i,4] = 1.0;
    
    update_progress((i+1)/n_t);
    
# Calculate ensemble averages and variances
for i in range(0,7):
    
    av_rg[i] = np.mean(rg_av[:,i]);
    std_rg[i] = np.std(rg_av[:,i]);
    
#%% Save to npz

np.savez_compressed('RG_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=rg_av, b=av_rg, c=std_rg, d=w, e=s3);