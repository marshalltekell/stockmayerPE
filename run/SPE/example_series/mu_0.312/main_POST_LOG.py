#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:34:22 2020

@author: marshalltekell
"""

#%% Import packages and scripts

# Import packages
import numpy as np #python array manipulation package

# Misc
from functions import parse_LAMMPs_out_int
from functions import update_progress

# Import values from parameters.py
from parameters import N, M, I, n_tlog, numloop;

#%% Additional parameters specific to this file

n_t = n_tlog;

#%% Parse the unfolded coordinates
        
# Initialize and allocate arrays for parse
names_in = ['ID','x','y','z']

# Initialize and allocate arrays for per-atom values

# Unfolded positions
ru_full = np.zeros((numloop*N*M,3*n_t));
ruc_full = np.zeros((numloop*I,3*n_t));
rua_full = np.zeros((numloop*I,3*n_t));

# Folded positions
r_full = np.zeros((numloop*N*M,3*n_t));
rc_full = np.zeros((numloop*I,3*n_t));
ra_full = np.zeros((numloop*I,3*n_t));

# Dipoles
mu_full = np.zeros((numloop*N*M,3*n_t));

print("Parsing output from simulation...")

# Parse per-atom dumps

for i in range(0,numloop):
    
    # Reset arrays (unfolded positions)
    ru_int = np.zeros((N*M,3*n_t));
    ruc_int = np.zeros((I,3*n_t));
    rua_int = np.zeros((I,3*n_t));
    
    # Reset arrays (folded positions)
    r_int = np.zeros((N*M,3*n_t));
    rc_int = np.zeros((I,3*n_t));
    ra_int = np.zeros((I,3*n_t));
    
    # Dipoles
    mu_int = np.zeros((N*M,3*n_t));
    
    for j in range(0,n_t):
        
        # Positions (folded)
        temp_df = parse_LAMMPs_out_int('r.{:d}.out'.format(i+1), M*N+2*I, 9, names_in, j)
        r_int[:,j*3:j*3+3] = temp_df[0:M*N,:];
        rc_int[:,j*3:j*3+3] = temp_df[M*N:M*N+I,:];
        ra_int[:,j*3:j*3+3] = temp_df[M*N+I:M*N+2*I,:];
        
        # Positions (unfolded)
        temp_df = parse_LAMMPs_out_int('ru.{:d}.out'.format(i+1), M*N+2*I, 9, names_in, j)
        ru_int[:,j*3:j*3+3] = temp_df[0:M*N,:];
        ruc_int[:,j*3:j*3+3] = temp_df[M*N:M*N+I,:];
        rua_int[:,j*3:j*3+3] = temp_df[M*N+I:M*N+2*I,:];
        
        # Dipoles
        temp_df = parse_LAMMPs_out_int("mu.{:d}.out".format(i+1), M*N+2*I, 9, names_in, j)
        mu_int[:,j*3:j*3+3] = temp_df[0:M*N,:];
        
        # Update progress to output
        update_progress((i*n_t + j)/(numloop*n_t))
        
        
    # Add data to master array (unfolded positions)
    ru_full[i*N*M:(i+1)*N*M,:] = ru_int;
    ruc_full[i*I:(i+1)*I,:] = ruc_int;
    rua_full[i*I:(i+1)*I,:] = rua_int;
    
    # Add data to master array (folded positions)
    r_full[i*N*M:(i+1)*N*M,:] = r_int;
    rc_full[i*I:(i+1)*I,:] = rc_int;
    ra_full[i*I:(i+1)*I,:] = ra_int;
    
    # Dipoles
    mu_full[i*N*M:(i+1)*N*M,:] = mu_int;
    
#%% Save all of the data to a compressed file

np.savez_compressed('r_ELEC_LOG', a=r_full, b=rc_full, c=ra_full);  
np.savez_compressed('ru_ELEC_LOG', a=ru_full, b=ruc_full, c=rua_full);
np.savez_compressed('mu_ELEC_LOG', a=mu_full);