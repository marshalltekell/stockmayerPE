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
from analysis_functions import update_progress
import cfunctions;
#import matplotlib.pyplot as plt;

# Import necessary parameters from master list specific to this simulation
from parameters import M, N, n_t, e, o, oij, k_FENE, R0, K, L, mu_mag;
           
#%% Load in data from Terremoto zip

# Position data
data_r = np.load('run_data_r_lammps_ELEC.npz');  # a=r_int, b=rc_int

# Assign position data
r_int = data_r['a'];
rc_int = data_r['b'];

# Dipole moment data
data_mu = np.load('run_data_mu_lammps_ELEC.npz');  # a=mu_int

# Assign dipole data
mu_int = data_mu['a'];

#%% Calculate potential energy for each time step

'''
print('Calculating potential energy (FENE bonds).')

# Initialize array
u = np.zeros((n_t), dtype=np.double);

# Loop over the number of times in the production run
for i in range(0, n_t):
    
    update_progress(i/n_t);
    
    # Get positions and dipoles at this time step
    r = r_int[:,i*3:(i+1)*3];
    
    # Flatten array
    r = r.flatten();
    
    # Make it a np.double
    r = r.astype(np.double);
    
    # Initialize energy
    u_t = np.array([0.0], dtype=np.double);
    
    # Calculate g(r) (updates g in place)
    cfunctions.pair_fene(u_t, r, e, k_FENE, R0, L, N, M);
    
    # Store in overall array
    u[i] = u_t;
    
print('Done.')


print('Calculating potential energy (angles).')

# Initialize array
u2 = np.zeros((n_t), dtype=np.double);

# Loop over the number of times in the production run
for i in range(0, n_t):
    
    update_progress(i/n_t);
    
    # Get positions and dipoles at this time step
    r = r_int[:,i*3:(i+1)*3];
    
    # Flatten array
    r = r.flatten();
    
    # Make it a np.double
    r = r.astype(np.double);
    
    # Initialize energy
    u_t = np.array([0.0], dtype=np.double);
    
    # Calculate g(r) (updates g in place)
    cfunctions.pair_angle(u_t, r, K, L, N*M, M);
    
    # Store in overall array
    u2[i] = u_t;
    
print('Done.')

print('Calculating potential energy (dipole-dipole interactions).')

# Initialize array
u3 = np.zeros((n_t), dtype=np.double);

# Loop over the number of times in the production run
for i in range(0, n_t):
    
    update_progress(i/n_t);
    
    # Get positions and dipoles at this time step
    r = r_int[:,i*3:(i+1)*3];
    mu = mu_int[:,i*3:(i+1)*3];
    
    # Flatten array
    r = r.flatten();
    mu = mu.flatten();
    
    # Make it a np.double
    r = r.astype(np.double);
    mu = mu.astype(np.double);
    
    # Initialize energy
    u_t = np.array([0.0], dtype=np.double);
    
    # Calculate g(r) (updates g in place)
    cfunctions.pair_lj_dipole(u_t, r, mu, e, o, L, N*M);
    
    # Store in overall array
    u3[i] = u_t;
    
print('Done.')

'''
#%%
print('Calculating potential energy (charge-dipole interactions).')

# Initialize array
u4 = np.zeros((n_t), dtype=np.double);

# Loop over the number of times in the production run
for i in range(0, n_t):
    
    update_progress(i/n_t);
    
    # Get positions and dipoles at this time step
    r = r_int[:,i*3:(i+1)*3];
    rc = rc_int[0,i*3:(i+1)*3];
    mu = mu_int[:,i*3:(i+1)*3];
    
    # Flatten array
    r = r.flatten();
    rc = rc.flatten();
    mu = mu.flatten();
    
    # Make it a np.double
    r = r.astype(np.double);
    rc = rc.astype(np.double);
    mu = mu.astype(np.double);
    
    # Initialize energy
    u_t = np.array([0.0], dtype=np.double);
    
    # Calculate g(r) (updates g in place)
    cfunctions.pair_lj_dipole_cat(u_t, r, mu, rc, e, oij[1], L, N*M, -1.0)
    
    # Store in overall array
    u4[i] = u_t;
    
print('Done.')


#print(u2)

#%% Plot result to debug
'''
plt.scatter(time,u)
plt.show()
plt.scatter(time,u2)
plt.show()
'''
np.savez_compressed('PE_onean_cut_{:.3f}'.format(mu_mag), a=u4);