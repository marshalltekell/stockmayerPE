#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:46:45 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package

# Import user-defined functions from separate scripts
from cython_functions_C1 import duim
from parameters import N, M, I, n_t, mu_mag;

#%% Load in data from Terremoto zip and calculate MSD for the cations
data_ru = np.load('run_data_ru_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
ruc_int = data_ru['b'];
rua_int = data_ru['c'];
ruI_int = np.concatenate((ruc_int, rua_int));

# Initialize arrays
xind = np.arange(0,3*n_t,3);
yind = np.arange(1,3*n_t,3);
zind = np.arange(2,3*n_t,3);

x_int = ruI_int[:,xind];
y_int = ruI_int[:,yind];
z_int = ruI_int[:,zind];

#%% Calculate MSD

# Calc MSD
print("Calculating the degree of un-correlated ion motion.")
r_av_n, r_av_d, a_t, count = duim(x_int, y_int, z_int, I, n_t);
print("Done with MSD.")

#%% Save data in compressed file

np.savez_compressed('DUIM_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=r_av_n, b=r_av_d, c=a_t[:,0]);