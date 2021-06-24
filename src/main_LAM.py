#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:46:45 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package
import cfunctions
import matplotlib.pyplot as plt;

# Import parameters for simulation
from parameters import N, M, I, n_t, mu_mag, time;

#%% Additional parameters

n_t_ext = np.int(50e6/1e4);
        
#%% Load in data

# First run
data_ru = np.load('run_data_ru_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
ru_int = data_ru['a'];
ruc_int = data_ru['b'];
rua_int = data_ru['c'];

# Second run
'''
data_ru2 = np.load('run_data_ru_ext_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
ru_int2 = data_ru2['a'];
ruc_int2 = data_ru2['b'];
rua_int2 = data_ru2['c'];

# Third run
data_ru3 = np.load('run_data_ru_ext2_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
ru_int3 = data_ru3['a'];
ruc_int3 = data_ru3['b'];
rua_int3 = data_ru3['c'];
'''

# Initialize arrays
xind = np.arange(0,3*(n_t),3);
yind = np.arange(1,3*(n_t),3);
zind = np.arange(2,3*(n_t),3);

# Stitch together
#ru_int = np.concatenate((ru_int, ru_int2, ru_int3), axis = 1);
#rua_int = np.concatenate((rua_int, rua_int2, rua_int3), axis = 1);
#ruc_int = np.concatenate((ruc_int, ruc_int2, ruc_int3), axis = 1);

# Polymer
xu_int = ru_int[:,xind];
yu_int = ru_int[:,yind];
zu_int = ru_int[:,zind];

# Cation
xuc_int = ruc_int[:,xind];
yuc_int = ruc_int[:,yind];
zuc_int = ruc_int[:,zind];

# Anion
xua_int = rua_int[:,xind];
yua_int = rua_int[:,yind];
zua_int = rua_int[:,zind];

# Ions
xui_int = np.concatenate((xuc_int, xua_int), axis=0);
yui_int = np.concatenate((yuc_int, yua_int), axis=0);
zui_int = np.concatenate((zuc_int, zua_int), axis=0);

#%% Create array of charges

qp = np.ones((I), dtype=np.intc);
qm = np.ones((I), dtype=np.intc);
qm = -qm;
q = np.concatenate((qp,qm));

#%% Create count array

ct = np.arange(1,n_t+1,1, dtype=np.intc);
ct = ct[::-1];
ct[0] = n_t-1;
ct = np.ascontiguousarray(ct, dtype=np.intc);

#ct = np.zeros((n_t), dtype=np.intc)

#%%

# Initialize arrays
rt1 = np.zeros((n_t), dtype=np.double);
rt2 = np.zeros((n_t), dtype=np.double);
rt3 = np.zeros((n_t), dtype=np.double);
rt4 = np.zeros((n_t), dtype=np.double);
r = np.zeros((n_t), dtype=np.double);
#ct = np.zeros((n_t), dtype=np.intc);
x1D = np.zeros((2*I*(n_t)), dtype=np.double);
y1D = np.zeros((2*I*(n_t)), dtype=np.double);
z1D = np.zeros((2*I*(n_t)), dtype=np.double);

# Make 2D position arrays 1D
for i in range(0, 2*I):

    for j in range(0, n_t):
        
        x1D[i*(n_t) + j] = xui_int[i,j];
        y1D[i*(n_t) + j] = yui_int[i,j];
        z1D[i*(n_t) + j] = zui_int[i,j];
        
#%%

# Calculate MSD (updates arrays rt#, rt, and ct in place)
print('Calculating ionic conductivity (full).')
cfunctions.lam(rt1, r, ct, x1D, y1D, z1D, q, 2*I, n_t)
print('Calculating ionic conductivity (i!=j).')
cfunctions.lam_cross1(rt2, r, ct, x1D, y1D, z1D, q, 2*I, n_t)
print('Calculating ionic conductivity (i!=j, qi = qj).')
cfunctions.lam_cross2(rt3, r, ct, x1D, y1D, z1D, q, 2*I, n_t)
print('Calculating ionic conductivity (i!=j, qi != qj).')
cfunctions.lam_cross3(rt4, r, ct, x1D, y1D, z1D, q, 2*I, n_t)
#r = np.divide(r, ct, out=np.zeros_like(r), where=ct!=0);
print('Done calculating ionic conductivity.')


#%% Save data in compressed file

np.savez_compressed('LAM_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=rt1, b=rt2, c=rt3, d=rt4);


