#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:06:48 2020

@author: marshalltekell

NVE MD with only Lennard-Jones interactions
"""
# Equation of state of the Lennard-Jones fluid

# Import packages
import numpy as np #python array manipulation package
from analysis_functions import gen_q
from analysis_functions import update_progress

# Import parameters for simulation
from parameters import N, M, I, n_t, mu_mag;
from param_NVT import L;
from classes import SSF

#%% Additional parameters needed for this specific calculation

# Prep for q vectors
dk = (2*np.pi)/L;
low = dk;
high = 5;
q_mag = np.arange(low,high,dk);
num_p = q_mag.shape[0];

#%% Load in data from Terremoto zip

data_r = np.load('run_data_r_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
r_int = data_r['a'];
rc_int = data_r['b'];
ra_int = data_r['c'];

r_full = np.concatenate((r_int, rc_int, ra_int));
rcb_int = np.concatenate((rc_int, r_int));
rca_int = np.concatenate((rc_int, ra_int));
rab_int = np.concatenate((ra_int, r_int));

# Initialize arrays
xind = np.arange(0,3*n_t,3);
yind = np.arange(1,3*n_t,3);
zind = np.arange(2,3*n_t,3);

# Full
x_full = r_full[:,xind];
x_full = x_full.flatten();
y_full = r_full[:,yind];
y_full = y_full.flatten();
z_full = r_full[:,zind];
z_full = z_full.flatten();

# Polymer
x_int = r_int[:,xind];
x_int = x_int.flatten();
y_int = r_int[:,yind];
y_int = y_int.flatten();
z_int = r_int[:,zind];
z_int = z_int.flatten();

# Cation
xc_int = rc_int[:,xind];
xc_int = xc_int.flatten();
yc_int = rc_int[:,yind];
yc_int = yc_int.flatten();
zc_int = rc_int[:,zind];
zc_int = zc_int.flatten();

# Anion
xa_int = ra_int[:,xind];
xa_int = xa_int.flatten();
ya_int = ra_int[:,yind];
ya_int = ya_int.flatten();
za_int = ra_int[:,zind];
za_int = za_int.flatten();

# Cation-bead
xcb_int = rcb_int[:,xind];
xcb_int = xcb_int.flatten();
ycb_int = rcb_int[:,yind];
ycb_int = ycb_int.flatten();
zcb_int = rcb_int[:,zind];
zcb_int = zcb_int.flatten();

# Cation-anion
xca_int = rca_int[:,xind];
xca_int = xca_int.flatten();
yca_int = rca_int[:,yind];
yca_int = yca_int.flatten();
zca_int = rca_int[:,zind];
zca_int = zca_int.flatten();

# Anion-bead
xab_int = rab_int[:,xind];
xab_int = xab_int.flatten();
yab_int = rab_int[:,yind];
yab_int = yab_int.flatten();
zab_int = rab_int[:,zind];
zab_int = zab_int.flatten();

#%% Generate array of q vectors

print('Generating arrays of q vectors for calculation.')

# Initialize arrays
num_q = np.zeros((q_mag.shape[0]),dtype=np.intc);
q_run = np.zeros((300,3*num_p));

# Loop through the q magnitudes and collect vectors at each mag
for i in range(0, num_p):
    
    update_progress((i+1)/num_p)
    
    q_vec, num_q[i] = gen_q(q_mag[i], L);
        
    # Scan in values to master
    for j in range(0, num_q[i]):
        
        for k in range(0,3):
            
            q_run[j,3*i+k] = q_vec[j,k];
            
print('Done generating vectors.')

# Flatten wave vectors
q_run  = q_run.flatten();

#%% Create instances of class

## KEY for S(q):
# S0(q) is for the whole system
# S1(q) is for polymer beads
# S2(q) is for cations
# S3(q) is for anions
# S4(q) is for cation/bead pairs
# S5(q) is for cation/anion pairs
# S6(q) is for anion/bead pairs
        
# Create instances
S_q = [];

S_q.append(SSF(x_full, y_full, z_full, q_run, num_q, N*M+2*I, n_t,  num_p));
S_q.append(SSF(x_int, y_int, z_int, q_run, num_q, N*M, n_t,  num_p));
S_q.append(SSF(xc_int, yc_int, zc_int, q_run, num_q, I, n_t,  num_p));
S_q.append(SSF(xa_int, ya_int, za_int, q_run, num_q, I, n_t,  num_p));
S_q.append(SSF(xcb_int, ycb_int, zcb_int, q_run, num_q, N*M, 10,  num_p));
S_q.append(SSF(xca_int, yca_int, zca_int, q_run, num_q, I, 10,  num_p));
S_q.append(SSF(xab_int, yab_int, zab_int, q_run, num_q, N*M, 10,  num_p));

# Initialize g(r) array
S = np.zeros((num_p,len(S_q)));

# Get g(r) for each class instance
S[:,0] = S_q[0].SSF1();
S[:,1] = S_q[1].SSF1();
S[:,2] = S_q[2].SSF1();
S[:,3] = S_q[3].SSF1();
S[:,4] = S_q[4].SSF2(I);
S[:,5] = S_q[5].SSF2(I);
S[:,6] = S_q[6].SSF2(I);

#%% Save to compressed file

np.savez_compressed('SSF_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=S);