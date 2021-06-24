#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:06:05 2020

@author: marshalltekell
"""

# Import packages and scripts
import numpy as np
from classes import vanHove

# Import parameters for simulation
from parameters import m, n_tlog, N, M, I, N_Av, numloop;
from param_NVT import L, V;

#%% Read in parameters that specify the conditions of the the run

# Component densities (g cm-3)
rhop = (10**24*N*M*m[0])/(N_Av*V);
rhoc = (10**24*I*m[1])/(N_Av*V);
rhoa = (10**24*I*m[2])/(N_Av*V);

# Choose number of bins
nbins= 1000;
r_RDF = np.linspace(0,L/2, nbins);

#%% Load in data from Terremoto zip

data_r = np.load('r_ELEC_LOG.npz');  # a=r_int, b=rc_int, c=ra_int

# Assign position data
r_full = data_r['a'];
rc_full = data_r['b'];
ra_full = data_r['c'];

# Indices
xind = np.arange(0,3*n_tlog,3);
yind = np.arange(1,3*n_tlog,3);
zind = np.arange(2,3*n_tlog,3);

# Polymer
x_full = r_full[:,xind];
y_full = r_full[:,yind];
z_full = r_full[:,zind];

# Cation
xc_full = rc_full[:,xind];
yc_full = rc_full[:,yind];
zc_full = rc_full[:,zind];

# Anion
xa_full = ra_full[:,xind];
ya_full = ra_full[:,yind];
za_full = ra_full[:,zind];

# Initialize arrays: CB
xcb_full = np.zeros(((I+N*M)*100,n_tlog));
ycb_full = np.zeros(((I+N*M)*100,n_tlog));
zcb_full = np.zeros(((I+N*M)*100,n_tlog));

# Initialize arrays: CA
xca_full = np.zeros(((I+I)*100,n_tlog));
yca_full = np.zeros(((I+I)*100,n_tlog));
zca_full = np.zeros(((I+I)*100,n_tlog));

# Initialize arrays: AB
xab_full = np.zeros(((I+N*M)*100,n_tlog));
yab_full = np.zeros(((I+N*M)*100,n_tlog));
zab_full = np.zeros(((I+N*M)*100,n_tlog));

for i in range(0, numloop):
    
    # CB
    xcb_full[i*(N*M+I):i*(N*M+I)+I,:] = xc_full[i*I:(i+1)*I,:];
    xcb_full[i*(N*M+I)+I:(i+1)*(N*M+I),:] = x_full[i*N*M:(i+1)*N*M,:];
    
    ycb_full[i*(N*M+I):i*(N*M+I)+I,:] = yc_full[i*I:(i+1)*I,:];
    ycb_full[i*(N*M+I)+I:(i+1)*(N*M+I),:] = y_full[i*N*M:(i+1)*N*M,:];
    
    zcb_full[i*(N*M+I):i*(N*M+I)+I,:] = zc_full[i*I:(i+1)*I,:];
    zcb_full[i*(N*M+I)+I:(i+1)*(N*M+I),:] = z_full[i*N*M:(i+1)*N*M,:];
    
    # CA
    xca_full[i*(I+I):i*(I+I)+I,:] = xc_full[i*I:(i+1)*I,:];
    xca_full[i*(I+I)+I:(i+1)*(I+I),:] = xa_full[i*I:(i+1)*I,:];
    
    yca_full[i*(I+I):i*(I+I)+I,:] = yc_full[i*I:(i+1)*I,:];
    yca_full[i*(I+I)+I:(i+1)*(I+I),:] = ya_full[i*I:(i+1)*I,:];
    
    zca_full[i*(I+I):i*(I+I)+I,:] = zc_full[i*I:(i+1)*I,:];
    zca_full[i*(I+I)+I:(i+1)*(I+I),:] = za_full[i*I:(i+1)*I,:];
    
    # AB
    xab_full[i*(N*M+I):i*(N*M+I)+I,:] = xa_full[i*I:(i+1)*I,:];
    xab_full[i*(N*M+I)+I:(i+1)*(N*M+I),:] = x_full[i*N*M:(i+1)*N*M,:];
    
    yab_full[i*(N*M+I):i*(N*M+I)+I,:] = ya_full[i*I:(i+1)*I,:];
    yab_full[i*(N*M+I)+I:(i+1)*(N*M+I),:] = y_full[i*N*M:(i+1)*N*M,:];
    
    zab_full[i*(N*M+I):i*(N*M+I)+I,:] = za_full[i*I:(i+1)*I,:];
    zab_full[i*(N*M+I)+I:(i+1)*(N*M+I),:] = z_full[i*N*M:(i+1)*N*M,:];

# Flatten
xcb_full = xcb_full.flatten();
ycb_full = ycb_full.flatten();
zcb_full = zcb_full.flatten();

xca_full = xca_full.flatten();
yca_full = yca_full.flatten();
zca_full = zca_full.flatten();

xab_full = xab_full.flatten();
yab_full = yab_full.flatten();
zab_full = zab_full.flatten();

#%% Create instance of class
### KEY for g(r,t):
        
# g0(r,t) is cation/bead
# g1(r,t) is cation/anion
# g2(r,t) is anion/bead
        
g_rt = [];

# Initialize g(r) array
grt = np.zeros((nbins*n_tlog));


g_rt.append(vanHove(xcb_full, ycb_full, zcb_full, I, N*M, nbins, n_tlog, L));
g_rt.append(vanHove(xca_full, yca_full, zca_full, I, I, nbins, n_tlog, L));
g_rt.append(vanHove(xab_full, yab_full, zab_full, I, N*M, nbins, n_tlog, L));

# Initialize g(r) array
grt = [];

# Get g(r,t) for each class instance
grt.append(g_rt[0].VH1(V, numloop));
grt.append(g_rt[1].VH1(V, numloop));
grt.append(g_rt[2].VH1(V, numloop));

#%% Save data to compressed file

np.savez_compressed('VH_log_{}_{}_{}'.format(N, M, I), a=grt);
                    
                
