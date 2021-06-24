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
from parameters import m, n_t, N, M, I, N_Av;
from param_NVT import L, V;

#%% Parameters needed for van Hove function calculation

# Componen density (g cm-3)
rhop = (10**24*N*M*m[0])/(N_Av*V);
rhoc = (10**24*I*m[1])/(N_Av*V);
rhoa = (10**24*I*m[2])/(N_Av*V);

# Choose number of bins
nbins= 1000;
r_RDF = np.linspace(0,L/2, nbins);

#%% Load in data from Terremoto zip

data_r = np.load('run_data_r_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int

# Assign position data
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

#%% Create instance of class
### KEY for g(r,t):
        
# g0(r,t) is cation/bead
# g1(r,t) is cation/anion
# g2(r,t) is anion/bead
        
g_rt = [];

g_rt.append(vanHove(xcb_int, ycb_int, zcb_int, I, N*M,  nbins, n_t, L));
g_rt.append(vanHove(xca_int, yca_int, zca_int, I, I,  nbins, n_t, L));
g_rt.append(vanHove(xab_int, yab_int, zab_int, I, N*M,  nbins, n_t, L));

# Initialize g(r) array
grt = [];

# Get g(r,t) for each class instance
grt.append(g_rt[0].VH2(V, 10));
grt.append(g_rt[1].VH2(V, 10));
grt.append(g_rt[2].VH2(V, 10));
#grt[:,2*10:3*10] = g_rt[2].VH0(m[0], rhop);

#g_rt[0].VH0(m[0], rhop);
#g_rt[1].VH0(m[2], rhoa);
#g_rt[2].VH0(m[0], rhop);

#%% Save data to compressed file

np.savez_compressed('VH_{}_{}_{}'.format(N, M, I), a=grt);
