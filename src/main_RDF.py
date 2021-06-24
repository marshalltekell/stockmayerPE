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
from classes import RDF
from analysis_functions import fcn_slow

# Import necessary parameters from master list specific to this simulation
from parameters import M, N, I, n_t, mu_mag;
from param_NVT import L, V;


#%% Parameters needed specifically for RDF calculations

# Density (by component, g cm-3)
#rhop = (10**24*N*M*m[0])/(N_Av*V);
rhop = (N*M)/V;
#rhoc = (10**24*I*m[1])/(N_Av*V);
rhoc = I/V;
#rhoa = (10**24*I*m[2])/(N_Av*V);
rhoa = I/V;
rho = (N*M+2*I)/V;

# Choose number of bins
nbins= 1000;
r_RDF = np.linspace(0, L/2, nbins);
            
#%% Load in data from Terremoto zip

data_r = np.load('run_data_r_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int

# Assign position data
r_int = data_r['a'];
rc_int = data_r['b'];
ra_int = data_r['c'];

# Concatenate data for functions
rcb_int = np.concatenate((rc_int, r_int));
rca_int = np.concatenate((rc_int, ra_int));
rab_int = np.concatenate((ra_int, r_int));
rt_int = np.concatenate((r_int, rca_int));

#%% Create class instances

### KEY for g(r):

# g0(r) is for polymer beads
# g1(r) is for polymer beads within the same chain (intra)
# g2(r) is for polymer beads not on the same chain (inter)
# g3(r) is for cations
# g4(r) is for anions
# g5(r) is for cation/bead pairs
# g6(r) is for cation/anion pairs
# g7(r) is for anion/bead pairs
# g8(r) is for the whole system

### KEY for CN(r):
# g5(r) is for cation/bead pairs
# g6(r) is for cation/anion pairs
# g7(r) is for anion/bead pairs

# Create instances

g_r = [];

g_r.append(RDF(N*M, L, nbins, r_int, rhop,  n_t));
g_r.append(RDF(N*M, L, nbins, r_int, rhop,  n_t));
g_r.append(RDF(N*M, L, nbins, r_int, rhop,  n_t));
g_r.append(RDF(I, L, nbins, rc_int, rhoc, n_t));
g_r.append(RDF(I, L, nbins, ra_int, rhoa, n_t));
g_r.append(RDF(I, L, nbins, rcb_int, rhop, n_t));
g_r.append(RDF(I, L, nbins, rca_int, rhoa, n_t));
g_r.append(RDF(I, L, nbins, rab_int, rhop, n_t));
g_r.append(RDF(N*M+2*I, L, nbins, rt_int, rho, n_t));

# Initialize g(r) array
g = np.zeros((nbins,len(g_r)));

# Get g(r) for each class instance
g[:,0] = g_r[0].RDF1();
g[:,1] = g_r[1].RDF2(M);
g[:,2] = g_r[2].RDF3(M);
g[:,3] = g_r[3].RDF1();
g[:,4] = g_r[4].RDF1();
g[:,5] = g_r[5].RDF4(N*M);
g[:,6] = g_r[6].RDF4(I);
g[:,7] = g_r[7].RDF4(N*M);
g[:,8] = g_r[8].RDF1();

#%% Get coordination number functions with g(r)

# Initialize array for coordination number
fCN = np.zeros((nbins,3));

# Compute coordination number for pair distributions
fCN[:,0] = fcn_slow(r_RDF, g[:,5], rhop, nbins, L);
fCN[:,1] = fcn_slow(r_RDF, g[:,6], rhoa, nbins, L);
fCN[:,2] = fcn_slow(r_RDF, g[:,7], rhop, nbins, L);


#%% Save to npz

np.savez_compressed('RDF_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=g, b=fCN);