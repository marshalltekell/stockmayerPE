#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:32:50 2020

@author: marshalltekell

Parameters defined from simulation master script, used in analysis files
"""

import numpy as np

#%% Set global variables

# LJ force-field parameters
e = !!!e; #kcal/mol
o = !!!o; #Angstrom

# Temperature
T = !!!T; #K

# Masses (g mol-1)
m = !!!m0;

# Diameters (Angstrom)
d = !!!d0;
N_Av = 6.0221409*(10**23);

# LJ parameters;
oij = !!!r11; #11

# Cutoffs (Angstrom)
rc = !!!rc11;

# Densities (g cm-3)
p = !!!p0;

# FENE bond parameters
l_FENE = !!!l_FENE; # Angstrom
bfl = !!!bfl;
R0 = !!!R0; # Angstrom
k_FENE = !!!k_FENE;# kcal mol-1 Ang-2

# Cosine potential parameters
kb = 1.38*10**(-23); # J K-1
K = !!!K; #kcal mol-1 

# System size and contents
M = !!!M; # Number of chains
N = !!!N; # Number of beads per chain
rho = !!!rho; # g cm3
V = !!!V; # Angstrom3
L = !!!L; # Angstrom
rccoul = !!!rccoul;

# Charges
e_c = 1.602*10**(-19);
mu_mag = !!!mu_mag; # q Angstrom


# Production run (with linear sampling) simulation variables
N_steps = !!!prod_run;
p2 = !!!prod_sf;
n_t = np.int(N_steps/p2);
time = np.arange(0,N_steps,p2);

# Production run (with log sampling) simulation variables
base = !!!base;
samp = !!!samp;
upmag = !!!upmag;
numloop = !!!loops;
n_tlog = upmag*samp + 1;

timelog = np.zeros((n_tlog));
timelog[0] = 0;

for i in range(0,upmag):
    
    for j in range(0,samp):
        
        timelog[samp*i+j + 1] = base**(i)*(j+1);
