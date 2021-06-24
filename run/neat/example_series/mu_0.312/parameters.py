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
e = 0.7014; #kcal/mol
o = 3.9000; #Angstrom

# Temperature
T = 353.0000; #K

# Masses (g mol-1)
m = 47.5500;

# Diameters (Angstrom)
d = 3.9000;
N_Av = 6.0221409*(10**23);

# LJ parameters;
oij = 3.9000; #11

# Cutoffs (Angstrom)
rc = 9.7500;

# Densities (g cm-3)
p = 2.5422;

# FENE bond parameters
l_FENE = 3.7590; # Angstrom
bfl = 3.9780;
R0 = 5.8500; # Angstrom
k_FENE = 1.3834;# kcal mol-1 Ang-2

# Cosine potential parameters
kb = 1.38*10**(-23); # J K-1
K = 0.7615; #kcal mol-1 

# System size and contents
M = 40; # Number of chains
N = 40; # Number of beads per chain
rho = 1.0600; # g cm3
V = 119182.8390; # Angstrom3
L = 49.2120; # Angstrom
rccoul = 10.0000;

# Charges
e_c = 1.602*10**(-19);
mu_mag = 0.3121; # q Angstrom


# Production run (with linear sampling) simulation variables
N_steps = 15000000;
p2 = 10000;
n_t = np.int(N_steps/p2);
time = np.arange(0,N_steps,p2);

# Production run (with log sampling) simulation variables
base = 10;
samp = 9;
upmag = 4;
numloop = 100;
n_tlog = upmag*samp + 1;

timelog = np.zeros((n_tlog));
timelog[0] = 0;

for i in range(0,upmag):
    
    for j in range(0,samp):
        
        timelog[samp*i+j + 1] = base**(i)*(j+1);
