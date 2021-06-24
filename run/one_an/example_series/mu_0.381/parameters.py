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
m = np.zeros((2));
m[0] = 47.5500;
m[1] = 144.9642;

# Diameters (Angstrom)
d = np.zeros(2);
N_Av = 6.0221409*(10**23);
d[0] = 3.9000;
d[1] = 5.8000;

# LJ parameters;
oij = np.zeros(3);
oij[0] = 3.9000; #11
oij[1] = 4.8500; #13
oij[2] = 5.8000; #22


# Cutoffs (Angstrom)
rc = np.zeros(3);
rc[0] = 9.7500;
rc[1] = 12.1250;
rc[2] = 14.5000;

# Densities (g cm-3)
p = np.zeros(2);
p[0] = 2.5422;
p[1] = 2.3563;

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
I = 1; # Number of ion PAIRS
rho = 1.0608; # g cm3
V = 119315.7414; # Angstrom3
L = 49.2303; # Angstrom
rccoul = 14.0000;

# Charges
e_c = 1.602*10**(-19);
qa = -1.0000; # Multiple of e
mu_mag = 0.3815; # q Angstrom


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
