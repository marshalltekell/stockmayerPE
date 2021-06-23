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
highT = !!!highT; #K

# Masses (g mol-1)
m = np.zeros((3));
m[0] = !!!m0;
m[1] = !!!m1;
m[2] = !!!m2;

# Diameters (Angstrom)
d = np.zeros(3);
N_Av = 6.0221409*(10**23);
d[0] = !!!d0;
d[1] = !!!d1;
d[2] = !!!d2;

# LJ parameters;
oij = np.zeros(6);
oij[0] = !!!r11; #11
oij[1] = !!!r12; #12
oij[2] = !!!r13; #13
oij[3] = !!!r22; #22
oij[4] = !!!r23; #23
oij[5] = !!!r33; #33


# Cutoffs (Angstrom)
rc = np.zeros(6);
rc[0] = !!!rc11;
rc[1] = !!!rc12;
rc[2] = !!!rc13;
rc[3] = !!!rc22;
rc[4] = !!!rc23;
rc[5] = !!!rc33;

# Densities (g cm-3)
p = np.zeros(3);
p[0] = !!!p0;
p[1] = !!!p1;
p[2] = !!!p2;

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
I = !!!I; # Number of ion PAIRS
rho = !!!rho; # g cm3
V = !!!V; # Angstrom3
L = !!!L; # Angstrom
rccoul = !!!rccoul;

# Charges
e_c = 1.602*10**(-19);
qc = !!!qc; # Multiple of e
qa = !!!qa # Multiple of e
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
