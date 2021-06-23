#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:34:22 2020

@author: marshalltekell
"""

#%% Import packages and scripts

# Import packagesvi 
from lammps import lammps
import numpy as np #python array manipulation package

# Set-up
from functions import gen_chainbox_lattice
from functions import gen_mu
from functions import lattice_positions
from functions import write2lammpsin_elec_angle
from functions import write2lammpsin_elec_angle_sphere

# Import values from parameters
from parameters import e, o, m, d, p, T, N, M, I, l_FENE, bfl;
from parameters import L, qc, qa, mu_mag;

#%% Initialize the system

# Create randomly distributed dipoles for each of the polymer beads
mu = np.zeros((N*M,3));
mu = gen_mu(N, M, mu_mag);

# Generate chains and place them in a simulation box much larger than the final size
r = np.zeros((M*N,3));
ru_init = np.zeros((M*N,3));

#%%

# Use gen_chainbox_lattice
r = gen_chainbox_lattice(M, r, ru_init, l_FENE, bfl, N, L, o, e, T);

#%%

# Place ions on lattice sites in the simulation box
rI = np.zeros((2*I,3));
rI = lattice_positions(2*I,L);

# Divide into cation and anion
rc = np.zeros((I,3));
ra = np.zeros((I,3));
temp = 0;
temp2 = 0;

for i in range(0,2*I):
    
    if (i%2 == 0):
        
        rc[temp,:] = rI[i,:];
        temp += 1;
        
    else:
        
        ra[temp2,:] = rI[i,:];
        temp2 += 1;

#%% Write initial atom positions to in file for lammps

write2lammpsin_elec_angle_sphere(N, M, I, L, "ELEC_PRE.in", m, p, d, r, rc, ra, mu, qc, qa)

#%% Run lammps serial with soft pairwise potential to move off overlaps

lmp = lammps()
lmp.file("in_PRE.ELEC")

data = lmp.gather_atoms("x",1,3)
r = np.array(data).reshape(-1,3)
rc = r[M*N:M*N+I,:]
ra = r[M*N+I:M*N+2*I,:]
r = r[0:M*N,:]

#%% Write atom positions to in file for lammps after NPT

write2lammpsin_elec_angle_sphere(N, M, I, L, "ELEC_EQUIL.in", m, p, d, r, rc, ra, mu, qc, qa)
