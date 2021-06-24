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
from functions import write2lammpsin_onean_angle_sphere

# Import values from parameters
from parameters import e, o, m, d, p, T, N, M, l_FENE, bfl;
from parameters import L, mu_mag;

#%% Initialize the system

# Create randomly distributed dipoles for each of the polymer beads
mu = np.zeros((N*M,3));
mu = gen_mu(N, M, mu_mag);

# Generate chains and place them in a simulation box much larger than the final size
r = np.zeros((M*N,3));
ru_init = np.zeros((M*N,3));

# Use gen_chainbox_lattice
r = gen_chainbox_lattice(M, r, ru_init, l_FENE, bfl, N, L, o, e, T);

# Place ion in center of box
ra = np.array([L/2, L/2, L/2], dtype=np.float);

#%% Write initial atom positions to in file for lammps

write2lammpsin_onean_angle_sphere(N, M, L, "ONEAN_PRE.in", m, p, d, r, ra, mu)

#%% Run lammps serial with soft pairwise potential to move off overlaps

lmp = lammps();
lmp.file("in_PRE.ONEAN");

data = lmp.gather_atoms("x",1,3);
r = np.array(data).reshape(-1,3);
ra = r[M*N:M*N+1,:];
r = r[0:M*N,:];
ra = np.array([ra[0,0], ra[0,1],ra[0,2]])

#%% Write atom positions to in file for lammps after NPT

write2lammpsin_onean_angle_sphere(N, M, L, "ONEAN_EQUIL.in", m, p, d, r, ra, mu)
