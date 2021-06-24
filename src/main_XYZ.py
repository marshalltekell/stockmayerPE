#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:00:02 2020

@author: marshalltekell
"""

#%% Import packages and scripts

# Import packages
import numpy as np #python array manipulation package

# Misc
#from functions import writeion2xyz
from analysis_functions import write2xyz
from analysis_functions import parse_LAMMPs_out_int

# Import from parameters.py
from parameters import N, M, I, n_t;

#%% Parse dump files to collect instantaneous positions, velocities, and forces

data_r = np.load('run_data_r_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int

# Assign position data
r_int = data_r['a'];
rc_int = data_r['b'];
ra_int = data_r['c'];

'''
names_in = ['ID','x','y','z']
temp_df = parse_LAMMPs_out_int("r.out", M*N+2*I, 9, names_in, 1)
r = temp_df[0:M*N,:];
rc = temp_df[M*N:M*N+I,:];
ra = temp_df[M*N+I:M*N+2*I,:];

    
#%% Write to .xyz file for VMD
'''
r =r_int[:, -3:];
rc = rc_int[:, -3:];
ra = ra_int[:, -3:];

write2xyz(N, 1, I, '40_40_100_1p0.xyz', r, rc, ra)