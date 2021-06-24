#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:46:45 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package
from classes import MSD

# Import parameters for simulation
from parameters import N, M, I, n_t, mu_mag
        
#%% Load in data

data_ru = np.load('run_data_ru_lammps_ELEC.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
ru_int = data_ru['a'];
ruc_int = data_ru['b'];
rua_int = data_ru['c'];

# Initialize arrays
xind = np.arange(0,3*n_t,3);
yind = np.arange(1,3*n_t,3);
zind = np.arange(2,3*n_t,3);

# Polymer
xu_int = ru_int[:,xind];
yu_int = ru_int[:,yind];
zu_int = ru_int[:,zind];

# Cation
xuc_int = ruc_int[:,xind];
yuc_int = ruc_int[:,yind];
zuc_int = ruc_int[:,zind];

# Anion
xua_int = rua_int[:,xind];
yua_int = rua_int[:,yind];
zua_int = rua_int[:,zind];

#%% Create class instances

### KEY
# r0 is polymer
# r1 is cations
# r2 is anions


rci = [];

rci.append(MSD(xu_int, yu_int, zu_int, N*M, n_t));
rci.append(MSD(xuc_int, yuc_int, zuc_int, I, n_t));
rci.append(MSD(xua_int, yua_int, zua_int, I, n_t));

# Initialize array
MSD = np.zeros((n_t, 6));

# Perform calculation
MSD[:, 0] =  rci[0].MSD0();
MSD[:, 1] = rci[1].MSD0();
MSD[:, 2] = rci[2].MSD0();
MSD[:, 3] = rci[0].MSD1();
MSD[:, 4] = rci[1].MSD1();
MSD[:, 5] = rci[2].MSD1();


#%% Save data in compressed file

np.savez_compressed('MSD_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=MSD);
