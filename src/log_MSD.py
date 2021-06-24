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
from parameters import N, M, I, n_tlog, mu_mag, numloop

#%% Load in data from Terremoto zip and calculate MSD for the cations

data_ru_full = np.load('ru_ELEC_LOG.npz');  # a=r_int, b=x_int, c=y_int,d=z_int
ru_full = data_ru_full['a'];
ruc_full = data_ru_full['b'];
rua_full = data_ru_full['c'];

# Initialize arrays
xind = np.arange(0,3*n_tlog,3);
yind = np.arange(1,3*n_tlog,3);
zind = np.arange(2,3*n_tlog,3);

# Polymer
xu_full = ru_full[:,xind];
yu_full = ru_full[:,yind];
zu_full = ru_full[:,zind];

# Cation
xuc_full = ruc_full[:,xind];
yuc_full = ruc_full[:,yind];
zuc_full = ruc_full[:,zind];

# Anion
xua_full = rua_full[:,xind];
yua_full = rua_full[:,yind];
zua_full = rua_full[:,zind];

#%% Create instances of MSD class

### KEY
# r0 is polymer
# r1 is cations
# r2 is anions

rci = [];

rci.append(MSD(xu_full, yu_full, zu_full, N*M, n_tlog));
rci.append(MSD(xuc_full, yuc_full, zuc_full, I, n_tlog));
rci.append(MSD(xua_full, yua_full, zua_full, I, n_tlog));

# Initialize array
MSD = np.zeros((n_tlog, 6));

# Perform calculation
MSD[:, 0] =  rci[0].MSD2(numloop);
MSD[:, 1] = rci[1].MSD2(numloop);
MSD[:, 2] = rci[2].MSD2(numloop);
MSD[:, 3] = rci[0].MSD3(numloop);
MSD[:, 4] = rci[1].MSD3(numloop);
MSD[:, 5] = rci[2].MSD3(numloop);

#%% Save data in compressed file

np.savez_compressed('MSD_log_{}_{}_{}_{:.3f}'.format(N, M, I, mu_mag), a=MSD);