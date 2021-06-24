#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 10:41:40 2020

@author: marshalltekell
"""
#%% Import packages

import numpy as np
import pandas as pd
import os
from sh import sed

# Import parameters
from parameters import m, I, N, M, N_Av;

#%% Find average volume from NPT run

pathname = 'equil_vol.out';
col_name = ['step','V'];
nrows = 1000;

# Import data from LAMMPS output and assign to numpy array for analysis
df = pd.read_csv(pathname,sep='\s+', skiprows=2, nrows=1000, skip_blank_lines = True, names = col_name)
df = df.to_numpy();

# Find average volume and calculate density
V = np.mean(df[:,1]);
L = V**(1/3);
rho = ((10**24)/(V*N_Av))*(I*(m[1]+m[2])+N*M*m[0]);


#%% Give this value to in_RUN and in_LOG_RUN

path = os.getcwd()
run_path = os.path.join(path, 'in_RUN.ELEC')
param_NVT_path = os.path.join(path, 'param_NVT.py')
sed(['-i', 's/!!!L2/{:.4f}/'.format(L), run_path])
sed(['-i', '',  's/!!!rho/{:.4f}/'.format(rho), param_NVT_path])
sed(['-i', '',  's/!!!V/{:.4f}/'.format(V), param_NVT_path])
sed(['-i', '',  's/!!!L/{:.4f}/'.format(L), param_NVT_path])