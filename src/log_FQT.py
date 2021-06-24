#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 21:11:45 2020

@author: marshalltekell

Calculate the self intermediate scattering function for log time sampling

"""

# Import packages and scripts
import numpy as np
import cfunctions
from analysis_functions import gen_q
from analysis_functions import update_progress

# Import parameters
from parameters import N, M, I, n_tlog, numloop;
from param_NVT import L;

#%% Read in parameters that specify the conditions of the the run

# Generate q array
dk = (2*np.pi)/L;
q_mag = [3*dk, 5*dk, 10*dk, 12*dk, 15*dk];
        
#%% Unpack data from compresssed file

data_ru_log = np.load('ru_ELEC_LOG.npz');  # a=ru_int, b=x_int, c=y_int,d=z_int
ru_full = data_ru_log['a'];
ruc_full = data_ru_log['b'];
rua_full = data_ru_log['c'];

#%% Calculate the self intermediate scattering function for the series of q values

# Initialize overall arrays
Fqt = [];

for i in range(0,3):
    
    t1 = np.zeros((n_tlog,len(q_mag)));
    Fqt.append(t1);

# Indices
xind = np.arange(0,3*n_tlog,3);
yind = np.arange(1,3*n_tlog,3);
zind = np.arange(2,3*n_tlog,3);

# Loop over q magnitudes
for i in range(0, len(q_mag)):
            
    # Generate q vectors
    q_vec, num_q = gen_q(q_mag[i], L);
    q_vec = q_vec.flatten();
            
    # Loop over the number of run loops
    for j in range(0,numloop):
        
        # Get fqt vector to alter in place
        fqt_0 = np.zeros((n_tlog));
        fqt_1 = np.zeros((n_tlog));
        fqt_2 = np.zeros((n_tlog));
                    
        # Get the data from that run
        ru_int = ru_full[j*N*M:(j+1)*N*M,:];
        ruc_int = ruc_full[j*I:(j+1)*I,:];
        rua_int = rua_full[j*I:(j+1)*I,:];
        
        # Split into x y and z
        xu_int = ru_int[:,xind];
        yu_int = ru_int[:,yind];
        zu_int = ru_int[:,zind];
        
        # Split into x y and z
        xuc_int = ruc_int[:,xind];
        yuc_int = ruc_int[:,yind];
        zuc_int = ruc_int[:,zind];
        
        # Split into x y and z
        xua_int = rua_int[:,xind];
        yua_int = rua_int[:,yind];
        zua_int = rua_int[:,zind];
        
        # Cast as 1D arrays (bead)
        xu_int = xu_int.flatten();
        yu_int = yu_int.flatten();
        zu_int = zu_int.flatten();
        
        # Cast as 1D arrays (cation)
        xuc_int = xuc_int.flatten();
        yuc_int = yuc_int.flatten();
        zuc_int = zuc_int.flatten();
        
        # Cast as 1D arrays (anion)
        xua_int = xua_int.flatten();
        yua_int = yua_int.flatten();
        zua_int = zua_int.flatten();
        
        # Update fqt in place
        cfunctions.logfqt(fqt_0, xu_int, yu_int, zu_int, q_vec, N*M, n_tlog, num_q);
        cfunctions.logfqt(fqt_1, xuc_int, yuc_int, zuc_int, q_vec, I, n_tlog, num_q);
        cfunctions.logfqt(fqt_2, xua_int, yua_int, zua_int, q_vec, I, n_tlog, num_q);
        
        # Add results to total
        Fqt[0][:,i] += fqt_0;
        Fqt[1][:,i] += fqt_1;
        Fqt[2][:,i] += fqt_2;
        
    # Divide the number of loops
    Fqt[0][:,i] = Fqt[0][:,i]/numloop;
    Fqt[1][:,i] = Fqt[1][:,i]/numloop;
    Fqt[2][:,i] = Fqt[2][:,i]/numloop;
    
    # Update progress
    update_progress((i+1)/len(q_mag))
    
    

#%% Save to compresssed file
    
np.savez_compressed('FQT_log_{}_{}_{}'.format(N, M, I), a=Fqt);
