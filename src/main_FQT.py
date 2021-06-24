#!/usr/bin/env python3 -u
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:06:48 2020

@author: marshalltekell

NVE MD with only Lennard-Jones interactions
"""
# Equation of state of the Lennard-Jones fluid

# Import packages
import numpy as np #python array manipulation package
from mpi4py import MPI

# Import user-defined functions from separate scripts
import cfunctions
from analysis_functions import gen_q

# Import simulation parameters
from parameters import n_t, N, M, I;
from param_NVT import L;

#%% Setup MPI conditions

comm = MPI.COMM_WORLD;
rank = comm.Get_rank();
proc = comm.Get_size();

#%% Additional parameters necessary for calculating F(q,t)

# Time array
#time = np.arange(0,N_steps,p);

# Generate q magnitudes
dk = (2*np.pi)/L;
q_mag = [3*dk, 5*dk, 10*dk, 12*dk, 15*dk];

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

# Flatten arrays (polymer)
xu_int = xu_int.flatten();
yu_int = yu_int.flatten();
zu_int = zu_int.flatten();

# Flatten arrays (cation)
xuc_int = xuc_int.flatten();
yuc_int = yuc_int.flatten();
zuc_int = zuc_int.flatten();

# Flatten arrays (anion)
xua_int = xua_int.flatten();
yua_int = yua_int.flatten();
zua_int = zua_int.flatten();
        
#%% Calculate the self intermediate scattering function

# Initialize Fqt
if (rank!=0):

    Fqt_s = np.zeros((n_t,3));
            
    # Get dummy arrays
    fqt_0 = np.zeros((n_t));
    fqt_1 = np.zeros((n_t));
    fqt_2 = np.zeros((n_t));
    
    # Generate q vectors at magnitude assigned to processor
    q_vec, num_q = gen_q(q_mag[rank-1], L);
    q_vec = q_vec.flatten();

    # Get normalization count
    ct = np.arange(1,n_t+1,1, dtype=np.intc);
    ct = ct[::-1];
    ct[0] = n_t-1;
    ct = np.ascontiguousarray(ct, dtype=np.intc);

    
    # Update dummy arrays in place
    print("Calculating self intermediate scattering function for polymer at q = {:.3f} on processor.".format(q_mag[rank-1], rank));
    cfunctions.fqt(fqt_0, ct, xu_int, yu_int, zu_int, q_vec, N*M, n_t, num_q);
    print("Done.")
    print("Calculating self intermediate scattering function for cation at q = {:.3f} on processor.".format(q_mag[rank-1], rank));
    cfunctions.fqt(fqt_1, ct, xuc_int, yuc_int, zuc_int, q_vec, I, n_t, num_q);
    print("Done.")
    print("Calculating self intermediate scattering function for anion at q = {:.3f} on processor.".format(q_mag[rank-1], rank));
    cfunctions.fqt(fqt_2, ct, xua_int, yua_int, zua_int, q_vec, I, n_t, num_q);
    print("Done.")
    
    # Give result to Fqt
    Fqt_s[:,0] = fqt_0;
    Fqt_s[:,1] = fqt_1;
    Fqt_s[:,2] = fqt_2;
    
    # Send the total over to the master process
    comm.send(Fqt_s, dest=0);
    
else:
    
    # Initialize array to hold values at different q mags
    Fqt = [];
    
    for i in range(0,3):
        
        Fqt.append(np.zeros((n_t,len(q_mag))));
    
    
    # Receive the messages and put them in the array
    for i in range(1,proc):
              
        Fqt_s = comm.recv(source=i);
        Fqt[0][:,i-1] = Fqt_s[:,0];
        Fqt[1][:,i-1] = Fqt_s[:,1];
        Fqt[2][:,i-1] = Fqt_s[:,2];
    
    # Save to compressed file
    np.savez_compressed('FQT_{}_{}_{}'.format(N, M, I), a=Fqt);
    
        
