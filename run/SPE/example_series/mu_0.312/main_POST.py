#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:34:22 2020

@author: marshalltekell
"""

#%% Import packages and scripts

# Import packages
import numpy as np #python array manipulation package

# Misc
from functions import parse_LAMMPs_out_int
from functions import parse_LAMMPs_out_fix
from functions import update_progress

# Import values from parameters.py
from parameters import N, M, I, n_t;

#%% Parse dump files to collect instantaneous positions, velocities, and forces

# Initialize and allocate arrays for parse
names_in = ['ID','x','y','z']
names_in_fix = ['step','T','U','K','P','MSD']
temp_df = np.zeros((N*M,3));
temp_df_fix = np.zeros((n_t,6));

# Initialize and allocate arrays for per-atom values
r_int = np.zeros((N*M,3*n_t));
rc_int = np.zeros((I,3*n_t));
ra_int = np.zeros((I,3*n_t));
ru_int = np.zeros((N*M,3*n_t));
ruc_int = np.zeros((I,3*n_t));
rua_int = np.zeros((I,3*n_t));
v_int = np.zeros((N*M,3*n_t));
f_int = np.zeros((N*M,3*n_t));
mu_int = np.zeros((N*M,3*n_t));

# Initialize and allocate arrays for global values
T_int = np.zeros((n_t,1));
U_int = np.zeros((n_t,1));
K_int = np.zeros((n_t,1));
P_int = np.zeros((n_t,1));
MSD_av = np.zeros((n_t,1));

print("Parsing output from simulation...")

# Parse global scalars from fix
temp_df_fix = parse_LAMMPs_out_fix("TEPMSD.out",n_t, names_in_fix)
T_int = temp_df_fix[:,1];
U_int = temp_df_fix[:,2];
K_int = temp_df_fix[:,3];
P_int = temp_df_fix[:,4];
MSD_av = temp_df_fix[:,5];

e_tot = U_int + K_int;

# Parse per-atom dumps
for i in range(0,n_t):
    
    # Positions (folded)
    temp_df = parse_LAMMPs_out_int("r.out", M*N+2*I, 9, names_in, i)
    r_int[:,i*3:i*3+3] = temp_df[0:M*N,:];
    rc_int[:,i*3:i*3+3] = temp_df[M*N:M*N+I,:];
    ra_int[:,i*3:i*3+3] = temp_df[M*N+I:M*N+2*I,:];
    
    # Positions (unfolded)
    temp_df = parse_LAMMPs_out_int("ru.out", M*N+2*I, 9, names_in, i)
    ru_int[:,i*3:i*3+3] = temp_df[0:M*N,:];
    ruc_int[:,i*3:i*3+3] = temp_df[M*N:M*N+I,:];
    rua_int[:,i*3:i*3+3] = temp_df[M*N+I:M*N+2*I,:];
    
    '''
    # Velocities
    temp_df = parse_LAMMPs_out_int("v.out", M*N+2*I, 9, names_in, i)
    v_int[:,i*3:i*3+3] = temp_df[0:M*N,:];
    
    # Forces
    temp_df = parse_LAMMPs_out_int("f.out", M*N+2*I, 9, names_in, i)
    f_int[:,i*3:i*3+3] = temp_df[0:M*N,:];
    '''
    # Dipoles
    temp_df = parse_LAMMPs_out_int("mu.out", M*N+2*I, 9, names_in, i)
    mu_int[:,i*3:i*3+3] = temp_df[0:M*N,:];
    
    update_progress(i/n_t)


#%% Save all of the data to a compressed file

np.savez_compressed('run_data_r_lammps_ELEC', a=r_int, b=rc_int, c=ra_int);
np.savez_compressed('run_data_ru_lammps_ELEC', a=ru_int, b=ruc_int, c=rua_int);
#np.savez_compressed('run_data_v_lammps_ELEC', a=v_int);
#np.savez_compressed('run_data_f_lammps_ELEC', a=f_int);
np.savez_compressed('run_data_e_lammps_ELEC', a=U_int, b=K_int, c=e_tot);
np.savez_compressed('run_data_TPMSD_lammps_ELEC', a=T_int, b=P_int, c=MSD_av);
np.savez_compressed('run_data_mu_lammps_ELEC', a=mu_int);
