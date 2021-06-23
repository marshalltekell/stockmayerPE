#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 19:10:22 2020

@author: marshalltekell
"""
# Import packages
import numpy as np #python array manipulation package
from random import randint
import os 
from shutil import copyfile
from sh import sed

#%% Set global variables

# LJ force-field parameters
e = 0.7014; #kcal/mol
o = 3.9; #Angstrom

# Temperature
T = 353; #K
highT = 1000; #K

# Masses (g mol-1)
m = np.zeros(3);
m[0] = 47.55;
m[1] = 6.941;
m[2] = 144.9642;

# Diameters (Angstrom)
d = np.zeros(3);
N_Av = 6.0221409*(10**23);
d[0] = 3.90;
d[1] = 1.80;
d[2] = 5.80;

# LJ paraemters;
oij = np.zeros(6);
oij[0] = d[0]; #11
#oij[1] = (1/2)*(d[0] + d[1]); #12
oij[1] = 2.5047; #13
oij[2] = (1/2)*(d[0] + d[2]); #13
oij[3] = d[1]; #22
oij[4] = (1/2)*(d[1] + d[2]); #23
oij[5] = d[2]; #33


# Cutoffs (Angstrom)
rc = np.zeros(6);
rc[0] = 2.5*oij[0];
rc[1] = 2.5*oij[1];
rc[2] = 2.5*oij[2];
rc[3] = 2.5*oij[3];
rc[4] = 2.5*oij[4];
rc[5] = 2.5*oij[5];

# Densities (g cm-3)
p = np.zeros(3);
p = (10**24/N_Av)*(m/((4/3)*np.pi*(d/2)**3));

# FENE bond parameters
l_FENE = 0.96384*o; # Angstrom
bfl = 1.02*o; # Angstrom
R0 = 1.5*o; # Angstrom
k_FENE = 30*(e/(o**2));# kcal mol-1 Ang-2

# Cosine potential parameters
kb = 1.38*10**(-23); # J K-1
Tref = 353; # K
K = 1.086*kb*Tref*N_Av*(1/4184); #kcal mol-1 

# System size and contents
M = 40; # Number of chains
N = 40; # Number of beads per chain
I = 100; # Number of ion PAIRS
rho = 1.50*((I*(m[1]+m[2]))/(I*(m[1]+m[2])+N*M*m[0])) + 1.06*((N*M*m[0])/(I*(m[1]+m[2])+N*M*m[0])); # g cm3
rho = rho*1.13;
V = ((10**24)/(rho*N_Av))*(I*(m[2]+m[1]) + N*M*m[0]); # Angstrom3
L = V**(1/3); # Angstrom
rccoul = 17.184;

# Charges
e_c = 1.602*10**(-19);
qc = 1.0;
qa = -1.0;
mu_mag = 1.125*(1.7*3.333564*10**(-30)*10**10)/(e_c); # q Angstrom

#%% Set simulation variables

# Pre (generate system and move-off overlaps)
pre_run = 10000;
pre_tf = 50;
pre_vseed = randint(0, 99999);
Amax = 100;

# Equilibration (relax the system into its lowest energy at constant volume)
equil_vseed = randint(0, 99999);
equil_run_1 = 1000;
equil_run_2 = 100000;
equil_run_3 = 10000000;
equil_tf = 1000;

# Production
prod_run = 15000000
prod_tf = 1000;
prod_sf = 10000;

# Log run
base = 10;
samp = 9;
upmag = 4;
numloop = 100;

#%% Set series variable (the variable that will be altered across the simulations)

SL = 6;
mu_PEO = (1.0*3.333564*10**(-30)*10**10)/(e_c);
mu_mag = np.zeros(SL);

for i in range(SL):
    
    temp =  1+ ((i+3)*(1/SL));
    mu_mag[i] = mu_PEO*temp;

#%% Loop over series variable

#pwd = "/moto/home/mct2180/MD/projects/SPE/series1"

pwd = os.getcwd()

# Open equilibration submission script
in_file_equil = open('sub_equil.sh',"w+")
in_file_equil.write('#!/bin/sh \n')
              
# Open production submission script
in_file = open('sub.sh',"w+")
in_file.write('#!/bin/sh \n')

# Open log production submission script
in_file_log = open('sub_log.sh',"w+")
in_file_log.write('#!/bin/sh \n')


for i in range(0,SL):
    
    # Get random numbers
    seed_equil_1 =  randint(0, 99999);
    seed_equil_2 =  randint(0, 99999);
    seed_equil_3 =  randint(0, 99999);
    seed_run_1 =  randint(0, 99999);
    seed_run_2 =  randint(0, 99999);
    seed_run_3 =  randint(0, 99999);
    seed_log_1 =  randint(0, 99999);
    seed_log_2 =  randint(0, 99999);
    seed_log_3 =  randint(0, 99999);
    
    # Create directory
    dicy = 'mu_{:.3f}'.format(mu_mag[i]);
    path = os.path.join(pwd, dicy);
    os.mkdir(path);

    # Copy functions
    tpath = os.path.join(path, 'functions.py')
    copyfile('functions.py',tpath)

    # Copy .py templates
    
    # Pre
    pre_path = os.path.join(path, 'main_PRE.py')
    copyfile('templates/main_PRE.py', pre_path)
    
    # Post
    post_path = os.path.join(path, 'main_POST.py')
    copyfile('templates/main_POST.py',post_path)
    
    # Post log
    postlog_path = os.path.join(path, 'main_POST_LOG.py')
    copyfile('templates/main_POST_LOG.py',postlog_path)

    # in_PRE
    in_pre_path = os.path.join(path, 'in_PRE.ELEC')
    copyfile('templates/in_PRE.ELEC',in_pre_path)
    
    # in_EQUIL
    in_equil_path = os.path.join(path, 'in_EQUIL.ELEC')
    copyfile('templates/in_EQUIL.ELEC',in_equil_path)
    
    # in_RUN
    in_run_path = os.path.join(path, 'in_RUN.ELEC')
    copyfile('templates/in_RUN.ELEC',in_run_path)
    
    # in_LOG_RUN
    in_log_path = os.path.join(path, 'in_LOG_RUN.ELEC')
    copyfile('templates/in_LOG_RUN.ELEC',in_log_path)
    
    # run_equil.sh
    run_equil_path = os.path.join(path, 'run_equil.sh')
    copyfile('templates/run_equil.sh',run_equil_path)
    
    # run.sh
    run_path = os.path.join(path, 'run.sh')
    copyfile('templates/run.sh',run_path)
    
    # run_log.sh
    run_log_path = os.path.join(path, 'run_log.sh')
    copyfile('templates/run_log.sh',run_log_path)
    
    # parameters
    param_path = os.path.join(path, 'parameters.py')
    copyfile('templates/parameters.py',param_path)
    
    # npt2nvt
    npt2nvt_path =  os.path.join(path, 'npt2nvt.py')
    copyfile('templates/npt2nvt.py', npt2nvt_path)
    
    # param_NVT
    param_NVT_path =  os.path.join(path, 'param_NVT.py')
    copyfile('templates/param_NVT.py', param_NVT_path)
    
    # analysis.tar
    analysis_path =  os.path.join(path, 'analysis.tar')
    copyfile('templates/analysis.tar', analysis_path)
    

    # Sed in_PRE, in_EQUIL, in_RUN, in_LOG_RUN, main_POST, main_POST_LOG, run.sh,
    
    # Sed in_PRE
    sed(['-i', '',  's/!!!k_FENE/{:.4f}/'.format(k_FENE), in_pre_path])
    sed(['-i', '',  's/!!!R0/{:.4f}/'.format(R0), in_pre_path])
    sed(['-i', '',  's/!!!e/{:.4f}/'.format(e), in_pre_path])
    sed(['-i', '',  's/!!!o/{:.4f}/'.format(o), in_pre_path])
    sed(['-i', '',  's/!!!T/{:.4f}/'.format(T), in_pre_path])
    sed(['-i', '',  's/!!!pre_vseed/{:d}/'.format(pre_vseed), in_pre_path])
    sed(['-i', '',  's/!!!r0/{:.4f}/'.format(d[0]), in_pre_path])
    sed(['-i', '',  's/!!!r1/{:.4f}/'.format(d[1]), in_pre_path])
    sed(['-i', '',  's/!!!r2/{:.4f}/'.format(d[2]), in_pre_path])
    sed(['-i', '',  's/!!!Amax/{:d}/'.format(Amax), in_pre_path])
    sed(['-i', '',  's/!!!pre_tf/{:d}/'.format(pre_tf), in_pre_path])
    sed(['-i', '',  's/!!!pre_run/{:d}/'.format(pre_run), in_pre_path])
    
    # Sed in_EQUIL
    sed(['-i', '',  's/!!!k_FENE/{:.4f}/'.format(k_FENE), in_equil_path])
    sed(['-i', '',  's/!!!R0/{:.4f}/'.format(R0), in_equil_path])
    sed(['-i', '',  's/!!!e/{:.4f}/'.format(e), in_equil_path])
    sed(['-i', '',  's/!!!o/{:.4f}/'.format(o), in_equil_path])
    sed(['-i', '',  's/!!!T/{:.4f}/'.format(T), in_equil_path])
    sed(['-i', '',  's/!!!highT/{:.4f}/'.format(highT), in_equil_path])
    sed(['-i', '',  's/!!!quil_vseed/{:d}/'.format(equil_vseed), in_equil_path])
    sed(['-i', '',  's/!!!r11/{:.4f}/'.format(oij[0]), in_equil_path])
    sed(['-i', '',  's/!!!r12/{:.4f}/'.format(oij[1]), in_equil_path])
    sed(['-i', '',  's/!!!r13/{:.4f}/'.format(oij[2]), in_equil_path])
    sed(['-i', '',  's/!!!r22/{:.4f}/'.format(oij[3]), in_equil_path])
    sed(['-i', '',  's/!!!r23/{:.4f}/'.format(oij[4]), in_equil_path])
    sed(['-i', '',  's/!!!r33/{:.4f}/'.format(oij[5]), in_equil_path])
    sed(['-i', '',  's/!!!rc11/{:.4f}/'.format(rc[0]), in_equil_path])
    sed(['-i', '',  's/!!!rc12/{:.4f}/'.format(rc[1]), in_equil_path])
    sed(['-i', '',  's/!!!rc13/{:.4f}/'.format(rc[2]), in_equil_path])
    sed(['-i', '',  's/!!!rc22/{:.4f}/'.format(rc[3]), in_equil_path])
    sed(['-i', '',  's/!!!rc23/{:.4f}/'.format(rc[4]), in_equil_path])
    sed(['-i', '',  's/!!!rc33/{:.4f}/'.format(rc[5]), in_equil_path])
    sed(['-i', '',  's/!!!rccoul/{:.4f}/'.format(rccoul), in_equil_path])
    sed(['-i', '',  's/!!!K/{:.4f}/'.format(K), in_equil_path])
    sed(['-i', '',  's/!!!quil_tf/{:d}/'.format(equil_tf), in_equil_path])
    sed(['-i', '',  's/!!!quil_run_1/{:d}/'.format(equil_run_1), in_equil_path])
    sed(['-i', '',  's/!!!quil_run_2/{:d}/'.format(equil_run_2), in_equil_path])
    sed(['-i', '',  's/!!!quil_run_3/{:d}/'.format(equil_run_3), in_equil_path])
    sed(['-i', '',  's/!!!seed_equil_1/{:d}/'.format(seed_equil_1), in_equil_path])
    sed(['-i', '',  's/!!!seed_equil_2/{:d}/'.format(seed_equil_2), in_equil_path])
    sed(['-i', '',  's/!!!seed_equil_3/{:d}/'.format(seed_equil_3), in_equil_path])
    
    # Sed in_RUN
    sed(['-i', '',  's/!!!k_FENE/{:.4f}/'.format(k_FENE), in_run_path])
    sed(['-i', '',  's/!!!R0/{:.4f}/'.format(R0), in_run_path])
    sed(['-i', '',  's/!!!e/{:.4f}/'.format(e), in_run_path])
    sed(['-i', '',  's/!!!o/{:.4f}/'.format(o), in_run_path])
    sed(['-i', '',  's/!!!T/{:.4f}/'.format(T), in_run_path])
    sed(['-i', '',  's/!!!r11/{:.4f}/'.format(oij[0]), in_run_path])
    sed(['-i', '',  's/!!!r12/{:.4f}/'.format(oij[1]), in_run_path])
    sed(['-i', '',  's/!!!r13/{:.4f}/'.format(oij[2]), in_run_path])
    sed(['-i', '',  's/!!!r22/{:.4f}/'.format(oij[3]), in_run_path])
    sed(['-i', '',  's/!!!r23/{:.4f}/'.format(oij[4]), in_run_path])
    sed(['-i', '',  's/!!!r33/{:.4f}/'.format(oij[5]), in_run_path])
    sed(['-i', '',  's/!!!rc11/{:.4f}/'.format(rc[0]), in_run_path])
    sed(['-i', '',  's/!!!rc12/{:.4f}/'.format(rc[1]), in_run_path])
    sed(['-i', '',  's/!!!rc13/{:.4f}/'.format(rc[2]), in_run_path])
    sed(['-i', '',  's/!!!rc22/{:.4f}/'.format(rc[3]), in_run_path])
    sed(['-i', '',  's/!!!rc23/{:.4f}/'.format(rc[4]), in_run_path])
    sed(['-i', '',  's/!!!rc33/{:.4f}/'.format(rc[5]), in_run_path])
    sed(['-i', '',  's/!!!rccoul/{:.4f}/'.format(rccoul), in_run_path])
    sed(['-i', '',  's/!!!K/{:.4f}/'.format(K), in_run_path])
    sed(['-i', '',  's/!!!prod_tf/{:d}/'.format(prod_tf), in_run_path])
    sed(['-i', '',  's/!!!prod_sf/{:d}/'.format(prod_sf), in_run_path])
    sed(['-i', '',  's/!!!prod_run/{:d}/'.format(prod_run), in_run_path])
    sed(['-i', '',  's/!!!seed_run_1/{:d}/'.format(seed_run_1), in_run_path])
    sed(['-i', '',  's/!!!seed_run_2/{:d}/'.format(seed_run_2), in_run_path])
    sed(['-i', '',  's/!!!seed_run_3/{:d}/'.format(seed_run_3), in_run_path])
    
    # Sed in_LOG_RUN
    sed(['-i', '',  's/!!!k_FENE/{:.4f}/'.format(k_FENE), in_log_path])
    sed(['-i', '',  's/!!!R0/{:.4f}/'.format(R0), in_log_path])
    sed(['-i', '',  's/!!!e/{:.4f}/'.format(e), in_log_path])
    sed(['-i', '',  's/!!!o/{:.4f}/'.format(o), in_log_path])
    sed(['-i', '',  's/!!!T/{:.4f}/'.format(T), in_log_path])
    sed(['-i', '',  's/!!!r11/{:.4f}/'.format(oij[0]), in_log_path])
    sed(['-i', '',  's/!!!r12/{:.4f}/'.format(oij[1]), in_log_path])
    sed(['-i', '',  's/!!!r13/{:.4f}/'.format(oij[2]), in_log_path])
    sed(['-i', '',  's/!!!r22/{:.4f}/'.format(oij[3]), in_log_path])
    sed(['-i', '',  's/!!!r23/{:.4f}/'.format(oij[4]), in_log_path])
    sed(['-i', '',  's/!!!r33/{:.4f}/'.format(oij[5]), in_log_path])
    sed(['-i', '',  's/!!!rc11/{:.4f}/'.format(rc[0]), in_log_path])
    sed(['-i', '',  's/!!!rc12/{:.4f}/'.format(rc[1]), in_log_path])
    sed(['-i', '',  's/!!!rc13/{:.4f}/'.format(rc[2]), in_log_path])
    sed(['-i', '',  's/!!!rc22/{:.4f}/'.format(rc[3]), in_log_path])
    sed(['-i', '',  's/!!!rc23/{:.4f}/'.format(rc[4]), in_log_path])
    sed(['-i', '',  's/!!!rc33/{:.4f}/'.format(rc[5]), in_log_path])
    sed(['-i', '',  's/!!!rccoul/{:.4f}/'.format(rccoul), in_log_path])
    sed(['-i', '',  's/!!!K/{:.4f}/'.format(K), in_log_path])
    sed(['-i', '',  's/!!!prod_tf/{:d}/'.format(prod_tf), in_log_path])
    sed(['-i', '',  's/!!!prod_sf/{:d}/'.format(prod_sf), in_log_path])
    sed(['-i', '',  's/!!!prod_run/{:d}/'.format(prod_run), in_log_path])
    sed(['-i', '',  's/!!!seed_log_1/{:d}/'.format(seed_log_1), in_log_path])
    sed(['-i', '',  's/!!!seed_log_2/{:d}/'.format(seed_log_2), in_log_path])
    sed(['-i', '',  's/!!!seed_log_3/{:d}/'.format(seed_log_3), in_log_path])
        
    # Sed run_equil.sh
    jobname = 'mu_{:.3f}'.format(mu_mag[i]);
    sed(['-i', '',  's/!!!jobname/{:s}/'.format(jobname), run_equil_path])
    
    # Sed run.sh
    jobname = 'mu_{:.3f}'.format(mu_mag[i]);
    sed(['-i', '',  's/!!!jobname/{:s}/'.format(jobname), run_path])
    
    # Sed run_log.sh
    jobname = 'mu_{:.3f}'.format(mu_mag[i]);
    sed(['-i', '',  's/!!!jobname/{:s}/'.format(jobname), run_log_path])
    
    # parameters
    sed(['-i', '',  's/!!!e/{:.4f}/'.format(e), param_path])
    sed(['-i', '',  's/!!!o/{:.4f}/'.format(o), param_path])
    sed(['-i', '',  's/!!!T/{:.4f}/'.format(T), param_path])
    sed(['-i', '',  's/!!!highT/{:.4f}/'.format(highT), param_path])
    sed(['-i', '',  's/!!!m0/{:.4f}/'.format(m[0]), param_path])
    sed(['-i', '',  's/!!!m1/{:.4f}/'.format(m[1]), param_path])
    sed(['-i', '',  's/!!!m2/{:.4f}/'.format(m[2]), param_path])
    sed(['-i', '',  's/!!!d0/{:.4f}/'.format(d[0]), param_path])
    sed(['-i', '',  's/!!!d1/{:.4f}/'.format(d[1]), param_path])
    sed(['-i', '',  's/!!!d2/{:.4f}/'.format(d[2]), param_path])
    sed(['-i', '',  's/!!!r11/{:.4f}/'.format(oij[0]), param_path])
    sed(['-i', '',  's/!!!r12/{:.4f}/'.format(oij[1]), param_path])
    sed(['-i', '',  's/!!!r13/{:.4f}/'.format(oij[2]), param_path])
    sed(['-i', '',  's/!!!r22/{:.4f}/'.format(oij[3]), param_path])
    sed(['-i', '',  's/!!!r23/{:.4f}/'.format(oij[4]), param_path])
    sed(['-i', '',  's/!!!r33/{:.4f}/'.format(oij[5]), param_path])
    sed(['-i', '',  's/!!!rc11/{:.4f}/'.format(rc[0]), param_path])
    sed(['-i', '',  's/!!!rc12/{:.4f}/'.format(rc[1]), param_path])
    sed(['-i', '',  's/!!!rc13/{:.4f}/'.format(rc[2]), param_path])
    sed(['-i', '',  's/!!!rc22/{:.4f}/'.format(rc[3]), param_path])
    sed(['-i', '',  's/!!!rc23/{:.4f}/'.format(rc[4]), param_path])
    sed(['-i', '',  's/!!!rc33/{:.4f}/'.format(rc[5]), param_path])
    sed(['-i', '',  's/!!!p0/{:.4f}/'.format(p[0]), param_path])
    sed(['-i', '',  's/!!!p1/{:.4f}/'.format(p[1]), param_path])
    sed(['-i', '',  's/!!!p2/{:.4f}/'.format(p[2]), param_path])
    sed(['-i', '',  's/!!!l_FENE/{:.4f}/'.format(l_FENE), param_path])
    sed(['-i', '',  's/!!!bfl/{:.4f}/'.format(bfl), param_path])
    sed(['-i', '',  's/!!!R0/{:.4f}/'.format(R0), param_path])
    sed(['-i', '',  's/!!!k_FENE/{:.4f}/'.format(k_FENE), param_path])
    sed(['-i', '',  's/!!!K/{:.4f}/'.format(K), param_path])
    sed(['-i', '',  's/!!!N/{:d}/'.format(N), param_path])
    sed(['-i', '',  's/!!!M/{:d}/'.format(M), param_path])
    sed(['-i', '',  's/!!!I/{:d}/'.format(I), param_path])
    sed(['-i', '',  's/!!!rho/{:.4f}/'.format(rho), param_path])
    sed(['-i', '',  's/!!!V/{:.4f}/'.format(V), param_path])
    sed(['-i', '',  's/!!!L/{:.4f}/'.format(L), param_path])
    sed(['-i', '',  's/!!!rccoul/{:.4f}/'.format(rccoul), param_path])
    sed(['-i', '',  's/!!!qc/{:.4f}/'.format(qc), param_path])
    sed(['-i', '',  's/!!!qa/{:.4f}/'.format(qa), param_path])
    sed(['-i', '',  's/!!!mu_mag/{:.4f}/'.format(mu_mag[i]), param_path])
    sed(['-i', '',  's/!!!prod_run/{:d}/'.format(prod_run), param_path])
    sed(['-i', '',  's/!!!prod_sf/{:d}/'.format(prod_sf), param_path])
    sed(['-i', '',  's/!!!base/{:d}/'.format(base), param_path])
    sed(['-i', '',  's/!!!samp/{:d}/'.format(samp), param_path])
    sed(['-i', '',  's/!!!upmag/{:d}/'.format(upmag), param_path])
    sed(['-i', '',  's/!!!loops/{:d}/'.format(numloop), param_path])

    # Write sub_equil.sh
    in_file_equil.write('cd {:s}\n'.format(dicy))
    in_file_equil.write('sbatch run_equil.sh\n')
    in_file_equil.write('cd ..\n')
    
    # Write sub.sh
    in_file.write('cd {:s}\n'.format(dicy))
    in_file.write('sbatch run.sh\n')
    in_file.write('cd ..\n')
    
    # Write sub_log.sh
    in_file_log.write('cd {:s}\n'.format(dicy))
    in_file_log.write('sbatch run_log.sh\n')
    in_file_log.write('cd ..\n')

# (Don't know why I have to do this) sub_equil.sh    
in_file_equil = open('sub_equil.sh',"w+")
in_file_equil.write('#!/bin/sh \n')
              
# (Don't know why I have to do this) sub.sh    
in_file = open('sub.sh',"w+")
in_file.write('#!/bin/sh \n')
              
# (Don't know why I have to do this) sub_log.sh    
in_file_log = open('sub_log.sh',"w+")
in_file_log.write('#!/bin/sh \n')
    
 