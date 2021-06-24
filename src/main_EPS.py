#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 15:06:48 2020

@author: marshalltekell

NVE MD with only Lennard-Jones interactions
"""
# Equation of state of the Lennard-Jones fluid

# Import packages
import matplotlib.pyplot as plt #python plotting package
import numpy as np #python array manipulation package
from numpy import linalg as LA
#from scipy import stats
#from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 

# Import parameters
from parameters import kb, T, V, n_t, M, N, e_c, time;

#%% Additional parameters

Vm = V*1e-30; # Ang3 to m3
e_0 = 8.854*1e-12; # vacuum permittivity

# Additional extension
#n_t_ext = np.int(50000000/10000);
#time = np.arange(0,np.int(2*50e+6 + 15e+6),np.int(1e+4));

#%% Load in data from Terremoto zip

# Original data
data_mu = np.load('run_data_mu_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int

# Extension data
#data_mu_ext = np.load('run_data_mu_ext_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int

# Extension data 2
#data_mu_ext2 = np.load('run_data_mu_ext2_lammps_ELEC.npz');  # a=r_int, b=rc_int, c=ra_int

M_int = np.zeros((n_t,3));
M1 = np.zeros((n_t,3));
M2 = np.zeros((n_t));

# Assign position data
mu_int = data_mu['a'];
#mu_int_ext = data_mu_ext['a'];
#mu_int_ext2 = data_mu_ext2['a'];

# Stitch together data sets
#mu_int = np.concatenate((mu_int, mu_int_ext), axis=1);
#mu_int = np.concatenate((mu_int, mu_int_ext2), axis=1);

#%%

# Get M
for i in range(0,n_t):
        
    mt = np.zeros((3));
    
    for j in range(0,N*M):
        
        mt += mu_int[j,3*i:3*(i+1)];
    
    M1[i] = mt;

# Get M**2
for i in range(0,n_t):
    
    t1 = 0.0;
    
    for j in range(0,N*M):
        
        t2 = 0.0;
        
        for k in range(0,3):
        
            t2 += mu_int[j,i*3+k]*mu_int[j,i*3+k];
        
        t1 += t2;
    
    M2[i] = t1;

# Get averages
M1_av = np.zeros((3));

for i in range(0,3):
    
    M1_av[i] = np.mean(M1[:,i]);

M1_av = np.dot(M1_av, M1_av);
M2_av = np.mean(M2);
Msc = (M2_av-M1_av**2)*(e_c**2)*(1e-20);

ep = 1 + ((4*np.pi)/(3*e_0*kb*T*Vm))*(Msc)
    
#%% Some prep for the upcoming figures

# Independent variables for plots
index = np.arange(0,N);

# Adjust the font sizes globally for all of the legends, titles, and labels

SMALL_SIZE = 8;
MEDIUM_SIZE = 10;
BIGGER_SIZE = 12;


plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


plt.rc('font',family='Arial')

#%% Figure 1

# Set figure height and width

width = 5.6;
height = 4.2;

fig1 = plt.figure(1,(width,height))

# Put inward-facing major and minor tickmarks with defined width
ax = fig1.add_subplot(111)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
  
ax.minorticks_on()
ax.tick_params(axis='x', which='minor', direction='in',width=1.2)
ax.tick_params(axis="x", which='major',direction="in",width=1.2)
ax.tick_params(axis='y', which='minor', direction='in',width=1.2)
ax.tick_params(axis="y", which='major',direction="in",width=1.2)

# Set the font name for axis tick labels to be Arial
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Write the labels for the x and y axes
plt.xlabel("\u03BC$_{z}$", fontname="Arial", fontsize=12);
plt.ylabel("#", fontname="Arial", fontsize=12);

# Set the limits for both axes
#plt.ylim(-0.50,1.50);
#plt.xlim(0, 5);

plt.hist(mu_int[:,-1],bins=100)

plt.show()


# Save the figure to a .pdf file
fig1.savefig('muz_hist.pdf')

plt.show()

#%% Figure 2

# Set figure height and width

width = 5.6;
height = 4.2;

fig2 = plt.figure(2,(width,height))

# Put inward-facing major and minor tickmarks with defined width
ax = fig2.add_subplot(111)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)
  
ax.minorticks_on()
ax.tick_params(axis='x', which='minor', direction='in',width=1.2)
ax.tick_params(axis="x", which='major',direction="in",width=1.2)
ax.tick_params(axis='y', which='minor', direction='in',width=1.2)
ax.tick_params(axis="y", which='major',direction="in",width=1.2)

# Set the font name for axis tick labels to be Arial
for tick in ax.get_xticklabels():
    tick.set_fontname("Arial")
for tick in ax.get_yticklabels():
    tick.set_fontname("Arial")

# Write the labels for the x and y axes
plt.xlabel("t", fontname="Arial", fontsize=12);
plt.ylabel("M$_z$", fontname="Arial", fontsize=12);

# Set the limits for both axes
#plt.ylim(-0.50,1.50);
#plt.xlim(0, 5);

# Plot data
plt.plot(time,M1[:,2])

plt.show()


# Save the figure to a .pdf file
fig2.savefig('Mz_int.pdf')

plt.show()
