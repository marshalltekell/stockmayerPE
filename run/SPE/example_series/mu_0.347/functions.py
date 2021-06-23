#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 20:38:57 2020

@author: marshalltekell
"""

# Import packages
import numpy as np #python array manipulation package
from numpy import linalg as LA
import pandas as pd
import sys

#%% Generate random walk using Monte Carlo alogrithm

def gen_chain(l_FENE, bfl, N, L, o, e, T):  

    # Initialize int to hold the number of beads added to the chain
    k = np.int(1);

    # Initialize energy
    en = 0;
    
    # Initialize array to hold the positions of the particles
    r = np.zeros((N,3));
    ru = np.zeros((N,3));
    
    
    # Run algorithm while the number of beads is less than N
    
    while k < N:
        
        # Create a number that indicates whether bead has been added to the chain.
        add = 0;
            
        # Add another bead to the chain
        while (add < 0.1):
        
            bead = 0;
            en = 0;
            
            # Calculate the purely repulsive LJ energy for the system as is
            for i in range(0,k-1):
                
                # Sum the forces and energy felt by the ith particle do to purely repulsive LJ interactions
                for j in range(i+1,k):
                    
                    d = r[i,:] - r[j,:];
                    
                    mag = 0;
                    
                    # Apply periodic boundary conditions 
                    for l in range(0,3):
                        
                        d[l] = d[l]-L*np.round(d[l]/L);
                        mag = mag + d[l]**2;
                        
                    mag = mag**(1/2);
                    
                    # Apply cut off (2^(1/6)o)
                    if (mag < 2**(1/6)*o):
                        
                        # Sum total energy
                        temp2 = 4*e*((o/mag)**12-(o/mag)**6) + e;
                        en = en + temp2;
                        
            # Store the starting energy and reset the energy
            en1 = en;
            en = 0;
            
            # Create a while loop on the condition that the added bead does not violate the backfolding position
            while (bead < 0.1):
                            
                move = np.zeros((1,3));
                move[0,:] = [np.random.rand() - (1/2), np.random.rand() - (1/2) , np.random.rand() - (1/2)];
                
                move_mag = np.sqrt(move[0,0]**2 + move[0,1]**2 + move[0,2]**2);
                a = l_FENE/move_mag;
                move = move*a;
                
                #print(k)
                
                r[k,:] = r[k-1,:] + move[0,:];
                ru[k,:] = ru[k-1,:] + move[0,:];
                

                # Folded coordinates
                for j in range(0,3):
                    
                    temp = np.floor(r[k,j]/L);
                    r[k,j] = r[k,j] - L*temp;
      
                # Check to see if the bead violates the backfolding condition
                if (np.int(k)==np.int(1)):
                    
                    bead = 1;
                
                else: 
                    
                    d3 = r[k,:] - r[k-2,:];   
                    
                    #Apply periodic boundary conditions
                    for j in range(0,3):
                        
                        d3[j] = d3[j]-L*np.round(d3[j]/L);
                    
                    mag3 = np.sqrt(d3[0]**2 + d3[1]**2 + d3[2]**2)
                    
                    if (mag3 > bfl):
                        bead = 1;
                        
    
            # Tentatively add the bead
            k += 1;
            
            # Calculate the new forces with the added bead
            for i in range(0,k-1):
                
                # Sum the forces and energy felt by the ith particle do to purely repulsive LJ interactions
                for j in range(i+1,k):
                    
                    d = r[i,:] - r[j,:];
                    
                    mag = 0;
                    
                    # Apply periodic boundary conditions 
                    for l in range(0,3):
                        
                        d[l] = d[l]-L*np.round(d[l]/L);
                        mag = mag + d[l]**2;
                        
                    mag = mag**(1/2);
                    
                    # Apply cut off (2^(1/6)o)
                    if (mag < 2**(1/6)*o):
                        
                        # Sum total energy
                        temp2 = 4*e*((o/mag)**12-(o/mag)**6) + e;
                        en = en + temp2;
            en2 = en;
            
            # Evaluate the Metropolis condition to accept the bead
            if (en2 - en1 < 0):
                add = 1;
                
            else:
                
                #Generate random number between 0 and 1
                temp = np.random.rand();
                temp2 = np.exp((en1-en2)/(T*e))
                
                #print(temp2)
                
                if (temp2 > temp):
                    add = 1;
                    
                else:
                    k = k-1;
            
    return r, ru

#%% Generate chains on spherical nanoparticle

def gen_chain_NP(l_FENE, bfl, N, L, o, e, T, R, sites, site_num):    

    # Initialize int to hold the number of beads added to the chain
    k = np.int(1);

    # Initialize energy
    en = 0;
    
    # Initialize array to hold the positions of the particles
    r = np.zeros((N,3));
    ru = np.zeros((N,3));
    
    # Set the first bead of the chain to be on a bonding site on the NP
    r[0,:] = sites[site_num,:];
    ru[0,:] = sites[site_num,:];
    
    # Calculate shits
    
    delta_NP = (1/2)*(2*R + 2*o) - 1;
    
    # Run algorithm while the number of beads is less than N
    
    while k < N:
        
        #print(k)
        
        # Create a number that indicates whether bead has been added to the chain.
        add = 0;
            
        # Add another bead to the chain
        while (add < 0.1):
        
            bead = 0;
            en = 0;
            
            # Calculate the purely repulsive LJ energy for the system as is
            
            # Contribution from bead-bead interactions
            for i in range(0,k-1):
                
                # Sum the forces and energy felt by the ith particle do to purely repulsive LJ interactions
                for j in range(i+1,k):
                    
                    d = r[i,:] - r[j,:];
                    
                    mag = 0;
                    
                    # Apply periodic boundary conditions 
                    for l in range(0,3):
                        
                        d[l] = d[l]-L*np.round(d[l]/L);
                        mag = mag + d[l]**2;
                        
                    mag = mag**(1/2);
                    
                    # Apply cut off (2^(1/6)o)
                    if (mag < 2**(1/6)*o):
                        
                        # Sum total energy
                        temp2 = 4*e*((o/mag)**12-(o/mag)**6) + e;
                        en = en + temp2;
                        
            # Contribution from bead-NP interactions
            
            # Only if not the first bead
            
            if (k!=0):
                
                for i in range(0,k):
                    
                    d = r[i,:];
                    
                    mag = 0;
                    
                    # Apply periodic boundary conditions 
                    for l in range(0,3):
                        
                        d[l] = d[l]-L*np.round(d[l]/L);
                        mag = mag + d[l]**2;
                        
                    mag = mag**(1/2);
                    
                    # Make make sure bead is outside hard-sphere radius
                    
                    if (mag > R+o):
                    
                        # Apply cut off (2^(1/6)o)
                        if (mag < 2**(1/6)*o + delta_NP):
                            
                            # Sum total energy
                            temp3 = 4*e*((o/(mag-delta_NP))**12-(o/(mag-delta_NP))**6);
                            en = en + temp3;
                    
                    else:
                        
                        en = en + 3.53*10**7;
                        
            # Store the starting energy and reset the energy
            en1 = en;
            en = 0;
            
            # Create a while loop on the condition that the added bead does not violate the backfolding position
            while (bead < 0.1):
                            
                move = np.zeros((1,3));
                move[0,:] = [np.random.rand() - (1/2), np.random.rand() - (1/2) , np.random.rand() - (1/2)];
                
                move_mag = np.sqrt(move[0,0]**2 + move[0,1]**2 + move[0,2]**2);
                a = l_FENE/move_mag;
                move = move*a;
                
                #print(k)
                
                r[k,:] = r[k-1,:] + move[0,:];
                ru[k,:] = ru[k-1,:] + move[0,:];
                

                # Folded coordinates
                for j in range(0,3):
                    
                    temp = np.floor(r[k,j]/L);
                    r[k,j] = r[k,j] - L*temp;
      
                # Check to see if the bead violates the backfolding condition
                if (np.int(k)==np.int(1)):
                    
                    bead = 1;
                
                else: 
                    
                    d3 = r[k,:] - r[k-2,:];   
                    
                    #Apply periodic boundary conditions
                    for j in range(0,3):
                        
                        d3[j] = d3[j]-L*np.round(d3[j]/L);
                    
                    mag3 = np.sqrt(d3[0]**2 + d3[1]**2 + d3[2]**2)
                    
                    if (mag3 > bfl):
                        bead = 1;
                        
    
            # Tentatively add the bead
            k += 1;
            
            # Calculate the new forces with the added bead
            for i in range(0,k-1):
                
                # Sum the forces and energy felt by the ith particle do to purely repulsive LJ interactions
                for j in range(i+1,k):
                    
                    d = r[i,:] - r[j,:];
                    
                    mag = 0;
                    
                    # Apply periodic boundary conditions 
                    for l in range(0,3):
                        
                        d[l] = d[l]-L*np.round(d[l]/L);
                        mag = mag + d[l]**2;
                        
                    mag = mag**(1/2);
                    
                    # Apply cut off (2^(1/6)o)
                    if (mag < 2**(1/6)*o):
                        
                        # Sum total energy
                        temp2 = 4*e*((o/mag)**12-(o/mag)**6) + e;
                        en = en + temp2;
                        
            # Contribution from bead-NP interactions
                        
            for i in range(0,k):
                
                d = r[i,:];
                
                mag = 0;
                
                # Apply periodic boundary conditions 
                for l in range(0,3):
                    
                    d[l] = d[l]-L*np.round(d[l]/L);
                    mag = mag + d[l]**2;
                    
                mag = mag**(1/2);
                
                # Make make sure bead is outside hard-sphere radius
                    
                if (mag > R+o):
                
                    # Apply cut off (2^(1/6)o)
                    if (mag < 2**(1/6)*o + delta_NP):
                        
                        # Sum total energy
                        temp3 = 4*e*((o/(mag-delta_NP))**12-(o/(mag-delta_NP))**6);
                        en = en + temp3;
                        
                else:
                        en = en + 3.53*10**7;
            
            
            # Assign new energy to variable for comparison
            en2 = en;
            
            # Evaluate the Metropolis condition to accept the bead
            if (en2 - en1 < 0):
                add = 1;
                
            else:
                
                #Generate random number between 0 and 1
                temp = np.random.rand();
                temp2 = np.exp((en1-en2)/(T*e))
                
                #print(temp2)
                
                if (temp2 > temp):
                    add = 1;
                    
                else:
                    k = k-1;
            
    return r, ru

#%% Generate approximately Gaussian chains using MC algorithm, place them randomly in box

def gen_chainbox(M, r, ru_init, l_FENE, bfl, N, L, o, e, T):

    for i in range(0,M):
        
        # Use gen chain function to create M random walks
        r[i*N:(i+1)*N,:], ru_init[i*N:(i+1)*N,:] = gen_chain(l_FENE, bfl, N, L, o, e, T);
        
        # Generate a random vector
        rand3 = [np.random.rand()*L, np.random.rand()*L, np.random.rand()*L];
        
        # Add the random vector to all of the points in the chain
        for j in range(0,N):
            r[i*N+j,:] = r[i*N+j,:] + rand3;
            
        # Print the current status
        update_progress(i/M)
        
    # Restrict the chains to the simulation box
    for i in range(0,M*N):
        for j in range(0,3):
            
            temp = np.floor(r[i,j]/L);
            r[i,j] = r[i,j] - L*temp;
            
    return r;

#%% Generate approximately Gaussian chains using MC algorithm, place on lattice sites in box

def gen_chainbox_lattice(M, r, ru_init, l_FENE, bfl, N, L, o, e, T):

    # Call lattice function to generate lattice sites for each of the chains
    l = lattice_positions_sites(M, L);
    
    # Create a vector of lattice sites chosen at random without replacement
    temp = np.ceil(M**(1/3));
    rand_site = np.random.choice(np.int(temp**3), M, replace=False, p=None);
    
    
    for i in range(0,M):
        
        # Use gen chain function to create M random walks
        r[i*N:(i+1)*N,:], ru_init[i*N:(i+1)*N,:] = gen_chain(l_FENE, bfl, N, L, o, e, T);
        
        # Generate a random vector
        rand3 = l[rand_site[i],:];
        
        # Add the random vector to all of the points in the chain
        for j in range(0,N):
            r[i*N+j,:] = r[i*N+j,:] + rand3;
            
        # Print the current status
        update_progress(i/M)
        
    # Restrict the chains to the simulation box
    for i in range(0,M*N):
        for j in range(0,3):
            
            temp = np.floor(r[i,j]/L);
            r[i,j] = r[i,j] - L*temp;
            
    return r;

#%% Create randomly distributed dipoles for each of the polymer beads

def gen_mu(N, M, mu_mag_LJ):
    
    mu = np.zeros((N*M,3));
    
    for i in range(0,M*N):
        
        mu_vec = np.zeros(3);
        
        for j in range(0,3):
            
            mu_vec[j] = np.random.uniform() - 0.5;
        
        z = mu_mag_LJ/LA.norm(mu_vec);
        mu[i,:] = z*mu_vec;
        
    return mu;

#%% This function creates evenly spaced cubic lattice sites in cartesian coordinates.

def lattice_positions(N, L):
    
    temp = np.ceil(N**(1/3));
    temp2 = L/(temp);
    
    # Initialize an array that holds all of the lattice positions

    l = np.zeros((np.int(temp**3), 3));
    
    
    # Generate the x, y, z coordinates of the lattice sites

    for j in range(0,np.int(temp**3)):
        l[j,0] = temp2*(np.remainder(j,temp)+1/2);
        l[j,1] = temp2*(np.floor(j/temp)- np.floor(j/temp**2)*temp + 1/2);
        l[j,2] = temp2*((np.floor(j/temp**2)*np.floor(j/temp**2))**(1/2) + 1/2);

    # Assign N particles randomly to temp**3 lattice sites
    temp3 = np.random.choice(np.int(temp**3), N,replace=False,p=None);
    
    r = np.zeros((N, 3));
    
    for j in range(0,N):
        r[j,:] = l[temp3[j],:];
        
    return r

#%% This function creates evenly spaced cubic lattice sites in cartesian coordinates, just the sites

def lattice_positions_sites(N, L):
    
    temp = np.ceil(N**(1/3));
    temp2 = L/(temp);
    
    # Initialize an array that holds all of the lattice positions
    l = np.zeros((np.int(temp**3), 3));
    
    # Generate the x, y, z coordinates of the lattice sites
    for j in range(0,np.int(temp**3)):
        l[j,0] = temp2*(np.remainder(j,temp)+1/2);
        l[j,1] = temp2*(np.floor(j/temp)- np.floor(j/temp**2)*temp + 1/2);
        l[j,2] = temp2*((np.floor(j/temp**2)*np.floor(j/temp**2))**(1/2) + 1/2);
        
    return l

#%% For one set of 500 rows
def parse_LAMMPs_out(pathname, N, h, names_in):

    # Import data from LAMMPS output and assign to numpy array for analysis
    df = pd.read_csv(pathname,sep = '\s+', nrows = 500, skiprows = h, skip_blank_lines = True, names = names_in)
    df = df.to_numpy();
    
    return df

#%% For per atom values dumped at each time step

def parse_LAMMPs_out_int(pathname, N, h, names_in, i):

    # Import data from LAMMPS output and assign to numpy array for analysis
    df = pd.read_csv(pathname,sep = '\s+', nrows = N, skiprows = h + i*(N+h), skip_blank_lines = True, names = names_in)
    df = df.to_numpy();
    df = df[:,1:4]
    
    return df

#%% For global values fixed at each time step
    
def parse_LAMMPs_out_fix(pathname, n_t, names_in):

    # Import data from LAMMPS output and assign to numpy array for analysis
    df = pd.read_csv(pathname,sep = '\s+', nrows = n_t, skiprows = 2, skip_blank_lines = True, names = names_in)
    df = df.to_numpy();    
    return df

#%% Updates progress of a for loop

# Source: https://stackoverflow.com/questions/3160699/python-progress-bar

# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

#%% Write to read_data file for Kremer-Grest melt with one atom type and 1 bond type



def write2lammpsin_KG(N,M,L,name,r):

    in_file = open(name,"w+")
    
    # Debug
    
    '''
    N = 50;
    M = 10;
    L = 8.378836055370968;
    data = np.load('melt.npz'); 
    
    r_temp = data['a'];
    
    temp = np.zeros((M*N,3));
    
    for i in range(0,M*N):
        for j in range(0,3):
            temp[i,j] = r_temp[i,j];
            
    r = temp;
    '''
     
    # Write file
    
    in_file.write(r"LAMMPS Description")
    in_file.write("\n")
    in_file.write("\n")
    in_file.write('{:d} atoms\n'.format(N*M))
    in_file.write('{:d} bonds\n'.format(N*M-M))
    in_file.write("\n")
    in_file.write('1 atom types\n')
    in_file.write('1 bond types\n')
    in_file.write("\n")
    in_file.write('0.0 {:f} xlo xhi\n'.format(L))
    in_file.write('0.0 {:f} ylo yhi\n'.format(L))
    in_file.write('0.0 {:f} zlo zhi\n'.format(L))
    in_file.write("\n")
    in_file.write('Masses\n')
    in_file.write("\n")
    in_file.write('1 1\n')
    in_file.write("\n")
    in_file.write("Atoms\n")
    in_file.write("\n")
    
    for i in range(1, M*N+1):
        in_file.write('{:d} {:d} {:d} {:f} {:f} {:f}\n'.format(i,np.int(np.floor((i-1)/N)+1), 1, r[i-1,0], r[i-1,1], r[i-1,2]))
    
    in_file.write("\n")  
    in_file.write("Bonds\n")
    in_file.write("\n")
    
    count = 0;
    
    for i in range(1, M*N+1):
        
        if ((i)/N != np.round(i/N)):
            
            in_file.write('{:d} {:d} {:d} {:d}\n'.format(np.int(i-count), 1, i, np.int(i+1)))
        else:
            
            count +=1


#%% Write to read_data file for Kremer-Grest melt with freely rotating dipoles and salt
            
def write2lammpsin_elec(N, M, I, L, name, masses, r, rc, ra, mu, q_c, q_a):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write(r"LAMMPS Description")
    in_file.write("\n")
    in_file.write("\n")
    in_file.write('{:d} atoms\n'.format(N*M+2*I))
    in_file.write('{:d} bonds\n'.format(N*M-M))
    in_file.write("\n")
    in_file.write('3 atom types\n')
    in_file.write('1 bond types\n')
    in_file.write("\n")
    in_file.write('0.0 {:f} xlo xhi\n'.format(L))
    in_file.write('0.0 {:f} ylo yhi\n'.format(L))
    in_file.write('0.0 {:f} zlo zhi\n'.format(L))
    in_file.write("\n")
    in_file.write('Masses\n')
    in_file.write("\n")
    in_file.write('1 {:f}\n'.format(masses[0]))
    in_file.write('2 {:f}\n'.format(masses[1]))
    in_file.write('3 {:f}\n'.format(masses[2]))
    in_file.write("\n")
    in_file.write("Atoms\n")
    in_file.write("\n")
    
    # The atom attributes have the following format: atom_ID atom_type x y z diameter density q mux muy muz molecule_ID
    
    # Polymer beads
    for i in range(1, M*N+1):
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:f} {:d} {:.30f} {:.30f} {:.30f} {:d}\n'.format(i, 1, r[i-1,0], r[i-1,1], r[i-1,2], 0, 1, 0, mu[i-1,0], mu[i-1,1], mu[i-1,2], np.int(np.floor((i-1)/N)+1)))
    
    # Cations
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:f} {:.19f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 2, rc[i-1,0], rc[i-1,1], rc[i-1,2], 0, 1, q_c, 0, 0, 0, molecule_ID))
    
    # Anions
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + I + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + I + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:f} {:.19f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 3, ra[i-1,0], ra[i-1,1], ra[i-1,2], 0, 1, q_a, 0, 0, 0, molecule_ID))
    
    # Bonds
    in_file.write("\n")  
    in_file.write("Bonds\n")
    in_file.write("\n")
    
    count = 0;
    
    for i in range(1, M*N+1):
        
        if ((i)/N != np.round(i/N)):
            
            in_file.write('{:d} {:d} {:d} {:d}\n'.format(np.int(i-count), 1, i, np.int(i+1)))
        
        else:
            
            count +=1

#%%
def write2lammpsin_elec2(N, M, I, boxlo, boxhi, name, masses, r, rc, ra, mu, q_c, q_a):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write(r"LAMMPS Description")
    in_file.write("\n")
    in_file.write("\n")
    in_file.write('{:d} atoms\n'.format(N*M+2*I))
    in_file.write('{:d} bonds\n'.format(N*M-M))
    in_file.write("\n")
    in_file.write('3 atom types\n')
    in_file.write('1 bond types\n')
    in_file.write("\n")
    in_file.write('{:f} {:f} xlo xhi\n'.format(boxlo[0],boxhi[0]))
    in_file.write('{:f} {:f} ylo yhi\n'.format(boxlo[1],boxhi[1]))
    in_file.write('{:f} {:f} zlo zhi\n'.format(boxlo[2],boxhi[2]))
    in_file.write("\n")
    in_file.write('Masses\n')
    in_file.write("\n")
    in_file.write('1 {:f}\n'.format(masses[0]))
    in_file.write('2 {:f}\n'.format(masses[1]))
    in_file.write('3 {:f}\n'.format(masses[2]))
    in_file.write("\n")
    in_file.write("Atoms\n")
    in_file.write("\n")
    
    # The atom attributes have the following format: atom_ID atom_type x y z diameter density q mux muy muz molecule_ID
    
    # Polymer beads
    for i in range(1, M*N+1):
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:d} {:d} {:d} {:.3f} {:.3f} {:.3f} {:d}\n'.format(i, 1, r[i-1,0], r[i-1,1], r[i-1,2], 0, 1, 0, mu[i-1,0], mu[i-1,1], mu[i-1,2], np.int(np.floor((i-1)/N)+1)))
    
    # Cations
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:d} {:d} {:.3f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 2, rc[i-1,0], rc[i-1,1], rc[i-1,2], 0, 1, q_c, 0, 0, 0, molecule_ID))
    
    # Anions
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + I + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + I + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:d} {:d} {:.3f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 3, ra[i-1,0], ra[i-1,1], ra[i-1,2], 0, 1, q_a, 0, 0, 0, molecule_ID))
    
    # Bonds
    in_file.write("\n")  
    in_file.write("Bonds\n")
    in_file.write("\n")
    
    count = 0;
    
    for i in range(1, M*N+1):
        
        if ((i)/N != np.round(i/N)):
            
            in_file.write('{:d} {:d} {:d} {:d}\n'.format(np.int(i-count), 1, i, np.int(i+1)))
        
        else:
            
            count +=1
            
#%% Write to read_data file for Kremer-Grest melt with freely rotating dipoles and salt
            
def write2lammpsin_elec_angle(N, M, I, L, name, masses, r, rc, ra, mu, q_c, q_a):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write(r"LAMMPS Description")
    in_file.write("\n")
    in_file.write("\n")
    in_file.write('{:d} atoms\n'.format(N*M+2*I))
    in_file.write('{:d} bonds\n'.format(N*M-M))
    in_file.write('{:d} angles\n'.format(M*(N-2)))
    in_file.write("\n")
    in_file.write('3 atom types\n')
    in_file.write('1 bond types\n')
    in_file.write('1 angle types\n')
    in_file.write("\n")
    in_file.write('0.0 {:f} xlo xhi\n'.format(L))
    in_file.write('0.0 {:f} ylo yhi\n'.format(L))
    in_file.write('0.0 {:f} zlo zhi\n'.format(L))
    in_file.write("\n")
    in_file.write('Masses\n')
    in_file.write("\n")
    in_file.write('1 {:f}\n'.format(masses[0]))
    in_file.write('2 {:f}\n'.format(masses[1]))
    in_file.write('3 {:f}\n'.format(masses[2]))
    in_file.write("\n")
    in_file.write("Atoms\n")
    in_file.write("\n")
    
    # The atom attributes have the following format: atom_ID atom_type x y z diameter density q mux muy muz molecule_ID
    
    # Polymer beads
    for i in range(1, M*N+1):
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:d} {:d} {:.3f} {:.3f} {:.3f} {:d}\n'.format(i, 1, r[i-1,0], r[i-1,1], r[i-1,2], 0, 1, 0, mu[i-1,0], mu[i-1,1], mu[i-1,2], np.int(np.floor((i-1)/N)+1)))
    
    # Cations
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:d} {:.1f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 2, rc[i-1,0], rc[i-1,1], rc[i-1,2], 0, 1, q_c, 0, 0, 0, molecule_ID))
    
    # Anions
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + I + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + I + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:d} {:.1f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 3, ra[i-1,0], ra[i-1,1], ra[i-1,2], 0, 1, q_a, 0, 0, 0, molecule_ID))
    
    # Bonds
    in_file.write("\n")  
    in_file.write("Bonds\n")
    in_file.write("\n")
    
    count = 0;
    
    for i in range(1, M*N+1):
        
        if ((i)/N != np.round(i/N)):
            
            in_file.write('{:d} {:d} {:d} {:d}\n'.format(np.int(i-count), 1, i, np.int(i+1)))
        
        else:
            
            count +=1
            
    # Angles
    in_file.write("\n")  
    in_file.write("Angles\n")
    in_file.write("\n")
    
    count = 0;
    
    for i in range(1, M*N+1):
        
        if ((i)/N != np.round(i/N) and (i+1)/N != np.round((i+1)/N)):
            
            in_file.write('{:d} {:d} {:d} {:d} {:d}\n'.format(np.int(i-count), 1, i, i+1, i+2))
        
        else:

            count +=1

#%% Write to read_data file for Kremer-Grest melt with freely rotating dipoles, salt, angles and spheres
            
def write2lammpsin_elec_angle_sphere(N, M, I, L, name, masses, dens, diam, r, rc, ra, mu, q_c, q_a):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write(r"LAMMPS Description")
    in_file.write("\n")
    in_file.write("\n")
    in_file.write('{:d} atoms\n'.format(N*M+2*I))
    in_file.write('{:d} bonds\n'.format(N*M-M))
    in_file.write('{:d} angles\n'.format(M*(N-2)))
    in_file.write("\n")
    in_file.write('3 atom types\n')
    in_file.write('1 bond types\n')
    in_file.write('1 angle types\n')
    in_file.write("\n")
    in_file.write('0.0 {:f} xlo xhi\n'.format(L))
    in_file.write('0.0 {:f} ylo yhi\n'.format(L))
    in_file.write('0.0 {:f} zlo zhi\n'.format(L))
    in_file.write("\n")
    in_file.write('Masses\n')
    in_file.write("\n")
    in_file.write('1 {:f}\n'.format(masses[0]))
    in_file.write('2 {:f}\n'.format(masses[1]))
    in_file.write('3 {:f}\n'.format(masses[2]))
    in_file.write("\n")
    in_file.write("Atoms\n")
    in_file.write("\n")
    
    # The atom attributes have the following format: atom_ID atom_type x y z diameter density q mux muy muz molecule_ID
    
    # Polymer beads
    for i in range(1, M*N+1):
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:f} {:d} {:.3f} {:.3f} {:.3f} {:d}\n'.format(i, 1, r[i-1,0], r[i-1,1], r[i-1,2], diam[0], dens[0], 0, mu[i-1,0], mu[i-1,1], mu[i-1,2], np.int(np.floor((i-1)/N)+1)))
    
    # Cations
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:f} {:.1f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 2, rc[i-1,0], rc[i-1,1], rc[i-1,2], diam[1], dens[1], q_c, 0, 0, 0, molecule_ID))
    
    # Anions
    for i in range(1, I+1):
        
        # Find atom_ID
        temp = i + I + M*N;
        atom_ID = np.int(temp);
        
        # Find molecule_ID
        temp = i + I + M;
        molecule_ID = np.int(temp);
        
        in_file.write('{:d} {:d} {:f} {:f} {:f} {:f} {:f} {:.1f} {:d} {:d} {:d} {:d}\n'.format(atom_ID, 3, ra[i-1,0], ra[i-1,1], ra[i-1,2], diam[2], dens[2], q_a, 0, 0, 0, molecule_ID))
    
    # Bonds
    in_file.write("\n")  
    in_file.write("Bonds\n")
    in_file.write("\n")
    
    count = 0;
    
    for i in range(1, M*N+1):
        
        if ((i)/N != np.round(i/N)):
            
            in_file.write('{:d} {:d} {:d} {:d}\n'.format(np.int(i-count), 1, i, np.int(i+1)))
        
        else:
            
            count +=1
            
    # Angles
    in_file.write("\n")  
    in_file.write("Angles\n")
    in_file.write("\n")
    
    count = 0;
    
    for i in range(1, M*N+1):
        
        if ((i)/N != np.round(i/N) and (i+1)/N != np.round((i+1)/N)):
            
            in_file.write('{:d} {:d} {:d} {:d} {:d}\n'.format(np.int(i-count), 1, i, i+1, i+2))
        
        else:

            count +=1
            