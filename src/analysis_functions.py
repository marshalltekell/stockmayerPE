#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 12:27:33 2020

@author: marshalltekell
"""
import numpy as np #python array manipulation package
from numpy import linalg as LA
import sys
import pandas as pd

#%% General RDF for LJ fluid or KG fluid in LJ units

def RDF(N, L, nbins, r, rho,o):
    
    # Intialize RDF
    g = np.zeros((nbins,1))
    
    # Calculate the bin width over the range 0 to L/2
    dg = L/(2*(nbins-1));
    
    for i in range(0,N-1):
        for j in range(i+1,N):
            
            d = r[i,:] - r[j,:];
            
            
            # Apply periodic boundary conditions
            
            mag = 0;
            
            
            for k in range(0,3):
                
                d[k] = d[k]-L*np.round(d[k]/L)
    
                mag = mag + d[k]**2;
                
            mag = mag**(1/2);
            
            # Apply cut off (half the box size) and add tally to bin if magnitude of interparticle distance falls within bin
            
            if (mag < L/2):
                
                k = np.int(np.round(mag/dg));
                g[k] = g[k] + 2;
    
    # Normalize RDF to ideal gas
    
    for i in range(0,nbins):
        
        dV = ((i+1)**3-i**3)*dg**3;
        rho_id = ((4/3)*3.1495*dV*rho)/(o**(3));
        g[i] = g[i]/(N*rho_id);

    # Return the RDF
    
    return g

#%% RDF taking into account real units for polymer portion of SPE

def RDF_p(N, L, nbins, r, rho, o, m):
    
    # Intialize RDF
    g = np.zeros((nbins,1))
    
    # Calculate the bin width over the range 0 to L/2
    dg = L/(2*(nbins-1));
    
    for i in range(0,N-1):
        for j in range(i+1,N):
            
            d = r[i,:] - r[j,:];
            
            
            # Apply periodic boundary conditions
            
            mag = 0;
            
            
            for k in range(0,3):
                
                d[k] = d[k]-L*np.round(d[k]/L)
    
                mag = mag + d[k]**2;
                
            mag = mag**(1/2);
            
            # Apply cut off (half the box size) and add tally to bin if magnitude of interparticle distance falls within bin
            
            if (mag < L/2):
                
                k = np.int(np.round(mag/dg));
                g[k] = g[k] + 2;
    
    # Normalize RDF to ideal gas
    
    for i in range(0,nbins):
        
        dV = (4/3)*np.pi*((i+1)**3-i**3)*dg**3;
        g[i] = g[i]*m*((10**24)/(6.022*10**(23)))*(1/dV)*(1/(N*rho))

    # Return the RDF
    
    return g

#%% RDF taking into account real units for cations of SPE

def RDF_c(N, L, nbins, r, rho, o, m):
    
    # Intialize RDF
    g = np.zeros((nbins,1))
    
    # Calculate the bin width over the range 0 to L/2
    dg = L/(2*(nbins-1));
    
    for i in range(0,N-1):
        for j in range(i+1,N):
            
            d = r[i,:] - r[j,:];
            
            
            # Apply periodic boundary conditions
            
            mag = 0;
            
            
            for k in range(0,3):
                
                d[k] = d[k]-L*np.round(d[k]/L)
    
                mag = mag + d[k]**2;
                
            mag = mag**(1/2);
            
            # Apply cut off (half the box size) and add tally to bin if magnitude of interparticle distance falls within bin
            
            if (mag < L/2):
                
                k = np.int(np.round(mag/dg));
                g[k] = g[k] + 2;
    
    # Normalize RDF to ideal gas
    
    for i in range(0,nbins):
        
        dV = (4/3)*np.pi*((i+1)**3-i**3)*dg**3;
        g[i] = g[i]*m*((10**24)/(6.022*10**(23)))*(1/dV)*(1/(N*rho))

    # Return the RDF
    
    return g

#%% RDF taking into account real units for cations of SPE

def RDF_ca(I, L, nbins, r, rho, m):
    
    # Intialize RDF
    g = np.zeros((nbins,1))
    
    # Calculate the bin width over the range 0 to L/2
    dg = L/(2*(nbins-1));
    
    for i in range(0,I):
        for j in range(I,2*I):
            
            d = r[i,:] - r[j,:];
            
            # Apply periodic boundary conditions
            
            mag = 0;
            
            for k in range(0,3):
                
                d[k] = d[k]-L*np.round(d[k]/L)
    
                mag = mag + d[k]**2;
                
            mag = mag**(1/2);
            
            # Apply cut off (half the box size) and add tally to bin if magnitude of interparticle distance falls within bin
            
            if (mag < L/2):
                
                k = np.int(np.round(mag/dg));
                g[k] = g[k] + 1;
    
    # Normalize RDF to ideal gas
    
    for i in range(0,nbins):
        
        dV = (4/3)*np.pi*((i+1)**3-i**3)*dg**3;
        g[i] = (g[i]*m*10**24)/(I*rho*dV*6.022*10**(23));

    # Return the RDF
    
    return g


#%% RDF for intra chain

def RDF_intra(N, L, nbins, r, rho, o, m, Nc):
    
    # Intialize RDF
    g = np.zeros((nbins,1))
    
    # Calculate the bin width over the range 0 to L/2
    dg = L/(2*(nbins-1));
    
    for i in range(0,N-1):
        for j in range(i+1,N):
            
            #check to see if on the same chain
            mol_i = np.floor(i/Nc);
            mol_j = np.floor(j/Nc);
            
            if (mol_i == mol_j):
                
            
                d = r[i,:] - r[j,:];
                
                
                # Apply periodic boundary conditions
                
                mag = 0;
                
                
                for k in range(0,3):
                    
                    d[k] = d[k]-L*np.round(d[k]/L)
        
                    mag = mag + d[k]**2;
                    
                mag = mag**(1/2);
                
                # Apply cut off (half the box size) and add tally to bin if magnitude of interparticle distance falls within bin
                
                if (mag < L/2):
                    
                    k = np.int(np.round(mag/dg));
                    g[k] = g[k] + 2;
    
    # Normalize RDF to ideal gas
    
    for i in range(0,nbins):
        
        dV = (4/3)*np.pi*((i+1)**3-i**3)*dg**3;
        g[i] = g[i]*m*((10**24)/(6.022*10**(23)))*(1/dV)*(1/(N*rho))

    # Return the RDF
    
    return g

#%% RDF for inter chain

def RDF_inter(N, L, nbins, r, rho, o, m, Nc):
    
    # Intialize RDF
    g = np.zeros((nbins,1))
    
    # Calculate the bin width over the range 0 to L/2
    dg = L/(2*(nbins-1));
    
    for i in range(0,N-1):
        for j in range(i+1,N):
            
            #check to see if on the same chain
            mol_i = np.floor(i/Nc);
            mol_j = np.floor(j/Nc);
            
            if (mol_i != mol_j):
            
                d = r[i,:] - r[j,:];
                
                
                # Apply periodic boundary conditions
                
                mag = 0;
                
                
                for k in range(0,3):
                    
                    d[k] = d[k]-L*np.round(d[k]/L)
        
                    mag = mag + d[k]**2;
                    
                mag = mag**(1/2);
                
                # Apply cut off (half the box size) and add tally to bin if magnitude of interparticle distance falls within bin
                
                if (mag < L/2):
                    
                    k = np.int(np.round(mag/dg));
                    g[k] = g[k] + 2;
    
    # Normalize RDF to ideal gas
    
    for i in range(0,nbins):
        
        dV = (4/3)*np.pi*((i+1)**3-i**3)*dg**3;
        g[i] = g[i]*m*((10**24)/(6.022*10**(23)))*(1/dV)*(1/(N*rho))

    # Return the RDF
    
    return g

#%% Given the unfolded polymer chains, calculate the gyration tensor and the end to end distance

def calc_RG(N, M, L, r):
    
    # Initialize arrays
    rcm = np.zeros((M,3));
    rg = np.zeros((M,7));
    d = np.zeros((3));
    
    # Find the center of mass
    for i in range(0,M):
        
        sum1 = [0, 0, 0];
        
        for j in range(0,N):
            
            sum1 += r[i*N + j,:];
            
        rcm[i,:] = sum1/N;
    
    # For each of the chains, calculate the gyration tensor
    for i in range(0,M):
        
        # Initialize arrays and sums
        s1 = 0;
        s2 = 0;
        s3 = np.zeros((3,3));
        
        
        for j in range(0,N):
            
            mag = 0;
            
            for k in range(0,3):

                d[k] = r[i*N + j,k] - rcm[i,k];
                mag += d[k]**2;
            
            s1 += mag;
            s2 += mag**2;
            
            # Find gyration tensor
            for k in range(0,3):
            
                for l in range(0,3):
                    
                    s3[k,l] += d[k]*d[l];
        
        # Find eigenvalues of gyration tensor matrix
        w, _ = LA.eig(s3);
        w = np.sort(w);
            
        # Store values    
        rg[i,0] = s1/N;
        rg[i,1] = s2/N;
        rg[i,4:7] = w;
                
    # Find the average squared and quartic end-to-end distance of the chain    
    for i in range(0,M):
        
        mag = 0;
          
        for j in range(0, 3):
            
            d[j] = r[i*N,j] - r[i*N + N-1,j];
            mag += d[j]**2;
            
        rg[i,2] = mag;
        rg[i,3] = mag**2;

    return rg, w, s3;

#%% Find where the polymer chain has been folded into the box and record image flags in count

def calc_RG2(N, temp, L, M):
        
    count = np.zeros((N*M,3));
    
    for i in range(0,M):
        
        for k in range(0,N-1):
            
            for j in range(0,3):
                
                d = temp[i*N + k + 1, j] - temp[i*N + k, j];
                
                if (d < -(1/2)*L):
                    
                    count[i*N + k +1:(i+1)*N,j] += 1.0;
                    
                if (d > (1/2)*L):
                    
                    count[i*N + k +1:(i+1)*N,j] -= 1.0;
        

    return count

#%% Calculate the MSD 

def calc_MSD3(xu_int, yu_int, zu_int, N, n_t, L):
    
        
    # Initialize arrays
    dx = np.zeros((N,n_t));
    dy = np.zeros((N,n_t));
    dz = np.zeros((N,n_t));
    count = np.zeros((n_t,1));
    x_av = np.zeros((n_t,1));
    y_av = np.zeros((n_t,1));
    z_av = np.zeros((n_t,1));
    
    # Start the first loop over the number of particles   
    for i in range(0,N):
        
        # After choosing a particle, loop over the number of times
        for j in range(1,n_t):
            
            dx[i,j] = xu_int[i,j]-xu_int[i,j-1];
            dy[i,j] = yu_int[i,j]-yu_int[i,j-1];
            dz[i,j] = zu_int[i,j]-zu_int[i,j-1];
            
    #Loop over particles and times again        
    for i in range(0,N):

        for j in range(0,n_t):
            
            #Initialize total displacements from this specific time origin
            tsx = 0.0;
            tsy = 0.0;
            tsz = 0.0;
            
            # For particular time, loop over all successive times
            for k in range(j+1, n_t-1):
                
                # Count the number of times included from each time origin
                count[k-j,0] += 1;
                
                 # Square the cumulative displacement  
                tsx += dx[i,k];
                tsy += dy[i,k];
                tsz += dz[i,k];

                # Assign the cumulative displacement to the MSD at the difference between this time and the time origin
                x_av[k-j,0] += tsx**2;
                y_av[k-j,0] += tsy**2;
                z_av[k-j,0] += tsz**2;
       
        # Update progress bar for MSD
        update_progress(i/N)
          
    # Divide the average by the number of points and the number of origins
    
    x_av = np.divide(x_av, count, out=np.zeros_like(x_av), where=count!=0);
    y_av = np.divide(y_av, count, out=np.zeros_like(y_av), where=count!=0);
    z_av = np.divide(z_av, count, out=np.zeros_like(z_av), where=count!=0);                 
    
  
    return x_av, y_av, z_av, count

#%% Calculate the MSD (^4)

def calc_MSD4(xu_int, yu_int, zu_int, N, n_t, L):
    
        
    # Initialize arrays
    dx = np.zeros((N,n_t));
    dy = np.zeros((N,n_t));
    dz = np.zeros((N,n_t));
    count = np.zeros((n_t,1));
    x_av = np.zeros((n_t,1));
    y_av = np.zeros((n_t,1));
    z_av = np.zeros((n_t,1));
    
    # Start the first loop over the number of particles   
    for i in range(0,N):
        
        # After choosing a particle, loop over the number of times
        for j in range(1,n_t):
            
            dx[i,j] = xu_int[i,j]-xu_int[i,j-1];
            dy[i,j] = yu_int[i,j]-yu_int[i,j-1];
            dz[i,j] = zu_int[i,j]-zu_int[i,j-1];
            
    #Loop over particles and times again        
    for i in range(0,N):

        for j in range(0,n_t):
            
            #Initialize total displacements from this specific time origin
            tsx = 0.0;
            tsy = 0.0;
            tsz = 0.0;
            
            # For particular time, loop over all successive times
            for k in range(j+1, n_t-1):
                
                # Count the number of times included from each time origin
                count[k-j,0] += 1;
                
                 # Square the cumulative displacement  
                tsx += dx[i,k];
                tsy += dy[i,k];
                tsz += dz[i,k];

                # Assign the cumulative displacement to the MSD at the difference between this time and the time origin
                x_av[k-j,0] += tsx**4;
                y_av[k-j,0] += tsy**4;
                z_av[k-j,0] += tsz**4;
       
        # Update progress bar for MSD
        update_progress(i/N)
          
    # Divide the average by the number of points and the number of origins
    
    x_av = np.divide(x_av, count, out=np.zeros_like(x_av), where=count!=0);
    y_av = np.divide(y_av, count, out=np.zeros_like(y_av), where=count!=0);
    z_av = np.divide(z_av, count, out=np.zeros_like(z_av), where=count!=0);                 
    
  
    return x_av, y_av, z_av, count

#%% Updates progress of a for loop

# Source: 

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

#%% 
def write2xyz(N, M, I, name, r, rc, ra):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write('{:d}\n'.format(N+2*I))
    in_file.write("Test\n")

    
    # Polymer beads
    for i in range(1, N+1):
        
        in_file.write('C {:f} {:f} {:f}\n'.format(r[i-1,0], r[i-1,1], r[i-1,2]))
    
    # Cations
    for i in range(1, I+1):
                
        in_file.write('Li {:f} {:f} {:f}\n'.format(rc[i-1,0], rc[i-1,1], rc[i-1,2]))
    
    # Anions
    for i in range(1, I+1):
                
        in_file.write('Cl {:f} {:f} {:f}\n'.format(ra[i-1,0], ra[i-1,1], ra[i-1,2]))     

#%% Write to .xyz file (just the ions)
            
def writeion2xyz(N, M, I, name, rc, ra):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write('{:d}\n'.format(2*I))
    in_file.write("Test\n")

    '''
    # Polymer beads
    for i in range(1, 41):
        
        in_file.write('C {:f} {:f} {:f}\n'.format(r[i-1,0], r[i-1,1], r[i-1,2]))
    '''
    
    # Cations
    for i in range(1, I+1):
                
        in_file.write('Li {:f} {:f} {:f}\n'.format(rc[i-1,0], rc[i-1,1], rc[i-1,2]))
    
    # Anions
    for i in range(1, I+1):
                
        in_file.write('Cl {:f} {:f} {:f}\n'.format(ra[i-1,0], ra[i-1,1], ra[i-1,2]))
        
#%% Write to .xyz file (GNP)
            
def write2xyzGNP(N, M, name, r, L):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write('{:d}\n'.format(M*N+1))
    in_file.write("Test\n")
    
    in_file.write('O {:f} {:f} {:f}\n'.format(L/2, L/2, L/2))

    
    # Polymer beads
    for i in range(1, M*N+1):
        
        in_file.write('C {:f} {:f} {:f}\n'.format(r[i-1,0], r[i-1,1], r[i-1,2]))
        
#%% Write to .xyz (GNP with ions)
    
def write2xyzGNPion(N, M, I, name, r, L, rc, ra):

    in_file = open(name,"w+")
     
    # Write file
    
    in_file.write('{:d}\n'.format(M*N+2*I+1))
    in_file.write("Test\n")
    
    in_file.write('O {:f} {:f} {:f}\n'.format(L/2, L/2, L/2))

    
    # Polymer beads
    for i in range(1, M*N+1):
        
        in_file.write('C {:f} {:f} {:f}\n'.format(r[i-1,0], r[i-1,1], r[i-1,2]))
        
    # Cations
    for i in range(1, I+1):
                
        in_file.write('Li {:f} {:f} {:f}\n'.format(rc[i-1,0], rc[i-1,1], rc[i-1,2]))
    
    # Anions
    for i in range(1, I+1):
                
        in_file.write('Cl {:f} {:f} {:f}\n'.format(ra[i-1,0], ra[i-1,1], ra[i-1,2]))
        
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

#%% Degree of uncorrelated ion motion
    
def duim(xu_int, yu_int, zu_int, N, n_t, L):
    
    # Initialize arrays
    dx = np.zeros((N,n_t));
    dy = np.zeros((N,n_t));
    dz = np.zeros((N,n_t));
    count = np.zeros((n_t,1));
    x_av = np.zeros((n_t,1));
    y_av = np.zeros((n_t,1));
    z_av = np.zeros((n_t,1));
    
    # Start the first loop over the number of particles   
    for i in range(0,N):
        
        # After choosing a particle, loop over the number of times
        for j in range(1,n_t):
            
            dx[i,j] = xu_int[i,j]-xu_int[i,j-1];
            dy[i,j] = yu_int[i,j]-yu_int[i,j-1];
            dz[i,j] = zu_int[i,j]-zu_int[i,j-1];
            
    #Loop over particles and times again        
    for i in range(0,N):
        
        for j in range(0,N):

            for k in range(0,n_t):
                
                #Initialize total displacements from this specific time origin
                tsxi = 0.0;
                tsyi = 0.0;
                tszi = 0.0;
                tsxj = 0.0;
                tsyj = 0.0;
                tszj = 0.0;
                
                # For particular time, loop over all successive times
                for l in range(k+1, n_t-1):
                    
                    # Count the number of times included from each time origin
                    count[l-k,0] += 1;
                    
                     # Square the cumulative displacement  
                    tsxi += dx[i,l];
                    tsyi += dy[i,l];
                    tszi += dz[i,l];
                    
                    tsxj += dx[j,l];
                    tsyj += dy[j,l];
                    tszj += dz[j,l];
    
                    # Assign the cumulative displacement to the MSD at the difference between this time and the time origin
                    x_av[l-k,0] += tsxi*tsxj;
                    y_av[l-k,0] += tsyi*tsyj;
                    z_av[l-k,0] += tszi*tszj;
       
        
        # Update progress bar for MSD
        #update_progress(i/N)
          
    # Divide the average by the number of points and the number of origins
    
    x_av = np.divide(x_av, count, out=np.zeros_like(x_av), where=count!=0);
    y_av = np.divide(y_av, count, out=np.zeros_like(y_av), where=count!=0);
    z_av = np.divide(z_av, count, out=np.zeros_like(z_av), where=count!=0);                 
    
  
    return x_av, y_av, z_av, count

#%% Generate set of random q vectors for static structure factor
    
def gen_q(low, high, q_mag, num_p, num_q):
    
    
    q_run = np.zeros((num_p,3*num_q));
    q_vec = np.zeros((1,3));
    for i in range(0,num_p):
        
        for j in range(0,num_q):
            
            q_vec[0,:] = [np.random.uniform() - 0.5, np.random.uniform() - 0.5, np.random.uniform() - 0.5];
            z = q_mag[i]/LA.norm(q_vec[0,:]);
            q_vec[0,:] = z*q_vec[0,:];
            
            q_run[i,3*j:3*j+3] = q_vec[0,:]
            
    return q_run

#%% Calculates the auto-correlation function for the dipole moment in the system
    
def musum(M_int, N, n_t):
    
    M_acf = np.zeros((n_t,1));
    M_t = np.zeros((3,1));
    M_0 = np.zeros((3,1));
    count = np.zeros((n_t,1));
    
    # Loop over the number of particles, the number of time origins, and the successive times
    
    for i in range(0,n_t):
        
        for j in range(i, n_t-1):
            
            count[j-i,0] += 1;
            
            M_0 = M_int[i,:];
            M_t = M_int[j,:];            
            M_acf[j-i,0] += np.dot(M_0,M_t);
                
    M_acf = np.divide(M_acf, count, out=np.zeros_like(M_acf), where=count!=0);
    
    return M_acf;

#%% Calculates the auto-correlation function for the dipole moment in the system--log time sampling (no time origins)
    
def musum_log(M_int, N, n_t):
    
    # Initialize arrays
    M_acf = np.zeros((n_t));
    M_t = np.zeros((3,1));
    M_0 = np.zeros((3,1));
    
    # Get the initial overall dipole
    M_0 = M_int[0,:];
    
    # Loop over the number of particles, the number of time origins, and the successive times
    
    
    for i in range(0,n_t):
        
        M_t = M_int[i,:];
        M_acf[i] = np.dot(M_0,M_t);
        
    return M_acf;

#%% Coordination number

def fcn_slow(r, g, rho, nbins, L):
    
    fCN = np.zeros((nbins));
    dr = L/(2*nbins);
    
    for i in range(0, nbins):
        
        tsum = 0;
        
        for j in range(0,i):
            
             tsum += g[j]*r[j]**2*dr;
             
        tsum = tsum*4*np.pi*rho;
        fCN[i] = tsum; 
    
    return fCN;

#%% Generate n^3 integer multiples of q vectors

def gen_q(q_mag, L):
        
    # Get the spacing
    dk = (2*np.pi)/L;
    
    # Choose the upper limit of n
    n = np.int(np.ceil(q_mag/dk));
    ct = 0;
    qv = np.zeros(((2*n+1)**3,4));
    
    # Generate q vectors
    
    for j in range(-n, n+1):
        
        for k in range(-n, n+1):
            
            for l in range(-n, n+1):
                
                mag = (j**2 + k**2 + l**2)**(1/2);
                
                
                if (q_mag <= mag*dk < q_mag + dk):
                                    
                    qv[ct] = [j, k, l, mag*dk];
                    ct += 1;
    
    # Get rid of the zeros at the end of the array
    qv = qv[:ct];
                   
    # Sort the q_vectors by magnitude
    ind = np.argsort(qv[:,-1]);
    qv = qv[ind];
    
    # If the size is greater than 300, pick 300
    
    sz = qv.shape[0];
    
    if (sz > 60):
        
        ind = np.random.choice(sz, 60, replace=False);
        
        q_vec = np.zeros((60,3));
        q_vec = qv[ind,0:3]*dk;
        
    else:
        
        qv = qv*dk;        
        q_vec = qv[:,0:3];
        
    num_q = q_vec.shape[0];

    return q_vec, num_q;
