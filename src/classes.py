#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 14:04:46 2020

Class definitions for RDF, SSF, MSD, and VH functions that are used elsewhere.
Class definitions call functions from cfunctions, a shared library written in C and wrapped using SWIG.

author: marshalltekell
email: mct2180@columbia.edu

"""

# Import packages
import numpy as np #python array manipulation package

# Import shared libraries
import cfunctions
from analysis_functions import update_progress

#%% Pair distribution function class

class RDF:
    
    # Instance of class is created with number of particles, length of box, number of bins, instantaneous position data, density, and mass of specific particle
    def  __init__(self, N, L, nbins, r_int, rho, n_t):
        
        self.N = N;
        self.L = L;
        self.nbins = nbins;
        self.r_int = r_int;
        self.rho = rho;
        self.n_t = n_t;
        
    # g(r) (can be used for any particle of the same type or entire system)
    def RDF1(self):
        
        # Initialize arrays
        g_t = np.zeros((self.nbins), dtype=np.double);
        r1D = np.zeros((self.N*3), dtype=np.double);
        
        print('Calculating g(r) (particles of same type).')
        
        for i in range(0,self.n_t):
            
            update_progress(i/self.n_t);
            
            # Initialize arrays
            g = np.zeros((self.nbins), dtype=np.double);
            
            # Get positions for this time step
            r = self.r_int[:,i*3:(i+1)*3];
            
            
            # Make 2D array 1D
            for i in range(0,self.N):
                
                for j in range(0,3):
                    
                    r1D[i*3+j] = r[i,j];
            
            # Calculate g(r) (updates g in place)
            cfunctions.RDF1(g, r1D, self.nbins, self.N, self.rho, self.L);
            g_t += g;
            
        
        # Find average
        g_t = g_t/self.n_t;
        
        print('Done calculating g(r) (particles of same type).')
        
        return g_t;
    
    # Polymer g(r) (intra) (receives the number of beads per chain)
    def RDF2(self, Nc):
        
        # Initialize arrays
        g_t = np.zeros((self.nbins), dtype=np.double);
        r1D = np.zeros((self.N*3), dtype=np.double);
        
        print('Calculating g(r) (polymer intra).')
        
        for i in range(0,self.n_t):
            
            update_progress(i/self.n_t);
            
            # Initialize arrays
            g = np.zeros((self.nbins), dtype=np.double);
            
            # Get positions for this time step
            r = self.r_int[:,i*3:(i+1)*3];
            
            
            # Make 2D array 1D
            for i in range(0,self.N):
                
                for j in range(0,3):
                    
                    r1D[i*3+j] = r[i,j];
            
            # Calculate g(r) (updates g in place)
            cfunctions.RDF2(g, r1D, self.nbins, self.N, self.rho, self.L, Nc);
            g_t += g;
        
        # Find average
        g_t = g_t/self.n_t;
        
        print('Done calculating g(r) (polymer intra).')
        
        return g_t;
    
    # Polymer g(r) (inter) (receives the number of chains)
    def RDF3(self, Nc):
        
        # Initialize arrays
        g_t = np.zeros((self.nbins), dtype=np.double);
        r1D = np.zeros((self.N*3), dtype=np.double);
        
        print('Calculating g(r) (polymer inter).')
        
        for i in range(0,self.n_t):
            
            update_progress(i/self.n_t);
            
            # Initialize arrays
            g = np.zeros((self.nbins), dtype=np.double);
            
            # Get positions for this time step
            r = self.r_int[:,i*3:(i+1)*3];
            
            
            # Make 2D array 1D
            for i in range(0,self.N):
                
                for j in range(0,3):
                    
                    r1D[i*3+j] = r[i,j];
            
            # Calculate g(r) (updates g in place)
            cfunctions.RDF3(g, r1D, self.nbins, self.N, self.rho, self.L, Nc);
            g_t += g;
        
        # Find average
        g_t = g_t/self.n_t;
        
        print('Done calculating g(r) (polymer inter).')
        
        return g_t;
        
    # Beads of different types (concatenated into the same initial position array) 
    def RDF4(self, N2):
        
        # Initialize arrays
        g_t = np.zeros((self.nbins), dtype=np.double);
        r1D = np.zeros(((self.N+N2)*3), dtype=np.double);
        
        print('Calculating g(r) (particles of different type).')
        
        for i in range(0,self.n_t):
            
            update_progress(i/self.n_t);
            
            # Initialize arrays
            g = np.zeros((self.nbins), dtype=np.double);
            
            # Get positions for this time step
            r = self.r_int[:,i*3:(i+1)*3];
            
            
            # Make 2D array 1D
            for i in range(0,self.N+N2):
                
                for j in range(0,3):
                    
                    r1D[i*3+j] = r[i,j];
            
            # Calculate g(r) (updates g in place)
            cfunctions.RDF4(g, r1D, self.nbins, self.N, self.rho, self.L, N2);
            g_t += g;
        
        # Find average
        g_t = g_t/self.n_t;
        
        print('Done calculating g(r) (particles of different type).')
        
        return g_t;
    
#%% Static structure factor class
        
class SSF:
    
    def __init__(self, x, y, z, q_run, num_q, N, n_t, num_p):
        
        self.x = x;
        self.y = y;
        self.z = z;
        self.q_run = q_run;
        self.num_q = num_q;
        self.N = N;
        self.n_t = n_t;
        self.num_p = num_p;
        
        
    # Calculate the static structure factor for atoms of the same type
    def SSF1(self):
        
        # Initialize array
        Sq = np.zeros((self.num_p), dtype=np.double);
        
        print('Calculating S(q) (particles of same type).')                
        cfunctions.SSF1(Sq, self.x, self.y, self.z, self.q_run, self.num_q, self.N, self.n_t, self.num_p);
        print('Done with calculation.')
        
        return Sq;
    
    # Calculate the partial static structure factor for atoms of different types.
    
    def SSF2(self, N2):
        
        # Initialize array
        Sq = np.zeros((self.num_p), dtype=np.double);
                
        print('Calculating S(q) (particles of different types).')
        cfunctions.SSF2(Sq, self.x, self.y, self.z, self.q_run, self.num_q, self.N, N2, self.n_t, self.num_p)
        print('Done with calculation.')
        
        return Sq;
    
    
#%% van Hove function class

class vanHove:
    
    def __init__(self, x, y, z, N1, N2, nbins, n_t, L):
        
        self.x = x;
        self.y = y;
        self.z = z;
        self.N1 = N1;
        self.N2 = N2;
        self.nbins = nbins;
        self.n_t = n_t;
        self.L= L;
    
    # Calculate the van Hove function for linear time sampling    
    def VH0(self, m0, rhop):
        
        # Initialize arrays and variables
        dg = self.L/(2*(self.nbins-1));
        grt = np.zeros((self.nbins,self.n_t), dtype=np.double);
        grt_t = np.zeros((self.nbins*self.n_t), dtype=np.double);
        
        # Create count array
        ct = np.arange(1,self.n_t+1,1, dtype=np.intc);
        ct = ct[::-1];
        ct[0] = self.n_t-1;
        ct = np.ascontiguousarray(ct, dtype=np.intc);
                
        print('Calculating van Hove function with linear time sampling.')
        print('Printing the sample number.')
        
        # Get van Hove function (missing prefactors)
        cfunctions.vanHove(grt_t, ct, self.x, self.y, self.z, self.N1, self.N2, self.nbins, self.n_t, self.L);
        
        print('Done calculating van Hove function with linear time sampling.')
        
        # Use van Hove function to get g(r,t)
        for i in range(0, self.n_t):
            
            for j in range(0, self.nbins):
                
                dV = (1.333)*np.pi*((j+1)**3-j**3)*dg**3;
                grt[j,i] = (grt_t[i*self.nbins+j]*m0*10**24)/(self.n_t*self.N2*rhop*dV*6.022*10**(23));
                
              
        return grt;
    
    # Calculate the van Hove function for multiple loops with a log time sampling
    def VH1(self, V, numloop):
                
        # Initialize arrays and variables
        dg = self.L/(2*(self.nbins-1));
        grt = np.zeros((self.nbins,self.n_t), dtype=np.double);
        grt_t = np.zeros((self.nbins*self.n_t), dtype=np.double);
        grt_t2 = np.zeros((self.nbins*self.n_t), dtype=np.double);

        print('Calculating van Hove function with log time sampling.')
        
        update_progress((0)/numloop);
        
        # Loop over the number of runs
        for i in range(0, numloop):
            
            # Initialize array
            grt_t = np.zeros((self.nbins*self.n_t), dtype=np.double);
            
            # Get positions for this particular run
            x1D = self.x[i*(self.N1 + self.N2)*self.n_t:(i+1)*(self.N1 + self.N2)*self.n_t];
            y1D = self.y[i*(self.N1 + self.N2)*self.n_t:(i+1)*(self.N1 + self.N2)*self.n_t];
            z1D = self.z[i*(self.N1 + self.N2)*self.n_t:(i+1)*(self.N1 + self.N2)*self.n_t];
            
            # Perform calculation
            cfunctions.logvanHove(grt_t, x1D, y1D, z1D, self.N1, self.N2, self.nbins, self.n_t, self.L);
            
            # Add total
            grt_t2 += grt_t;
            
            update_progress((i+1)/numloop);
        
        #Take average
        grt_t2 = grt_t2/numloop;

        # Use van Hove function to get g(r,t)
        for i in range(0, self.n_t):
            
            for j in range(0, self.nbins):
                
                dV = (1.333)*np.pi*((j+1)**3-j**3)*dg**3;
                grt[j,i] = (grt_t2[i*self.nbins+j]*V)/(self.N1*self.N2*dV);
                
        return grt;
    
    # Calculate the van Hove function for linear time sampling    
    def VH2(self, V, n_t2):
        
        # Initialize arrays and variables
        dg = self.L/(2*(self.nbins-1));
        grt = np.zeros((self.nbins,self.n_t), dtype=np.double);
        grt_t = np.zeros((self.nbins*self.n_t), dtype=np.double);
                        
        print('Calculating van Hove function with linear time sampling.')
        
        # Get van Hove function (missing prefactors)
        cfunctions.vanHove_short(grt_t, self.x, self.y, self.z, self.N1, self.N2, self.nbins, self.n_t, n_t2, self.L);
        
        print('Done calculating van Hove function with linear time sampling.')
        
        # Use van Hove function to get g(r,t)
        for i in range(0, self.n_t):
            
            for j in range(0, self.nbins):
                
                dV = (1.333)*np.pi*((j+1)**3-j**3)*dg**3;
                grt[j,i] = (grt_t[i*self.nbins+j]*V)/(self.N1*self.N2*dV);
        
        return grt;
    
#%% Calculate mean squared and mean quartic displacement

class MSD:
    
    def __init__(self, x, y, z, N, n_t):
        
        self.x = x;
        self.y = y;
        self.z = z;
        self.N = N;
        self.n_t = n_t;

    # Calculate the mean squared displacement    
    def MSD0(self):
        
        # Initialize arrays
        r = np.zeros((self.n_t), dtype=np.double);
        ct = np.zeros((self.n_t), dtype=np.intc);
        x1D = np.zeros((self.N*self.n_t), dtype=np.double);
        y1D = np.zeros((self.N*self.n_t), dtype=np.double);
        z1D = np.zeros((self.N*self.n_t), dtype=np.double);
        
        # Make 2D position arrays 1D
        for i in range(0, self.N):
    
            for j in range(0, self.n_t):
                
                x1D[i*self.n_t + j] = self.x[i,j]
                y1D[i*self.n_t + j] = self.y[i,j]
                z1D[i*self.n_t + j] = self.z[i,j]
        
        # Calculate MSD (updates arrays r and ct in place)
        print('Calculating mean squared displacement.')
        cfunctions.calcMSD(r, ct, x1D, y1D, z1D, self.N, self.n_t)
        r = np.divide(r, ct, out=np.zeros_like(r), where=ct!=0);
        print('Done calculating mean squared displacement.')
        
        return r;
    
    # Calculate the mean quartic displacement
    def MSD1(self):
        
        # Initialize arrays
        r = np.zeros((self.n_t), dtype=np.double);
        ct = np.zeros((self.n_t), dtype=np.intc);
        x1D = np.zeros((self.N*self.n_t), dtype=np.double);
        y1D = np.zeros((self.N*self.n_t), dtype=np.double);
        z1D = np.zeros((self.N*self.n_t), dtype=np.double);
        
        # Make 2D position arrays 1D
        for i in range(0, self.N):
    
            for j in range(0, self.n_t):
                
                x1D[i*self.n_t + j] = self.x[i,j]
                y1D[i*self.n_t + j] = self.y[i,j]
                z1D[i*self.n_t + j] = self.z[i,j]
        
        # Calculate MSD (updates arrays r and ct in place)
        print('Calculating mean quartic displacement.')
        cfunctions.calcM4D(r, ct, x1D, y1D, z1D, self.N, self.n_t)
        r = np.divide(r, ct, out=np.zeros_like(r), where=ct!=0);
        print('Done calculating mean quartic displacement.')
        
        return r;
    
    # Calculate the mean squared displacement for looped log time sampling
    def MSD2(self, numloop):
        
        print("Calculating mean squared displacement.")
        # Initialize arrays
        x1D = np.zeros((self.N*self.n_t), dtype=np.double);
        y1D = np.zeros((self.N*self.n_t), dtype=np.double);
        z1D = np.zeros((self.N*self.n_t), dtype=np.double);
        r_t = np.zeros((self.n_t), dtype=np.double);
        
        for i in range(0, numloop):
            
            r = np.zeros((self.n_t), dtype=np.double);
            
            # Get data for current run
            for j in range(0, self.N):
                
                for k in range(0, self.n_t):
                    
                    x1D[j*self.n_t + k] = self.x[i*self.N+j, k];
                    y1D[j*self.n_t + k] = self.y[i*self.N+j, k];
                    z1D[j*self.n_t + k] = self.z[i*self.N+j, k];
                
            # Calculate MSD (square)
            cfunctions.logMSD(r, x1D, y1D, z1D, self.N, self.n_t)
            r_t += r/self.N;
                                
        # Find the average
        r_t = r_t/numloop;
        print('Done calculating mean squared displacement.')
        
        return r_t;
    
    # Calculate the mean quartic displacement for looped log time sampling
    def MSD3(self, numloop):
        
        print("Calculating mean quartic displacement.")
        # Initialize arrays
        x1D = np.zeros((self.N*self.n_t), dtype=np.double);
        y1D = np.zeros((self.N*self.n_t), dtype=np.double);
        z1D = np.zeros((self.N*self.n_t), dtype=np.double);
        r_t = np.zeros((self.n_t), dtype=np.double);
        
        for i in range(0, numloop):
            
            r = np.zeros((self.n_t), dtype=np.double);
            
            # Get data for current run
            for j in range(0, self.N):
                
                for k in range(0, self.n_t):
                    
                    x1D[j*self.n_t + k] = self.x[i*self.N+j, k];
                    y1D[j*self.n_t + k] = self.y[i*self.N+j, k];
                    z1D[j*self.n_t + k] = self.z[i*self.N+j, k];
                
            # Calculate MSD (square)
            cfunctions.logM4D(r, x1D, y1D, z1D, self.N, self.n_t)
            r_t += r/self.N;
                                
        # Find the average
        r_t = r_t/numloop;
        print('Done calculating mean quartic displacement.')
        
        return r_t;
    
