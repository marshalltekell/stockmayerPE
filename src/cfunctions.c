/***********************************************************************************************

INFORMATION      

author: Marshall Tekell
email: mct2180@columbia.edu

These are a set of functions that are intended to be called from Python using SWIG.
As such, each function returns void and, instead, operates on arrays in-place.
This means that the array is given to the function in the Python script that calls it.
These are designed to be used with one-dimensional numpy arrays and the specific syntax
is shown in the accompanying cfunctions.i file. So, each array is given along with its 1D size according to the syntax given in cfunctions.i.

These functions are intended to be used for computational efficiency but they are also useful illustrations of algorithms.
                
***********************************************************************************************/

#include "cfunctions.h"
#include  <stdlib.h>
#include  <stdio.h>
#include  <math.h>
#define PI 3.14159

/***********************************************************************************************

MEAN-SQUARED AND MEAN-QUARTIC DISPLACEMENT        

calcMSD averages over multiple time origins
logMSD averages over a single time origin
calcM4D calculates the mean-quartic displacement at multiple time origins
logM4D is mean-quartic displacement at a single time origin
                
***********************************************************************************************/

/* Calculate the mean squared displacement for linear time sampling (multipe time origin) */
void calcMSD(double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t) {

    int i, j, k;
    double mag;
    double d[3];
    
    for (i=0; i<N; i++){
    
        for (j=0; j<n_t-1; j++){
    
            for (k=j; k<n_t; k++){
        
                ct[k-j]++;
                mag = 0;
                d[0] = x[i*n_t+k] - x[i*n_t+j];
                d[1] = y[i*n_t+k] - y[i*n_t+j];
                d[2] = z[i*n_t+k] - z[i*n_t+j];
                mag = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                r[k-j] += mag;
            
            }
        }
    }
}
    
/* Calculate the mean  quartic displacement for linear time sampling (multipe time origin) */
void calcM4D(double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t) {

    int i, j, k;
    double mag;
    double d[3];
    
    for (i=0; i<N; i++){
    
        for (j=0; j<n_t-1; j++){
    
            for (k=j; k<n_t; k++){
        
                ct[k-j]++;
                mag = 0;
                d[0] = x[i*n_t+k] - x[i*n_t+j];
                d[1] = y[i*n_t+k] - y[i*n_t+j];
                d[2] = z[i*n_t+k] - z[i*n_t+j];
                mag = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
                mag = pow(mag,2);
                r[k-j] += mag;
            
            }
        }
    }
}
        
/* Calculate the mean squared displacement for log time sampling (single time origin) */
void logMSD(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t) {

    int i, j;
    double mag;
    double d[3];
    
    for (i=0; i<N; i++){
    
        for (j=0; j<n_t; j++){
            
            mag = 0;
            d[0] = x[i*n_t+j] - x[i*n_t];
            d[1] = y[i*n_t+j] - y[i*n_t];
            d[2] = z[i*n_t+j] - z[i*n_t];
            mag = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            r[j] += mag;
            
        }
    }
}
            
/* Calculate the mean squared displacement for log time sampling (single time origin) */
void logM4D(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t) {

    int i, j;
    double mag;
    double d[3];
    
    for (i=0; i<N; i++){
    
        for (j=0; j<n_t; j++){
            
            mag = 0;
            d[0] = x[i*n_t+j] - x[i*n_t];
            d[1] = y[i*n_t+j] - y[i*n_t];
            d[2] = z[i*n_t+j] - z[i*n_t];
            mag = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            mag = pow(mag,2);
            r[j] += mag;
            
        }
    }
}
            
/***********************************************************************************************

RADIAL DISTRIBUTION FUNCTION              
                            
RDF1 considers the whole system
RDF2 considers pairs on the same molecule (assumes M atoms on each molecule)
RDF2 considers pairs on the different molecules (assumes M atoms on each molecule)
RDF1 considers two particles of different types
fCN is just a numerical integration with density normalization to get the coordination number.
                
***********************************************************************************************/
            
/* Calculate the pair distribution function for particles of the same type (or entire system) */
void RDF1(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L) {

    int i, j, k, l;
    double dg = L/(2*(nbins-1));
    double d[3], mag, dV;
    
    for (i=0; i<N-1; i++){
        
        for (j=i+1; j<N; j++){
            
            mag = 0;
            
            for (k=0; k<3; k++){
                
                d[k] = r[i*3+k] - r[j*3+k];
                d[k] = d[k]-L*round(d[k]/L);
                mag = mag + d[k]*d[k];
            }
             
            mag = sqrt(mag);
            
            if (mag < L/2){
                
                l = lrint(mag/dg);
                g[l] = g[l] + 2;
            }           
        }
    }
           
    for (i=0; i<nbins; i++){
            
        dV = 1.333*PI*(pow((i+1),3)-pow(i,3))*pow(dg,3);
        g[i] = (g[i])/(dV*N*rho);
    
    }
}

/* Calculate the pair distribution function for polymer beads on the same chain */
void RDF2(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L, int Nc) {

    int i, j, k, l, mol_i, mol_j;
    double dg = L/(2*(nbins-1));
    double d[3], mag, dV;
    
    for (i=0; i<N-1; i++){
        
        for (j=i+1; j<N; j++){
    
            /* Check to see if on the same chain */
            mol_i = floor(i/Nc);
            mol_j = floor(j/Nc);
            
            mag = 0;
    
            if (mol_i == mol_j){
            
                for (k=0; k<3; k++){
                    
                    d[k] = r[i*3+k] - r[j*3+k];
                    d[k] = d[k]-L*round(d[k]/L);
                    mag = mag + d[k]*d[k];
                }
                 
                mag = sqrt(mag);
                
                if (mag < L/2){
                    
                    l = lrint(mag/dg);
                    g[l] = g[l] + 2;
                } 
            }          
        }
    }
           
    for (i=0; i<nbins; i++){
            
        dV = 1.333*PI*(pow((i+1),3)-pow(i,3))*pow(dg,3);
        g[i] = (g[i])/(dV*N*rho);
    
    }
}

/* Calculate the pair distribution function for polymer beads NOT on the same chain */
void RDF3(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L, int Nc) {

    int i, j, k, l, mol_i, mol_j;
    double dg = L/(2*(nbins-1));
    double d[3], mag, dV;
    
    for (i=0; i<N-1; i++){
        
        for (j=i+1; j<N; j++){
    
            /* Check to see if on the same chain */
            mol_i = floor(i/Nc);
            mol_j = floor(j/Nc);
            
            mag = 0;
    
            if (mol_i != mol_j){
            
                for (k=0; k<3; k++){
                    
                    d[k] = r[i*3+k] - r[j*3+k];
                    d[k] = d[k]-L*round(d[k]/L);
                    mag = mag + d[k]*d[k];
                }
                 
                mag = sqrt(mag);
                
                if (mag < L/2){
                    
                    l = lrint(mag/dg);
                    g[l] = g[l] + 2;
                } 
            }          
        }
    }
           
    for (i=0; i<nbins; i++){
            
        dV = 1.333*PI*(pow((i+1),3)-pow(i,3))*pow(dg,3);
        g[i] = (g[i])/(dV*N*rho);
    
    }
}
    
    
/* Calculate the pair distribution function for particles of different type */
void RDF4(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L, int N2) {

    int i, j, k, l;
    double dg = L/(2*(nbins-1));
    double d[3], mag, dV;
    
    for (i=0; i<N; i++){
        
        for (j=N; j<N+N2; j++){
            
            mag = 0;
            
            for (k=0; k<3; k++){
                
                d[k] = r[i*3+k] - r[j*3+k];
                d[k] = d[k]-L*round(d[k]/L);
                mag = mag + d[k]*d[k];
            }
             
            mag = sqrt(mag);
            
            if (mag < L/2){
                
                l = lrint(mag/dg);
                g[l] = g[l] + 1;
            }           
        }
    }
           
    for (i=0; i<nbins; i++){
            
        dV = 1.333*PI*(pow((i+1),3)-pow(i,3))*pow(dg,3);
        g[i] = (g[i])/(dV*N*rho);
    
    }
}
           
/* Calculate the coordination number given a g(r) */
void fcn(double* f, int nbins3, double* r_RDF, int nbins2, double* g, int nbins1, int nbins, double L, double rho){
        
    double dr = L/(2*nbins);
    int i, j;
    double t1;

    for (i=0; i<nbins; i++){

        t1 = 0.0;
    
        for (j=0; j<i; j++){
    
            t1 += g[j]*r_RDF[j]*r_RDF[j]*dr;

        }

        t1 = t1*4*PI*rho;
        f[i] = t1;
    }      
    
}
            
/***********************************************************************************************

IONIC CONDUCTIVITY               
                            
Calculates the ionic conductivity for a polymer electrolyte simulation.
lam calculates the whole Einstein expression with multiple time origins.
loglam calculates the whole Einstein expression with a single time origins.
(log)lamcross1 calculate the expression for i!=j.
(log)lamcross2 calculate the expression for i!=j, zi==zj.
(log)lamcross2 calculate the expression for i!=j, zi!=zj.
                
***********************************************************************************************/
        
/* Calculate the ionic conductivity */
void lam(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k, l;
    double mag;
    double di[3];
    double dj[3];
    int t1;
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
    
            t1 = (i*N + j)/(N*N);
    
            printf("\rProgress: %3.5f%%", 0.0025*(i*N + j));
            /*fflush(stdout);*/
            
            /* Loop over time origins */
            for (k=0; k<n_t-1; k++){
        
                /* Loop over time origins */
                for (l=k; l<n_t; l++){
                    
                    /* Reset magnitude */
                    mag = 0;
        
                    /* Get MSD */
                    di[0] = x[i*n_t+l] - x[i*n_t+k];
                    di[1] = y[i*n_t+l] - y[i*n_t+k];
                    di[2] = z[i*n_t+l] - z[i*n_t+k];
                    dj[0] = x[j*n_t+l] - x[j*n_t+k];
                    dj[1] = y[j*n_t+l] - y[j*n_t+k];
                    dj[2] = z[j*n_t+l] - z[j*n_t+k];
                    mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                    
                    /* Assign to sum */
                    r[l-k] += mag;

                }       
            }
            
            for (k=0; k<n_t; k++){
                    
                /* Take average of the number of time origins */
                r[k] = r[k]/ct[k];
                rt[k] += q[i]*q[j]*r[k];
                r[k] = 0.0;
                    
            }
        }
    }
}
            
/* Calculate the ionic conductivity (single time origin) */
void loglam(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k;
    double mag;
    double di[3];
    double dj[3];
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
                
            /* Loop over time origins */
            for (k=0; k<n_t; k++){
                    
                /* Reset magnitude */
                mag = 0;
    
                /* Get MSD */
                di[0] = x[i*n_t+k] - x[i*n_t];
                di[1] = y[i*n_t+k] - y[i*n_t];
                di[2] = z[i*n_t+k] - z[i*n_t];
                dj[0] = x[j*n_t+k] - x[j*n_t];
                dj[1] = y[j*n_t+k] - y[j*n_t];
                dj[2] = z[j*n_t+k] - z[j*n_t];
                mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                
                /* Assign to sum */
                r[k] += q[i]*q[j]*mag;
     
            }
        }
    }
}
            
/* Calculate the ionic conductivity--just the cross terms*/
void lam_cross1(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k, l;
    double mag;
    double di[3];
    double dj[3];
    int t1;
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
    
            t1 = (i*N + j)/(N*N);
    
            printf("\rProgress: %3.5f%%", 0.0025*(i*N + j));
            /*fflush(stdout);*/
    
            if (i!=j){
                
                /* Loop over time origins */
                for (k=0; k<n_t-1; k++){
            
                    /* Loop over time origins */
                    for (l=k; l<n_t; l++){
                        
                        /* Reset magnitude */
                        mag = 0;
            
                        /* Get MSD */
                        di[0] = x[i*n_t+l] - x[i*n_t+k];
                        di[1] = y[i*n_t+l] - y[i*n_t+k];
                        di[2] = z[i*n_t+l] - z[i*n_t+k];
                        dj[0] = x[j*n_t+l] - x[j*n_t+k];
                        dj[1] = y[j*n_t+l] - y[j*n_t+k];
                        dj[2] = z[j*n_t+l] - z[j*n_t+k];
                        mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                        
                        /* Assign to sum */
                        r[l-k] += mag;
    
                    }       
                }
            
                for (k=0; k<n_t; k++){
                        
                    /* Take average of the number of time origins */
                    r[k] = r[k]/ct[k];
                    rt[k] += q[i]*q[j]*r[k];
                    r[k] = 0.0;
                        
                }
            }
        }
    }
}
                        
/* Calculate the ionic conductivity--just the cross terms, single time origin*/
void loglam_cross1(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k;
    double mag;
    double di[3];
    double dj[3];
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
        
            if (i!=j){

                /* Loop over time origins */
                for (k=0; k<n_t; k++){
                                    
                    /* Reset magnitude */
                    mag = 0.0;
        
                    /* Get MSD */
                    di[0] = x[i*n_t+k] - x[i*n_t];
                    di[1] = y[i*n_t+k] - y[i*n_t];
                    di[2] = z[i*n_t+k] - z[i*n_t];
                    dj[0] = x[j*n_t+k] - x[j*n_t];
                    dj[1] = y[j*n_t+k] - y[j*n_t];
                    dj[2] = z[j*n_t+k] - z[j*n_t];
                    mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                    
                    /* Assign to sum */
                    r[k] += q[i]*q[j]*mag;
                }
            }
        }
    }
}
                
/* Calculate the ionic conductivity--just the cross terms where zi = zj*/
void lam_cross2(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k, l;
    double mag;
    double di[3];
    double dj[3];
    int t1;
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
    
            t1 = (i*N + j)/(N*N);
    
            printf("\rProgress: %3.5f%%", 0.0025*(i*N + j));
            /*fflush(stdout);*/
    
            if (i!=j && q[i]==q[j]){
                
                /* Loop over time origins */
                for (k=0; k<n_t-1; k++){
            
                    /* Loop over time origins */
                    for (l=k; l<n_t; l++){
                        
                        /* Reset magnitude */
                        mag = 0;
            
                        /* Get MSD */
                        di[0] = x[i*n_t+l] - x[i*n_t+k];
                        di[1] = y[i*n_t+l] - y[i*n_t+k];
                        di[2] = z[i*n_t+l] - z[i*n_t+k];
                        dj[0] = x[j*n_t+l] - x[j*n_t+k];
                        dj[1] = y[j*n_t+l] - y[j*n_t+k];
                        dj[2] = z[j*n_t+l] - z[j*n_t+k];
                        mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                        
                        /* Assign to sum */
                        r[l-k] += mag;
    
                    }       
                }
            
                for (k=0; k<n_t; k++){
                        
                    /* Take average of the number of time origins */
                    r[k] = r[k]/ct[k];
                    rt[k] += q[i]*q[j]*r[k];
                    r[k] = 0.0;
                        
                }
            }
        }
    }
}
                        
/* Calculate the ionic conductivity--just the cross terms, single time origin*/
void loglam_cross2(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k;
    double mag;
    double di[3];
    double dj[3];
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
        
            if (i!=j && q[i]==q[j]){
                
                /* Loop over time origins */
                for (k=0; k<n_t; k++){
                                    
                    /* Reset magnitude */
                    mag = 0;
        
                    /* Get MSD */
                    di[0] = x[i*n_t+k] - x[i*n_t];
                    di[1] = y[i*n_t+k] - y[i*n_t];
                    di[2] = z[i*n_t+k] - z[i*n_t];
                    dj[0] = x[j*n_t+k] - x[j*n_t];
                    dj[1] = y[j*n_t+k] - y[j*n_t];
                    dj[2] = z[j*n_t+k] - z[j*n_t];
                    mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                    
                    /* Assign to sum */
                    r[k] += q[i]*q[j]*mag;
                }
            }
        }
    }
}
    
/* Calculate the ionic conductivity--just the cross terms where zi != zj*/
void lam_cross3(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k, l;
    double mag;
    double di[3];
    double dj[3];
    int t1;
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
    
            t1 = (i*N + j)/(N*N);
    
            printf("\rProgress: %3.5f%%", 0.0025*(i*N + j));
            /*fflush(stdout);*/
    
            if (i!=j && q[i]!=q[j]){
                
                /* Loop over time origins */
                for (k=0; k<n_t-1; k++){
            
                    /* Loop over time origins */
                    for (l=k; l<n_t; l++){
                        
                        /* Reset magnitude */
                        mag = 0;
            
                        /* Get MSD */
                        di[0] = x[i*n_t+l] - x[i*n_t+k];
                        di[1] = y[i*n_t+l] - y[i*n_t+k];
                        di[2] = z[i*n_t+l] - z[i*n_t+k];
                        dj[0] = x[j*n_t+l] - x[j*n_t+k];
                        dj[1] = y[j*n_t+l] - y[j*n_t+k];
                        dj[2] = z[j*n_t+l] - z[j*n_t+k];
                        mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                        
                        /* Assign to sum */
                        r[l-k] += mag;
    
                    }       
                }
            
                for (k=0; k<n_t; k++){
                        
                    /* Take average of the number of time origins */
                    r[k] = r[k]/ct[k];
                    rt[k] += q[i]*q[j]*r[k];
                    r[k] = 0.0;
                        
                }
            }
        }
    }
}
                        
/* Calculate the ionic conductivity--just the cross terms, single time origin*/
void loglam_cross3(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t) {

    int i, j, k;
    double mag;
    double di[3];
    double dj[3];
    
    /* Loop over particles */
    for (i =0; i<N; i++){
    
        /* Loop over particles */
        for (j=0; j<N; j++){
        
            if (i!=j && q[i]!=q[j]){
                
                /* Loop over time origins */
                for (k=0; k<n_t; k++){
                                    
                    /* Reset magnitude */
                    mag = 0;
        
                    /* Get MSD */
                    di[0] = x[i*n_t+k] - x[i*n_t];
                    di[1] = y[i*n_t+k] - y[i*n_t];
                    di[2] = z[i*n_t+k] - z[i*n_t];
                    dj[0] = x[j*n_t+k] - x[j*n_t];
                    dj[1] = y[j*n_t+k] - y[j*n_t];
                    dj[2] = z[j*n_t+k] - z[j*n_t];
                    mag = di[0]*dj[0] + di[1]*dj[1] + di[2]*dj[2];
                    
                    /* Assign to sum */
                    r[k] += q[i]*q[j]*mag;
                }
            }
        }
    }
}

/***********************************************************************************************

SELF-INTERMEDIATE SCATTERING FUNCTION               
                            
Calculates the self-intermediate scattering functions on log (logfqt) and lin (fqt) samplings.  static structure factor. The wave-vectors are generated elsewhere.
Logfqt loops over times from a single time origin, averaging across time origins happens elswhere.
Fqt loops over different time origins.
                
***********************************************************************************************/
                    
                    
/* Calculate the self-intermediate scattering function */       
void fqt(double* fqt, int n_t4, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int N, int n_t, int num_q){
                                    
    int i, j, k, l;
    double ts[3] = {0};
    double tsq, t1;
    
    /* Loop over particles and times again */   
    for (i=0; i<N; i++){
                    
        printf("\r Progress: %.4f %%",100.0*((1.0*i)/(1.0*N)));
        fflush(stdout); 
        
        /* Loop over time origins */
        for (j=0; j<n_t-1; j++){
                                
            /* For particular time, loop over all successive times */
            for (k=j; k<n_t; k++){
                
                /* Get displacement */    
                ts[0] = x[i*n_t+k] - x[i*n_t+j];
                ts[1] = y[i*n_t+k] - y[i*n_t+j];
                ts[2] = z[i*n_t+k] - z[i*n_t+j];
                   
                /* Reset sum */     
                tsq = 0.0;
                
                for (l=0; l<num_q; l++){
                    
                    t1 = qv[3*l]*ts[0] + qv[3*l+1]*ts[1] + qv[3*l+2]*ts[2];
                    tsq += cos(t1);
    
                }

            /* Assign value to F(q,t) */
            fqt[k-j] += tsq;
                                    
            }                            
        }  
    } 
                    
    /* Take the average of the number of time origins, wave vectors, and particles */
    for (k=0; k<n_t; k++){
        
        fqt[k] = fqt[k]/(ct[k]*num_q*N);
    
    }
                                                       
}
            
/* Calculate the self-intermediate scattering function--single time origin */       
void logfqt(double* fqt, int n_t4, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int N, int n_t, int num_q){
                                    
    int i, j, k;
    double ts[3] = {0};
    double tsq, t1;
    
    /* Loop over particles and times */   
    for (i=0; i<N; i++){
        
        /* Loop over time origins */
        for (j=0; j<n_t; j++){
            
            /* Get displacement */                                        
            ts[0] = x[i*n_t+j] - x[i*n_t];
            ts[1] = y[i*n_t+j] - y[i*n_t];
            ts[2] = z[i*n_t+j] - z[i*n_t];
                   
            /* Reset sum */     
            tsq = 0.0;
                
            for (k=0; k<num_q; k++){
                
                t1 = qv[3*k]*ts[0] + qv[3*k+1]*ts[1] + qv[3*k+2]*ts[2];
                tsq += cos(t1);

            }
    
            /* if (j==0){printf("tsq: %.6f",tsq);} */

            /* Assign value to F(q,t) */
            fqt[j] += tsq/(num_q*N);
            
            /* if (j==0){printf("fqt: %.6f",fqt[j]);} */
                               
        }  
    }                                    
}
            
/***********************************************************************************************

STATIC STRUCTURE FACTOR                
                            
Functions to calculate the static structure factor. These use wave vectors generated elsewhere.
There are two functions: SSF1 and SSF2. SSF1 considers particles of the same type (fast). SSF2 considers particles of different types (slow).
                
***********************************************************************************************/
        
/* Calculate the static structure factor (same type) */
void SSF1(double* ssf, int n_k, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int* nq, int sz_nq, int N, int n_t, int num_p){
    
    double t1, t2, t3;
    int i, j, k, l;
                                                        
    /* Loop over the number of times */
    for (i=0; i<n_t; i++){

        printf("\rProgress: %3.5f%%", 100.0*((1.0*i)/(1.0*n_t)));
        /*fflush(stdout);*/                                                        
        
        /* Loop over the number of points (L/(2*dk)) */
        for (j=0; j<num_p; j++){
            
            /* Reset sums */
            t2 = 0.0;
            t3 = 0.0;
            
            /* Loop over the number of particles */
            for (k=0; k<N; k++){
    
                /* Loop over the number of q-vectors at each magnitude */
                for (l=0; l<nq[j]; l++){
                    
                    /* Perform calculation, see ref. Zhang arXiv */
                    t1 = x[k*n_t+i]*qv[3*num_p*l+3*j] + y[k*n_t+i]*qv[3*num_p*l+3*j+1] + z[k*n_t+i]*qv[3*num_p*l+3*j+2];
                    t2 += cos(t1);
                    t3 += sin(t1);
                    
                }
            }
                    
            /* Average over N, nq, and n_t */
            ssf[j] += (pow(t2,2.0) + pow(t3,2.0))/(nq[j]*N*n_t);     
                
        }      
    }                        
}                         

/* Calculate the partial static structure factor (different types) */
void SSF2(double* ssf, int n_k, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int* nq, int sz_nq, int N1, int N2, int n_t, int num_p){

    double t1, t2, t3;
    int i, j, k, l, m;
    double d[3] = {0};
                 
    /* Loop over the number of times */
    for (i=0; i<n_t; i++){
    
        printf("\rProgress: %3.5f%%", 100.0*((1.0*i)/(1.0*n_t)));
                    
        /* Loop over the number of points (L/(2*dk)) */
        for (j=0; j<num_p; j++){
            
            /* Reset sums */
            t2 = 0.0;
            t3 = 0.0;
            
            /* Loop over particles of the first type */
            for (k=0; k<N1; k++){

                /* Loop over particles of the second type */
                for (l=N1; l<N1+N2; l++){

                    /* Get the distance between particles */
                    d[0] = x[n_t*k+i] - z[n_t*l+i];
                    d[1] = y[n_t*k+i] - y[n_t*l+i];
                    d[2] = z[n_t*k+i] - z[n_t*l+i];
    
                    /* Loop over the number of q-vectors */
                    for (m=0; m<nq[j]; m++){
                 
                        t1 = d[0]*qv[3*num_p*m+3*j] + d[1]*qv[3*num_p*m+3*j+1] + d[2]*qv[3*num_p*m+3*j+2];
                        t2 += cos(t1);
                    
                    }
                }
            }
            
            /* Average over (N*I)**(1/2), nq, n_t */
            ssf[j] += t2/(n_t*nq[j]*pow(N1*N2,0.5));
            
        }
    }
}

    
/***********************************************************************************************

VAN HOVE FUNCTION               
                            
Calculate the van Hove function for linear (multiple time origins) (vanHove) and log (single time origin) (logvanHove) 
                
***********************************************************************************************/
    
/* Calculate the time-dependent van Hove function between the cation and the bead */
void vanHove(double* g, int nbins1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N1, int N2, int nbins, int n_t, double L){
    
    int t2 = 0;
    double mag = 0.0;
    double d[3] = {0};
    double dg = 1.0*(L/(2*(nbins-1)));
    int *g1 = (int *) malloc((nbins)*sizeof(int));
    int i, j, k, l, m, t1;
    
    for (i=0; i<n_t-1; i++){
        
        /*printf("\rProgress: %d", i); */
        
        for (j=i; j<n_t; j++){
        
            t2 += 1;
            
            printf("\rProgress: %3.5f%%", 100.0*((2.0*t2)/((n_t-1)*n_t)));
                            
            for (k=0; k<N1; k++){
                
                for (l=N1; l<N1+N2; l++){
                    
                    /* Get the distance between particles */
                    d[0] = x[n_t*k+i] - z[n_t*l+i];
                    d[1] = y[n_t*k+i] - y[n_t*l+i];
                    d[2] = z[n_t*k+i] - z[n_t*l+i];
                    
                    /* Reset magnitude */
                    mag = 0.0;
                    
                    /* Apply periodic boundary conditions */
                    for (m=0; m<3; m++){
                        
                        d[m] += -L*round(d[m]/L);
                        mag += pow(d[m],2);
                            
                    }
                                        
                    /* Get norm */
                    mag = pow(mag,0.5);
                    
                    if (mag <= L/2){ 
                        
                        /* Get bin for magnitude */
                        t1 = lrint(mag/dg); 
            
                        /* Add it to count for this time */
                        g1[t1] += 1; 
            
                    }
                }   
            }
    
            /* Give this result to overall average */
            for (k=0; k<nbins; k++){ 
                
                g[(j-i)*nbins+k] += g1[k]/ct[j-i];

                /* Reset g1 */
                g1[k] = 0;  
                
            }    
        }     
    }
    free(g1);            
}

/* Calculate the time-dependent van Hove function between the cation and the bead */
void logvanHove(double* g, int nbins1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N1, int N2, int nbins, int n_t, double L){
    
    /* Initialize arrays */
    double mag = 0.0;
    double d[3] = {0};
    double dg = 1.0*(L/(2*(nbins-1)));
    int *g1 = (int *) malloc((nbins)*sizeof(int));
    int i, j, k, l, t1;
    
     /* Initialize g1 */
    for (i=0; i<nbins; i++){
        
        g1[i] = 0;
        
    }  
        
    for (i=0; i<n_t; i++){
        
        for (j=0; j<N1; j++){
                
            /*printf("\rProgress: %3.5f%%", 100.0*((1.0*i*n_t+j)/(2.0*n_t*n_t))); */
            
            for (k=N1; k<N1+N2; k++){
            
                /* Get the distance between particles */
                d[0] = x[n_t*j] - x[n_t*k+i];
                d[1] = y[n_t*j] - y[n_t*k+i];
                d[2] = z[n_t*j] - z[n_t*k+i];
                
                /* Reset magnitude */
                mag = 0.0;
                
                /* Apply periodic boundary conditions */
                for (l=0; l<3; l++){
                    
                    d[l] = d[l]-L*round(d[l]/L);
                    mag += pow(d[l],2);
                        
                }
                                    
                /* Get norm */
                mag = pow(mag,0.5);
            
                if (mag < L/2){
                    
                    /* Get bin for magnitude */
                    t1 = lrint(mag/dg);
        
                    /* Add it to count for this time */
                    g1[t1] += 1;
        
                }
            
            }
        }
        
        /* Give this result to overall average */
        for (j=0; j<nbins; j++){
            
            g[i*nbins+j] += g1[j];
    
            /* Reset g1 */
            g1[j] = 0;
        }     
    }
    free(g1);                                    
}
    
/* Calculate the time-dependent van Hove function between the cation and the bead */
void vanHove_short(double* g, int nbins1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N1, int N2, int nbins, int n_t, int n_t2, double L){
    
    int t2 = 0;
    double mag = 0.0;
    double d[3] = {0};
    double dg = 1.0*(L/(2*(nbins-1)));
    int *g1 = (int *) malloc((nbins)*sizeof(int));
    int i, j, k, l, m, t1;
    
    /* Initialize g1 */
    for (i=0; i<nbins; i++){
        
        g1[i] = 0;
        
    } 
    
    for (i=0; i<n_t-n_t2; i++){
        
        /*printf("\rProgress: %d", i); */
        
        for (j=i; j<i+n_t2; j++){
        
            t2 += 1;
            
            printf("\rProgress: %3.5f%%", 100.0*((1.0*t2)/((n_t-n_t2)*n_t2)));
                            
            for (k=0; k<N1; k++){
                
                for (l=N1; l<N1+N2; l++){
                    
                    /* Get the distance between particles */
                    d[0] = x[n_t*k+i] - x[n_t*l+j];
                    d[1] = y[n_t*k+i] - y[n_t*l+j];
                    d[2] = z[n_t*k+i] - z[n_t*l+j];
                    
                    /* Reset magnitude */
                    mag = 0.0;
                    
                    /* Apply periodic boundary conditions */
                    for (m=0; m<3; m++){
                        
                        d[m] = d[m] - L*round(d[m]/L);
                        mag += pow(d[m],2);
                            
                    }
                                        
                    /* Get norm */
                    mag = pow(mag,0.5);
                    
                    if (mag <= L/2){ 
                        
                        /* Get bin for magnitude */
                        t1 = lrint(mag/dg); 
            
                        /* Add it to count for this time */
                        g1[t1] += 1; 
            
                    }
                }   
            }
    
            /* Give this result to overall average */
            for (k=0; k<nbins; k++){
                
                g[(j-i)*nbins+k] += 1.0*g1[k]*((1.0)/(n_t-n_t2));

                /* Reset g1 */
                g1[k] = 0;  
                
            }     
                
        }     
    }
                                                                   
    free(g1);            
}
                        
/***********************************************************************************************

DIPOLE-DIPOLE AUTOCORRELATION FUNCTION               
                            
Calculate the dipole-dipole autocorrelation function for linear (multiple time origins) (MUAC) and log (single time origin) (MUAClog) 
                
***********************************************************************************************/
    
/* Calculates the auto-correlation function for the dipole moment in the system */

void MUAC(double* muac, int n_t5, double* mu_int, int n_t6, int* ct, int n_t2, int N,  int n_t){

    double t1 = 0.0;
    int i, j, k, l;
                        
    /* Loop over the number of particles, the number of time origins, and the successive times */
    
    for (i=0; i<N; i++){
    
        for (j=0; j<n_t-1; j++){
    
            for (k=0; k<n_t; k++){
                        
                t1 = 0.0;
                        
                /* Take dot product */
                for (l=0; l<3; l++){

                    t1 += mu_int[i*n_t*3 + 3*j +l]*mu_int[i*n_t*3 + 3*k+l];

                }
    
                /* Assign to function */
                muac[k-j] += t1;
            }
        }                                                           
    }
     
    /* Divide by number of times and particles */
    for (i=0; i<n_t; i++){
        
        muac[i] = muac[i]/(N*ct[i]);  
    }   

    
}

/* Calculates the auto-correlation function for the dipole moment in the system--log time sampling, single time origin */

void MUAC_log(double* muac, int n_t5, double* mu_int, int n_t6, int N,  int n_t){
    

    double t1 = 0.0;
    int i, j, k;
    
    /* Loop over the number of particles, the number of time origins, and the successive times */
    
    for (i=0; i<N; i++){

        for (j=0; j<n_t; j++){
    
            t1 = 0.0;
    
            /* Take dot product */
            for (k=0;  k<3; k++){

                t1 += mu_int[i*n_t*3 + k]*mu_int[i*n_t*3 + 3*j + k];

            }
            
            /* Assign to function */
            muac[j] += t1/N;
        
        }       
    }                    
}
 
/***********************************************************************************************

Coordination number statistics               
                            
Calculate the number of monomers within a certain cutoff--lo < rij < hi--for a given collection of particles at a sampled timestep
    
Calculate the number of chains in a single coordination structure
                
***********************************************************************************************/
           
  
/* Calculates number of N2 beads within lo < rij < hi for N1 beads */           
          
void pofn(int* p1, int sz_p, double* r, int sz_r, int N1, int N2, double lo, double hi, double L){
     
    int i, j, k;                   
    double d[3], mag;
    
    /* Loop over the number of cations */ 
    for (i=0; i<N1; i++){
                
        /* Loop over the number of PEO monomers */
        for (j=N1; j<N1+N2; j++){
                
            /* Reset magnitude */
            mag = 0.0;
            
             /* Loop over the number of dimensions */
            for (k=0; k<3; k++){
    
                d[k] = r[i*3+k] - r[j*3+k];
                d[k] = d[k]-L*round(d[k]/L);
                mag = mag + d[k]*d[k];
            
            }
            
            mag = sqrt(mag);

            /* Test if distance is within cutoff */ 

            if (mag >= lo && mag < hi){
    
                p1[i] += 1;

            }   
        }                  
    }             
}

            
/* Calculates number of unique chains within lo < rij < hi for N1 beads */           
        
void pofM(int* p1, int sz_p, double* r, int sz_r, int N1, int N2, double lo, double hi, double L, int Nc){
     
    int i, j, k, mol_j, sw1, sw2, ct1;                   
    double d[3], mag;
    int ch[6];

    /* Initialize ch */
    for (i=0; i<6; i++){

        ch[i]= -1;
                 
    }
        
    /* Loop over the number of cations */ 
    for (i=0; i<N1; i++){

        /* Switches */
        sw1 = 0;
        ct1 = 0;
                
        /* Loop over the number of PEO monomers */
        for (j=N1; j<N1+N2; j++){
           
            /* Reset second switch */
            sw2 = 0;

            /* Get the chain number for the jth particle */
            mol_j = floor((j-N1)/Nc);
                
            /* Reset magnitude */
            mag = 0.0;
            
             /* Loop over the number of dimensions */
            for (k=0; k<3; k++){
    
                d[k] = r[i*3+k] - r[j*3+k];
                d[k] = d[k]-L*round(d[k]/L);
                mag = mag + d[k]*d[k];
            
            }
            
            mag = sqrt(mag);

            /* Test if distance is within cutoff */ 
            if (mag >= lo && mag < hi){
        
                /* For first success store chain number */
                if (sw1 ==0){
                    
                    ch[0] = mol_j;
                    sw1 = 1;
                    ct1 += 1;
                    
                }
                
                /* For subsequent chain numbers, compare this chain to previous list */
                if (sw1 ==1){
                    
                    /* Check to see if you've already counted this chain */
                    for (k=0; k<6; k++){
                        
                        if (mol_j==ch[k]){
            
                            sw2 = 1;
                        } 
                    }
            
                    /* Only if new chain, add*/
                    if (sw2==0){

                        ch[ct1] = mol_j;
                        ct1 += 1;

                    }
                }
            }   
        }
        
        /* Store ct1 (the number of unique chains) */
        p1[i] = ct1;
                          
    }             
}

/***********************************************************************************************

Sum pair-wise interactions

***********************************************************************************************/

/* Calculates pairwise interactions between polymer beads in electrolyte */

void pair_lj_dipole(double* u, int sz_u, double* r, int sz_r, double* mu, int sz_mu, double ep, double sig, double L, int N) {

    int i, j, k;
    double mag, u_t, dot1, dot2, dot3, t1, t2;
    double d[3];
    
    /* Initialize internal energy for this timestep*/
    u_t = 0.0;
    
    /* Loop over the particles (int N = number of particles) */
    for (i=0; i<N; i++){
        
        /* Loop over the upper half of the matrix */
        for (j=i+1; j<N; j++){
            
            /* Initialize magnitude and energy*/
            mag = 0.0;
        
            /* Get distance between the particles */
            for (k=0; k<3; k++) {
            
                d[k] = r[i*3 + k] - r[j*3 + k];
                d[k] = d[k]-L*round(d[k]/L);
                mag += pow(d[k],2);
            }
                                    
            /* Get norm */
            mag = pow(mag,0.5);
            
            /* Apply cutoff for Lennard-Jones interaction */
            if (mag < 2.5*sig){
            
                /* Compute Lennard-Jones interaction */
                t1 = 4*ep*(pow(sig/mag,12) - pow(sig/mag,6));
                u_t += t1; /* [=] kcal mol-1*/
                /* printf("Mag: %3.5f \n", mag); */
                /* printf("LJ energy: %3.5f \n", t1); */
            
            }
            
            /* Apply cutoff for dipole-diple interactions (L/2) */
            if (mag < 1.0*L) {
            
                /* Get dot product of dipoles, mu_i * r_ij,  mu_j * r_ij  */
                dot1 = 0.0;
                dot2 = 0.0;
                dot3 = 0.0;
                
                for (k=0; k<3; k++){
                
                    dot1 += mu[i*3 + k] * mu[j*3 + k];
                    dot2 += mu[i*3 + k] * d[k];
                    dot3 += mu[j*3 + k] * d[k];
                
                }
                
                /* printf("Mag: %3.5f \n", mag); */
                /* printf("mu_i * mu_j: %3.5f \n", dot1); */
                /* printf("mu_i * r_ij: %3.5f \n", dot2); */
                /* printf("mu_j * r_ij: %3.5f \n", dot3); */
                
                /* Get dipole-dipole interaction */
                t2 = pow(mag,-3.0)*dot1 - 3*pow(mag, -5.0)*dot2*dot3;
                /* printf("Mag: %3.5f \n", mag); */
                /* printf("Dipole-dipole interaction: %3.5f \n", 332.07 *  t2); */
                u_t += 332.07 * t2; /* [=] kcal mol-1 */
        
            }  
        }
    }
    
    /* Assign overall energy */
    u[0] = u_t;
}

/* Calculates bonded interactions between polymer beads in electrolyte */

void pair_fene(double* u, int sz_u, double* r, int sz_r, double ep, double k_fene, double R0, double L, int N, int M) {

    int i, j;
    double mag, u_t, t1;
    double d[3];
    
    /* Initialize internal energy for this timestep*/
    u_t = 0.0;
    
    /* Loop over the particles (int N = number of particles) */
    for (i=0; i<N; i++){
    
        /* See if this is the last bead on the chain */
        if ( (i+1)%M != 0){
        
            /* Initialize magnitude and energy*/
            mag = 0.0;
            
            /* Get distance between the particles */
            for (j=0; j<3; j++) {
            
                d[j] = r[i*3 + j] - r[(i+1)*3 + j];
                d[j] = d[j]-L*round(d[j]/L);
                mag += pow(d[j],2);
            }
            
            /* Get norm */
            mag = pow(mag,0.5);
            
            /* Calculate bond energy */
            t1 = -0.5*k_fene*pow(R0,2)*log(1-pow(mag/R0,2)) + ep;
            /* printf("Mag: %3.5f \n", mag); */
            /* printf("FENE energy: %3.5f \n", t1); */
            u_t += t1; /* [=] kcal mol-1*/
            
            
        }
    }

    /* Assign overall energy */
    u[0] = u_t;
}

/* Calculates angle interactions (cosine) between polymer beads in electrolyte */

void pair_angle(double* u, int sz_u, double* r, int sz_r, double K, double L, int N, int M) {

    int i, j, ct;
    double maga, magb, magc, cos0, u_t, t1;
    double a[3], b[3], c[3];
    
    /* Initialize internal energy for this timestep*/
    u_t = 0.0;
    ct = 0;
    
    /* Loop over the particles (int N = number of particles) */
    for (i=0; i<N-2; i++){
    
        /* See if this is the second to last OR last bead on the chain */
        if ( (i+2)%M != 0 &&  (i+1)%M != 0){
        
            /* Initialize magnitude and energy*/
            maga = 0.0;
            magb = 0.0;
            magc = 0.0;
            /* printf("test 1: %d \n", (i+2)%M); */
            /* printf("test 1: %d \n", (i+1)%M); */
            
            /* Get distance between 1 and 2 (a) and between 2 and 3 (b) and between 1 and 3 (c) */
            for (j=0; j<3; j++) {
            
                a[j] = r[i*3 + j] - r[(i+1)*3 + j];
                b[j] = r[(i+1)*3 + j] - r[(i+2)*3 + j];
                c[j] = r[i*3 + j] - r[(i+2)*3 + j];
                a[j] = a[j]-L*round(a[j]/L);
                b[j] = b[j]-L*round(b[j]/L);
                c[j] = c[j]-L*round(c[j]/L);
                maga += pow(a[j],2);
                magb += pow(b[j],2);
                magc += pow(c[j],2);
            }
            
            /* Get norm */
            maga = pow(maga,0.5);
            magb = pow(magb,0.5);
            magc = pow(magc,0.5);
            
            /* Get cosine(theta) from law of cosines */
            
            cos0 = (pow(magc,2) - pow(maga,2) - pow(magb,2)) / (-2 * maga * magb);
            
            /* Calculate angle energy */
            ct += 1;
            /* printf("angle #: %d \n", ct); */
            t1 = K * (1 + cos0);
            /* printf("cos0: %3.5f \n", cos0); */
            /* printf("angle energy: %3.5f \n", t1); */
            u_t += t1; /* [=] kcal mol-1*/
            
            
            
        }
    }

    /* Assign overall energy */
    u[0] = u_t;
}

/* Calculates pairwise interactions between polymer beads and a single cation */

void pair_lj_dipole_cat(double* u, int sz_u, double* r, int sz_r, double* mu, int sz_mu, double* rc, int sz_rc, double ep, double sig, double L, int N, double q) {

    int i, j;
    double mag, u_t, dot1, t1, t2;
    double d[3];
    
    /* Initialize energy */
    u_t = 0.0;
        
    /* Loop over the polymer particles */
    for (i=0; i<N; i++){
        
        /* Initialize magnitude*/
        mag = 0.0;
    
        /* Get distance between the particles */
        for (j=0; j<3; j++) {
        
            d[j] = rc[j] - r[i*3 + j];
            d[j] = d[j]-L*round(d[j]/L);
            mag += pow(d[j],2);
        }
                                
        /* Get norm */
        mag = pow(mag,0.5);
        
        /* Apply cutoff for Lennard-Jones interaction */
        if (mag < 2.5*sig){
        
            /* Compute Lennard-Jones interaction */
            t1 = 4*ep*(pow(sig/mag,12) - pow(sig/mag,6));
            /* u_t += t1; [=] kcal mol-1*/
            /* printf("Mag: %3.5f \n", mag); */
            /* printf("LJ energy: %3.5f \n", t1); */
        
        }
        
        /* Apply cutoff for charge-diple interactions (L/2) */
        if (mag < 4.0) {
        
            /* Get dot product of dipoles, mu_i * r_ij,  mu_j * r_ij  */
            dot1 = 0.0;
            
            for (j=0; j<3; j++){
            
                dot1 += mu[i*3 + j] * d[j];
            
            }
            /* printf("Index: %d \n", i); */
            /* printf("Mag: %3.5f \n", mag); */
            /* printf("mu_i * r_ij: %3.5f \n", dot1); */
            /* printf("mu_i * r_ij: %3.5f \n", dot2); */
            /* printf("mu_j * r_ij: %3.5f \n", dot3); */
            
            /* Get dipole-charge interaction */
            t2 = q*pow(mag,-3.0)*dot1;
            /* printf("Mag: %3.5f \n", mag); */
            /* printf("Charge-dipole interaction: %3.5f \n", 332.07 *  t2); */
            u_t += 332.07 * t2; /* [=] kcal mol-1 */
            /* printf("Energy: %3.5f \n", u_t); */
    
        }  
    }
    
    /* Assign overall energy */
    u[0] = u_t;
    /* printf("Total charge-dipole interaction: %3.5f \n", u_t); */
}

/* Calculates pairwise interactions between polymer beads and a single cation */

void find_shell(int* index, int sz_ind, double* r, int sz_r, double* rc, int sz_rc, double L, int N) {

    int i, j;
    double mag;
    double d[3];
            
    /* Loop over the polymer particles to see who is in coordination shell */
    for (i=0; i<N; i++){
        
        /* Initialize magnitude*/
        mag = 0.0;
    
        /* Get distance between the particles */
        for (j=0; j<3; j++) {
        
            d[j] = rc[j] - r[i*3 + j];
            d[j] = d[j]-L*round(d[j]/L);
            mag += pow(d[j],2);
        }
                                
        /* Get norm */
        mag = pow(mag,0.5);
        
        if (mag < 4.0) {
        
            index[i] = 1;
        
        }
    }
}

/* Calculates bonded interactions between polymer beads in electrolyte */

void pair_fene_shell(double* u, int sz_u, double* r, int sz_r, int* index, int sz_ind, double ep, double k_fene, double R0, double L, int N, int M, int* bnd, int sz_bnd) {

    int i, j;
    double mag, u_t, t1;
    double d[3];
    
    /* Initialize internal energy for this timestep*/
    u_t = 0.0;
    
    /* Loop over the particles (int N = number of particles) */
    for (i=0; i<N; i++){
    
        /* Make sure this is not the last bead on the chain AND that it is in the coordination shell (index==1) */
        if ( (i+1)%M != 0 && index[i] > 0){
        
            bnd[0] += 1;
        
            /* printf("Index: %d \n", i); */
        
            /* Initialize magnitude and energy*/
            mag = 0.0;
            
            /* Get distance between the particles */
            for (j=0; j<3; j++) {
            
                d[j] = r[i*3 + j] - r[(i+1)*3 + j];
                d[j] = d[j]-L*round(d[j]/L);
                mag += pow(d[j],2);
            }
            
            /* Get norm */
            mag = pow(mag,0.5);
        
            
            /* Calculate bond energy */
            t1 = -0.5*k_fene*pow(R0,2)*log(1-pow(mag/R0,2)) + ep;
            /* printf("Mag: %3.5f \n", mag); */
            /* printf("FENE energy: %3.5f \n", t1); */
            u_t += t1; /* [=] kcal mol-1*/
            
            
        }
    }

    /* Assign overall energy */
    u[0] = u_t;
}

/* Calculates angle interactions (cosine) between polymer beads in electrolyte */
void pair_angle_shell(double* u, int sz_u, double* r, int sz_r, int* index, int sz_ind, double K, double L, int N, int M, int* ang, int sz_ang) {

    int i, j, ct;
    double maga, magb, magc, cos0, u_t, t1;
    double a[3], b[3], c[3];
    
    /* Initialize internal energy for this timestep*/
    u_t = 0.0;
    ct = 0;
    
    /* Loop over the particles (int N = number of particles) */
    for (i=0; i<N-2; i++){
    
        /* Make sure this is not the last bead on the chain, the second to last bead on the chain, and that the index is 1  */
        if ( (i+2)%M != 0 &&  (i+1)%M != 0 && index[i] > 0){
        
            ang[0] += 1;
        
            /* printf("Index: %d \n", i); */
        
            /* Initialize magnitude and energy*/
            maga = 0.0;
            magb = 0.0;
            magc = 0.0;
            /* printf("test 1: %d \n", (i+2)%M); */
            /* printf("test 1: %d \n", (i+1)%M); */
            
            /* Get distance between 1 and 2 (a) and between 2 and 3 (b) and between 1 and 3 (c) */
            for (j=0; j<3; j++) {
            
                a[j] = r[i*3 + j] - r[(i+1)*3 + j];
                b[j] = r[(i+1)*3 + j] - r[(i+2)*3 + j];
                c[j] = r[i*3 + j] - r[(i+2)*3 + j];
                a[j] = a[j]-L*round(a[j]/L);
                b[j] = b[j]-L*round(b[j]/L);
                c[j] = c[j]-L*round(c[j]/L);
                maga += pow(a[j],2);
                magb += pow(b[j],2);
                magc += pow(c[j],2);
            }
            
            /* Get norm */
            maga = pow(maga,0.5);
            magb = pow(magb,0.5);
            magc = pow(magc,0.5);
            
            /* Get cosine(theta) from law of cosines */
            
            cos0 = (pow(magc,2) - pow(maga,2) - pow(magb,2)) / (-2 * maga * magb);
            
            /* Calculate angle energy */
            ct += 1;
            /* printf("angle #: %d \n", ct); */
            t1 = K * (1 + cos0);
            /* printf("cos0: %3.5f \n", cos0); */
            /* printf("angle energy: %3.5f \n", t1); */
            u_t += t1; /* [=] kcal mol-1*/
            
            
            
        }
    }

    /* Assign overall energy */
    u[0] = u_t;
}

void pair_lj_dipole_shell(double* u, int sz_u, double* r, int sz_r, int* index, int sz_ind, double* mu, int sz_mu, double ep, double sig, double L, int N) {

    int i, j, k;
    double mag, u_t, dot1, dot2, dot3, t1, t2;
    double d[3];
    
    /* Initialize internal energy for this timestep*/
    u_t = 0.0;
    
    /* Loop over the particles (int N = number of particles) */
    for (i=0; i<N; i++){
        
        /* Loop over the upper half of the matrix */
        for (j=i+1; j<N; j++){
        
            /* See if BOTH i AND j are in the cluster around the ion*/
            
            if (index[i] > 0 && index[j]){
            
                /* printf("Index 1: %d \n", i); */
                /* printf("Index 2: %d \n", j); */
            
                /* Initialize magnitude and energy*/
                mag = 0.0;
            
                /* Get distance between the particles */
                for (k=0; k<3; k++) {
                
                    d[k] = r[i*3 + k] - r[j*3 + k];
                    d[k] = d[k]-L*round(d[k]/L);
                    mag += pow(d[k],2);
                }
                                        
                /* Get norm */
                mag = pow(mag,0.5);
                
                /* Apply cutoff for Lennard-Jones interaction */
                if (mag < 2.5*sig){
                
                    /* Compute Lennard-Jones interaction */
                    t1 = 4*ep*(pow(sig/mag,12) - pow(sig/mag,6));
                    u_t += t1; /* [=] kcal mol-1*/
                    /* printf("Mag: %3.5f \n", mag); */
                    /* printf("LJ energy: %3.5f \n", t1); */
                
                }
                
                /* Apply cutoff for dipole-diple interactions (L/2) */
                if (mag < 1.0*L) {
                
                    /* Get dot product of dipoles, mu_i * r_ij,  mu_j * r_ij  */
                    dot1 = 0.0;
                    dot2 = 0.0;
                    dot3 = 0.0;
                    
                    for (k=0; k<3; k++){
                    
                        dot1 += mu[i*3 + k] * mu[j*3 + k];
                        dot2 += mu[i*3 + k] * d[k];
                        dot3 += mu[j*3 + k] * d[k];
                    
                    }
                    
                    /* printf("Mag: %3.5f \n", mag); */
                    /* printf("mu_i * mu_j: %3.5f \n", dot1); */
                    /* printf("mu_i * r_ij: %3.5f \n", dot2); */
                    /* printf("mu_j * r_ij: %3.5f \n", dot3); */
                    
                    /* Get dipole-dipole interaction */
                    t2 = pow(mag,-3.0)*dot1 - 3*pow(mag, -5.0)*dot2*dot3;
                    /* printf("Mag: %3.5f \n", mag); */
                    /* printf("Dipole-dipole interaction: %3.5f \n", 332.07 *  t2); */
                    u_t += 332.07 * t2; /* [=] kcal mol-1 */
            
                } 
            
            } 
        }
    }
    
    /* Assign overall energy */
    u[0] = u_t;
}

