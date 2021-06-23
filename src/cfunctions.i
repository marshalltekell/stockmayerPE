%module cfunctions

%{
    #define SWIG_FILE_WITH_INIT
    #include "cfunctions.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (int* ARGOUT_ARRAY1, int DIM1) {(int* rangevec, int n)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* invec, int n)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* r, int n_t1)}
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* ct, int n_t2)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* rt, int n_t3)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* x, int sz_x)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* y, int sz_y)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* z, int sz_z)}
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* q, int sz_q)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* g, int nbins1)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* r, int sz_r)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* r_RDF, int nbins2)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* f, int nbins3)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* qv, int sz_qv)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* fqt, int n_t4)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* muac, int n_t5)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* mu_int, int n_t6)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* ssf, int n_k)}
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* nq, int sz_nq)}
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* p1, int sz_p)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* r, int sz_r)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* mu, int sz_mu)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* u, int sz_u)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* rc, int sz_rc)}
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* index, int sz_ind)}
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* bnd, int sz_bnd)}
%apply (int* INPLACE_ARRAY1, int DIM1) {(int* ang, int sz_ang)}

%include "cfunctions.h"