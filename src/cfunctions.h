/* mean squared displacement */
void calcMSD(double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t);
void logMSD(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t);
void calcM4D(double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t);
void logM4D(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N, int n_t);

/* radial distribution function */
void RDF1(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L);
void RDF2(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L, int Nc);
void RDF3(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L, int Nc);
void RDF4(double* g, int nbins1, double* r, int sz_r, int nbins, int N, double rho, double L, int N2);
void fcn(double* f, int nbins3, double* r_RDF, int nbins2, double* g, int nbins1, int nbins, double L, double rho);

/* ionic conductivity */
void lam(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);
void loglam(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);
void lam_cross1(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);
void loglam_cross1(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);
void lam_cross2(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);
void loglam_cross2(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);
void lam_cross3(double* rt, int n_t3, double* r, int n_t1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);
void loglam_cross3(double* r, int n_t1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int* q, int sz_q, int N, int n_t);

/* self-intermediate scattering */
void fqt(double* fqt, int n_t4, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int N, int n_t, int num_q);
void logfqt(double* fqt, int n_t4, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int N, int n_t, int num_q);

/* static structure factor */
void SSF1(double* ssf, int n_k, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int* nq, int sz_nq, int N, int n_t, int num_p);
void SSF2(double* ssf, int n_k, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, double* qv, int sz_qv, int* nq, int sz_nq, int N1, int N2, int n_t, int num_p);

/* van Hove function */
void logvanHove(double* g, int nbins1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N1, int N2, int nbins, int n_t, double L);
void vanHove(double* g, int nbins1, int* ct, int n_t2, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N1, int N2, int nbins, int n_t, double L);
void vanHove_short(double* g, int nbins1, double* x, int sz_x, double* y, int sz_y, double* z, int sz_z, int N1, int N2, int nbins, int n_t, int n_t2, double L);

/* dipole autocorrelation function */
void MUAC(double* muac, int n_t5, double* mu_int, int n_t6, int* ct, int n_t2, int N,  int n_t);
void MUAC_log(double* muac, int n_t5, double* mu_int, int n_t6, int N,  int n_t);

/* coordination number statistics */
void pofn(int* p1, int sz_p, double* r, int sz_r, int N1, int N2, double lo, double hi, double L);
void pofM(int* p1, int sz_p, double* r, int sz_r, int N1, int N2, double lo, double hi, double L, int Nc);

/* potential energy calculations */
void pair_lj_dipole(double* u, int sz_u,  double* r, int sz_r, double* mu, int sz_mu, double ep, double sig, double L, int N);
void pair_lj_dipole_cat(double* u, int sz_u, double* r, int sz_r, double* mu, int sz_mu, double* rc, int sz_rc, double ep, double sig, double L, int N, double q);
void pair_fene(double* u, int sz_u, double* r, int sz_r, double ep, double k_fene, double R0, double L, int N, int M);
void pair_angle(double* u, int sz_u, double* r, int sz_r, double K, double L, int N, int M);
void find_shell(int* index, int sz_ind, double* r, int sz_r, double* rc, int sz_rc, double L, int N);
void pair_fene_shell(double* u, int sz_u, double* r, int sz_r, int* index, int sz_ind, double ep, double k_fene, double R0, double L, int N, int M, int* bnd, int sz_bnd);
void pair_angle_shell(double* u, int sz_u, double* r, int sz_r, int* index, int sz_ind, double K, double L, int N, int M, int* ang, int sz_ang);
void pair_lj_dipole_shell(double* u, int sz_u, double* r, int sz_r, int* index, int sz_ind, double* mu, int sz_mu, double ep, double sig, double L, int N);