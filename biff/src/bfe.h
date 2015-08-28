extern void c_density(double *xyz, int K, double M, double r_s,
                      double *Snlm, double *Tnlm,
                      int nmax, int lmax, double *dens);

extern void c_potential(double *xyz, int K,
                        double G, double M, double r_s,
                        double *Snlm, double *Tnlm,
                        int nmax, int lmax, double *val);

extern void c_gradient(double *xyz, int K,
                       double G, double M, double r_s,
                       double *Snlm, double *Tnlm,
                       int nmax, int lmax, double *grad);

extern double scf_value(double t, double *pars, double *q);
extern void scf_gradient(double t, double *pars, double *q, double *grad);
