extern double rho_nl(double s, int n, int l);
extern double rho_nlm(double s, double phi, double X, int n, int l, int m);

extern double phi_nl(double s, int n, int l);
extern double phi_nlm(double s, double phi, double X, int n, int l, int m);

extern void sph_grad_phi_nlm(double s, double phi, double X, int n, int l, int m, double *sphgrad);
