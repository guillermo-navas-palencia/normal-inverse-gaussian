#ifndef nig_hpp
#define nig_hpp


double nig_x_eq_mu(double alpha, double beta, double delta);
double nig_series_beta_zero(double x, double alpha, double mu, double delta);
double nig_cdf(double x, double alpha, double beta, double mu, double delta);

double bessel_k0_scaled(const double x);
double bessel_k1_scaled(const double x);

#endif // nig_hpp