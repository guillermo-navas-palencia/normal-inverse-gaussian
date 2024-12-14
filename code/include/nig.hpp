#ifndef nig_hpp
#define nig_hpp

double bessel_k0_scaled(const double x);
double bessel_k1_scaled(const double x);

double nig_integration(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double eps,
  const int maxlevel
);

double nig_x_eq_mu(double alpha, double beta, double delta);
double nig_beta_eq_zero(double x, double alpha, double mu, double delta);
double nig_general(double x, double alpha, double beta, double mu, double delta);

double nig_cdf(double x, double alpha, double beta, double mu, double delta);

#endif // nig_hpp