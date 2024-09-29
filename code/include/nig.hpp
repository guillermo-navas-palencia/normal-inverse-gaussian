#ifndef nig_hpp
#define nig_hpp


double nig_x_eq_mu(double alpha, double beta, double delta);
double nig_series_beta_zero(double x, double alpha, double mu, double delta);
double nig_cdf(double x, double alpha, double beta, double mu, double delta);

double bessel_k0_scaled(const double x);
double bessel_k1_scaled(const double x);

int truncation(
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double eps);

double estimate_h(
  const double tau,
  const double tol,
  const int maxiter);

double nig_integration(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double eps,
  const int maxlevel
);

#endif // nig_hpp