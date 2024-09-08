#include <nig.hpp>


double nig_cdf(
  double x,
  double alpha,
  double beta,
  double mu,
  double delta
)
{
  if (x == mu)
    return nig_series_x_eq_mu(alpha, beta, delta);
  else if (beta == 0.0)
    return nig_series_beta_zero(x, alpha, mu, delta);

  return -2.0;
}