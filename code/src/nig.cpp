#include <nig.hpp>


double nig_cdf(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta
)
{
  if (x == mu)
    return nig_x_eq_mu(alpha, beta, delta);
  else if (beta == 0.0)
    return nig_beta_eq_zero(x, alpha, mu, delta);
  else
    return nig_general(x, alpha, beta, mu, delta);
}