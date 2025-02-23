/*  Normal Inverse Gaussian cumulative distribution function.
 *  
 *  Implementation combines two special cases and the general case.
 *  Use numerical integration as a backup.
 * 
 *  Guillermo Navas-Palencia <g.navas.palencia@gmail.com>
 *  Copyright (C) 2024
 */

#include <nig.hpp>


double nig_cdf(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta
)
{
  double cdf;

  if (x == mu)
    cdf = nig_x_eq_mu(alpha, beta, delta);
  else if (beta == 0.0)
    cdf = nig_beta_eq_zero(x, alpha, mu, delta);
  else
    cdf = nig_general(x, alpha, beta, mu, delta);

  if (cdf == -1.0)
    return nig_integration(x, alpha, beta, mu, delta, 1e-13, 14);
  else
    return cdf;
}