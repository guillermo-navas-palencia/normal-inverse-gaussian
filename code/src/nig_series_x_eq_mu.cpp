#include <cmath>

#include <nig.hpp>
#include <constants.hpp>
#include <specfun.hpp>

#include <iostream>
#include <iomanip>


double bessel_series(
  const double alpha,
  const double beta,
  const double delta,
  const size_t maxiter = 10000,
  const double eps = 5e-16
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constant
  const double ad = alpha * delta;
  const double dg = delta * gamma;
  const double doa = delta / alpha;
  const double z = beta * beta * doa;
  const double oad = 1.0 / ad;

  // Check if scaled version is required
  double C, k0, k1, t;
  const bool scaled = dg > 705.342;
  const double caux = beta * delta / constants::pi;

  if (scaled)
    C = std::fma(delta, gamma - alpha, std::log(caux));
  else
    C = -caux * std::exp(dg);

  // Series: compute the ratio of Bessel K_0 / K_1 for recursion
  if (scaled) {
    k0 = specfun::bessel_k0_scaled(ad);
    k1 = specfun::bessel_k1_scaled(ad);
    t = k0;
  } else {
    const double expad = std::exp(-ad);
    k0 = specfun::bessel_k0_scaled(ad) * expad;
    k1 = specfun::bessel_k1_scaled(ad) * expad;
    t = C * k0;
  }

  double rp = k0 / k1;
  double s = t;

  // Start recursion
  double sp = s;
  for (unsigned int k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion: r = 1.0 / rp + 2 * (k - 1) / ad
    double r = std::fma(2 * (k - 1), oad, 1.0 / rp);

    // New term
    t *= z / (2 * k + 1) * r;
    s += t;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps) {
      return (scaled) ? 0.5 - std::exp(C + std::log(s)) : 0.5 + s;
    } else {
      rp = r;
      sp = s;
    }
  }

  return -1.0;
}


double asymptotic_delta(
  const double alpha,
  const double beta,
  const double delta,
  const size_t maxiter = 100,
  const double eps = 5e-13  
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constants
  const double da = delta * alpha;
  const double oda = 1.0 / da;
  const double z = -2.0 * alpha / (beta * beta * delta);

  // Ratio of scaled Bessel functions recursion
  const double C = alpha / beta * constants::oneopi * std::exp(delta * gamma - da);
  const double k0 = specfun::bessel_k0_scaled(da);
  const double k1 = specfun::bessel_k1_scaled(da);
  double rp = k1 / k0;

  // Gamma recursion
  double t = C * k1;
  double s = t;

  // Start recursion
  double sp = s;
  for (unsigned k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion: r = 1.0 / rp + 2 * k / da
    double r = std::fma(2 * k, oda, 1.0 / rp);

    // New term
    t *= z * r * (k - 0.5);
    s += t;

    if (std::fabs(1.0 - sp / s) < eps) {
      return (beta > 0.0) ? s : 1.0 + s;
    } else {
      rp = r;
      sp = s;
    }
  }

  return -1.0;  
}



double nig_x_eq_mu(const double alpha, const double beta, const double delta)
{
  if (beta == 0.0)
    return 0.5;

  double rba = std::fabs(beta) / alpha;

  if ((alpha <= 10.0) & (delta <= 10.0) & (beta <= 1.5) & (rba <= 0.9))
    return bessel_series(alpha, beta, delta);
  else if ((rba >= 0.75) & (delta * alpha >= 300.0) & (delta >= 15.0))
    return asymptotic_delta(alpha, beta, delta);
  else
    return nig_integration(0.0, alpha, beta, 0.0, delta, 1e-13, 14);
}