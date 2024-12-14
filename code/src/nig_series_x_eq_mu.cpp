#include <cmath>

#include <nig.hpp>
#include <constants.hpp>
#include <specfun.hpp>


double bessel_series(
  const double alpha,
  const double beta,
  const double delta,
  const int maxiter = 10000,
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

  // Check if scaled version is required
  double C, k0, k1, t;
  const bool scaled = dg > 705.342;
  const double caux = beta * delta / constants::pi;

  if (scaled)
    C = std::log(caux) + delta * (gamma - alpha);
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
  for (int k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion
    double r = 1.0 / rp + 2 * (k - 1) / ad;

    // New term
    t *= z / (2 * k + 1) * r;
    s += t;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      if (scaled)
        return 0.5 - std::exp(C + std::log(s));
      else
        return 0.5 + s;
    }
    else {
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

  double rba = beta / alpha;

  if ((alpha <= 10.0) & (delta <= 10.0) & (beta <= 1.5) & (rba <= 0.9))
    return bessel_series(alpha, beta, delta);
  else
    return nig_integration(0.0, alpha, beta, 0.0, delta, 1e-13, 14);
}