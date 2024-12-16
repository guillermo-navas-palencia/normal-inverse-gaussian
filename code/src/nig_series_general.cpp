#include <cmath>

#include <nig.hpp>
#include <constants.hpp>
#include <specfun.hpp>

#include <iostream>
#include <iomanip>


double asymptotic_delta(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const int maxiter = 100,
  const double eps = 5e-13
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constants
  const double da = delta * alpha;
  const double xmu = x - mu;
  const double xmub = xmu * beta;
  const double xmub2 = xmub * xmub;
  const double oxmub = 1.0 / xmub;
  const double expxmub = std::exp(xmub);
  const double z = -0.5 * alpha / (beta * beta * delta);

  // Ratio of scaled Bessel functions recursion
  const double C = alpha / beta / constants::pi * std::exp(delta * gamma - da);
  const double k0 = specfun::bessel_k0_scaled(da);
  const double k1 = specfun::bessel_k1_scaled(da);
  double rp = k1 / k0;

  // Incomplete gamma recursion
  double q = expxmub;
  double rg = xmub2;

  double t = C * k1;
  double s = t * q;

  // Start recursion
  double sp = s;
  for (int k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion
    double r = 1.0 / rp + 2 * k / da;

    // Regularized incomplete gamma
    q += expxmub * rg * (0.5 / k - oxmub);

    // New term
    t *= z * (4.0 * k - 2.0) * r;
    s += t * q;

    std::cout << k << " " << std::setprecision(16) << std::fabs(t) << " " << q  << " " << s << " " << std::fabs(s - sp) << std::endl;

    // Check convergence
    bool check_1 = (beta < 0.0) & (std::fabs(sp - s) * 100 < eps);
    bool check_2 = (beta > 0.0) & (std::fabs(1.0 - s / sp) * 10 < eps);

    if (check_1 | check_2) {
      return (beta > 0.0) ? s : 1.0 + s;
    } else {
      rp = r;
      rg *= xmub2 / ((2 * k) * (2 * k + 1));
      sp = s;
    }
  }

  return -1.0;
}


double nig_general(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta
)
{
  return asymptotic_delta(x, alpha, beta, mu, delta);
}