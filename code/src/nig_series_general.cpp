#include <cmath>

#include <nig.hpp>
#include <constants.hpp>
#include <specfun.hpp>

#include <iostream>
#include <iomanip>


double bessel_series_xmu(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const size_t maxiter = 100,
  const double eps = 5e-15
)
{}


double incgamma_series_xmu(
  const double x,
  const double alpha,
  const double mu,
  const double delta,
  const size_t maxiter = 10000,
  const double eps = 5e-15
)
{}


double hermite_series_beta(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const size_t maxiter = 100,
  const double eps = 5e-15
)
{}


double asymptotic_delta(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const size_t maxiter = 100,
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
  const double oda = 1.0 / da;
  const double z = -0.5 * alpha / (beta * beta * delta);

  // Ratio of scaled Bessel functions recursion
  const double C = alpha / beta * constants::oneopi * std::exp(delta * gamma - da);
  const double k0 = specfun::bessel_k0_scaled(da);
  const double k1 = specfun::bessel_k1_scaled(da);
  double rp = k1 / k0;

  // Incomplete gamma recursion
  double qp = expxmub;
  double v = 1.0;

  double t = C * k1;
  double s = t * expxmub;

  // Start recursion
  double sp = s;
  for (size_t k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion: r = 1.0 / rp + 2 * k / da
    double r = std::fma(2 * k, oda, 1.0 / rp);

    // Regularized incomplete gamma
    size_t n = 2 * k;
    size_t m = n * (n - 1);
    v *= xmub2;
    double q = m * qp + expxmub * (v * (1.0 - n * oxmub));

    // New term
    t *= z * r / k;
    s += t * q;

    // Check convergence
    if (std::fabs(1.0 - s / sp) < eps) {
      return (beta > 0.0) ? s : 1.0 + s;
    } else {
      rp = r;
      qp = q;
      sp = s;
    }
  }

  return -1.0;
}


double asymptotic_xmu(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const size_t maxiter = 100,
  const double eps = 5e-13
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);  
  const double xmu = x - mu;
  const double omega = std::hypot(xmu, delta);

  // Constants
  const double xmu2 = xmu * xmu;
  const double gw = gamma * omega;
  const double xmub = xmu * beta;
  const double xmub2 = xmub * xmub;
  const double oxmub = 1.0 / xmub;
  const double expxmub = std::exp(xmub);
  const double ogw = 1.0 / gw;
  const double z = -0.5 * omega / (gamma * xmu2);

  // Ratio of scaled Bessel functions recursion
  const double C = -delta / xmu * constants::oneopi * std::exp(gamma * delta - gw);
  const double k0 = specfun::bessel_k0_scaled(gw);
  const double k1 = specfun::bessel_k1_scaled(gw);
  double rp = k0 / k1;

  // Incomplete gamma recursion
  double qp = expxmub;
  double v = 1.0;

  double t = C * k0;
  double s = t * expxmub;

  // Start recursion
  double sp = s;
  for (size_t k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion: r = 1.0 / rp + 2 * k / da
    double r = std::fma(2 * (k - 1), ogw, 1.0 / rp);

    // Regularized incomplete gamma
    size_t n = 2 * k;
    size_t m = n * (n - 1);
    v *= xmub2;
    double q = m * qp + expxmub * (v * (1.0 - n * oxmub));

    // New term
    t *= z * r / k;
    s += t * q;

    // Check convergence
    if (std::fabs(1.0 - s / sp) < eps) {
      return (beta > 0.0) ? s : 1.0 + s;
    } else {
      rp = r;
      qp = q;
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
  // return asymptotic_delta(x, alpha, beta, mu, delta);
  return asymptotic_xmu(x, alpha, beta, mu, delta);
}