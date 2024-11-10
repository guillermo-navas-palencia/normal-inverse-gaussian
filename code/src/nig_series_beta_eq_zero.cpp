#include <cmath>

#include <nig.hpp>
#include <constants.hpp>


double bessel_series(
  const double x,
  const double alpha,
  const double mu,
  const double delta,
  const int maxiter = 10000,
  const double eps = 5e-15
)
{
  // Parameters
  const double xmu = x - mu;
  const double xmu2 = xmu * xmu;
  const double omega = std::hypot(xmu, delta);

  // Constants
  const double aw = alpha * omega;
  const double aow = alpha / omega;
  const double da = delta * alpha;
  const double z = xmu2 * aow;

  // Check if scaled version is required
  double C, k0, k1, t;
  const bool scaled = da > 705.342;
  const double caux = delta * xmu * aow / constants::pi;

  if (scaled)
    C = std::log(caux) + alpha * (delta - omega);
  else 
    C = caux * std::exp(da);

  // Series: compute the ratio of Bessel K_1 / K_2 for recursion
  if (scaled) {
    k0 = bessel_k0_scaled(aw);
    k1 = bessel_k1_scaled(aw);
    t = k1;
  } else {
    k0 = std::cyl_bessel_k(0, aw);
    k1 = std::cyl_bessel_k(1, aw);
    t = C * k1;
  }

  double rp = k1 / k0;
  double s = t;

  // Start recursion
  double sp = s;
  for (int k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion
    double r = 1.0 / rp + 2 * k / aw;

    // New term
    t *= z / (2 * k + 1) * r;
    s += t;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      if (scaled)
        return 0.5 + std::exp(C + std::log(s));
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


double asymptotic_alpha(
  const double x,
  const double alpha,
  const double mu,
  const double delta,
  const int maxiter = 200,
  const double eps = 5e-13
)
{
  // Parameters
  const double xmu = x - mu;

  // Constants
  const double a = xmu / constants::sqrt2;
  const double b = delta / alpha;

  const double z = alpha * alpha / 2.0;
  const double ovz = 1.0 / z;
  const double xi = b * z;

  const double C = delta * std::sqrt(z) * constants::oneosqrttwopi;

  const double t0 = 2.0 * xi;
  const double t1 = 2.0 * a * a;
  const double t2 = 3.0 * b;
  const double t3 = 2.0 * b * b;
  const double t4 = xi * xi;
  const double t5 = 0.5 - t0;

  // Compute the first two terms of Phi((x-mu) / sqrt(t)) at t=b
  double c0 = std::erfc(-a / std::sqrt(b)) / 2.0;
  double c1 = -a / 2.0 * std::exp(-a * a / b) * constants::oneosqrtpi / std::pow(b, 1.5);

  // Compute first elements of binomial sum of Bessel functions recursion
  const double kh = std::sqrt(constants::pi / 2 / t0);
  const double khp2 = kh * (t0 + 1) / t0;

  double q0 = 2.0 / std::sqrt(xi) * kh;
  double q1 = 0.0;
  double q2 = 2.0 * std::pow(xi, 1.5) * (khp2 - kh);

  double s = c0 * q0;  
  double num = ovz;

  double sp = s;
  double qk;

  for (int k = 2; k < maxiter; k++)
  {
    int n = k - 2;
    int np1 = n + 1;
    double ck = (np1 * c1 * (t1 - 4.0 * b * n - t2) -
                 (2 * n * n + n) * c0) / (t3 * np1 * (n + 2));

    if (k >= 3) {
      qk = (n + t5) * q2 + xi * (2 * n + 0.5) * q1 + n * t4 * q0;
    } else {
      qk = q2;
    }

    // New term
    num *= ovz;
    s += num * ck * qk;

    if (std::fabs(1.0 - sp / s) < eps) {
      return C * s;
    } else {
      c0 = c1;
      c1 = ck;

      if (k >= 3) {
        q0 = q1;
        q1 = q2;
        q2 = qk;
      }

      sp = s;
    }
  }

  return -1.0;
}


double asymptotic_xmu(
  const double x,
  const double alpha,
  const double mu,
  const double delta,
  const int maxiter = 200,
  const double eps = 5e-13
)
{
  // Parameters
  const double xmu = x - mu;
  const double xmu2 = xmu * xmu;
  const double omega = std::hypot(xmu, delta);
  const double aw = alpha * omega;
  const double da = delta * alpha;
  const double woa = omega / alpha;
  const double z = woa / xmu2;

  // Check if scaled version is required
  double C, k0, k1;
  const bool scaled = da > 705.342;

  if (scaled) {
    C = -delta * std::exp(da - aw) / constants::pi / xmu;
    k0 = bessel_k0_scaled(aw);
    k1 = bessel_k1_scaled(aw);    
  } else {
    k0 = std::cyl_bessel_k(0, aw);
    k1 = std::cyl_bessel_k(1, aw);
    C = -delta * std::exp(da) / constants::pi / xmu;
  }
  
  double t = C * k0;

  double rp = k0 / k1;
  double s = t;

  // Start recursion
  double sp = s;
  for (int k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion
    double r = 1.0 / rp + 2 * (k - 1) / aw;

    // New term
    t *= -z * (2 * k - 1) * r;
    s += t;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      return s;
    } else {
      rp = r;
      sp = s;
    }
  }

  return -1.0;
}


double nig_beta_eq_zero(
  const double x,
  const double alpha,
  const double mu,
  const double delta
)
{
  const double xmu = x - mu;
  const double xmu2 = xmu * xmu;
  const double omega = std::hypot(xmu, delta);
  const double da = delta * alpha;
  const double aw = alpha * omega;
  const double aow = alpha / omega;

  const bool use_series_c1 = (std::fabs(xmu) <= 10) & (aow <= 0.25);
  const bool use_series_c2 = (xmu2 <= 1.25) & (aw <= 750.0);

  const bool use_asymp_a = (xmu2 <= 2.5) & (da >= 150.0) & (alpha >= 5.0);
  const bool use_asymp_xmu = (xmu2 >= 70) & (aow >= 1.0);

  if (use_series_c1 | use_series_c2) {
    // std::cout << "bessel_series" << std::endl;
    if (xmu > 0.0)
      return bessel_series(x, alpha, mu, delta);
    else
      return 1.0 - bessel_series(-x, alpha, -mu, delta);
  } else if (use_asymp_a){
    // std::cout << "asymptotic_alpha" << std::endl;
    return asymptotic_alpha(x, alpha, mu, delta);
  } else if (use_asymp_xmu) {
    // std::cout << "asymptotic_xmu" << std::endl;
    if (xmu < 0.0)
      return asymptotic_xmu(x, alpha, mu, delta);
    else
      return 1.0 - asymptotic_xmu(-x, alpha, -mu, delta);
  } else {
    // std::cout << "integration" << std::endl;
    return nig_integration(x, alpha, 0.0, mu, delta, 1e-13, 14);
  }
}