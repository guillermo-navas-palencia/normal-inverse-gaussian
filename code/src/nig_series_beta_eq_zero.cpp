/*  Normal Inverse Gaussian cumulative distribution function for beta = 0.
 *  
 *  Implementation combines the following methods:
 *    - bessel_series: series expansion for |x-mu| -> 0
 *    - asymptotic_alpha: uniform asymptotic expansion for alphas -> inf
 *    - asymptotic_mu: asymptotic expansion for |x-mu| -> inf
 *    - nig_integration: numerical integration using tanh-sinh quadrature
 * 
 *  Guillermo Navas-Palencia <g.navas.palencia@gmail.com>
 *  Copyright (C) 2024
 */

#include <cmath>

#include <constants.hpp>
#include <nig.hpp>
#include <specfun.hpp>


double bessel_series(
  const double x,
  const double alpha,
  const double mu,
  const double delta,
  const size_t maxiter = 10000,
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
  const double oaw = 1.0 / aw;

  // Check if scaled version is required
  double C, k0, k1, t;
  const bool scaled = da > 705.342;

  // Constant C = delta * alpha / omega / pi * e^(delta * alpha)
  const double caux = delta * aow * constants::oneopi;

  if (scaled)
    C = std::fma(alpha, delta - omega, std::log(caux));
  else 
    C = caux * std::exp(da);

  // Series: compute the ratio of Bessel K_1 / K_2 for recursion
  if (scaled) {
    k0 = specfun::bessel_k0_scaled(aw);
    k1 = specfun::bessel_k1_scaled(aw);
    t = k1;
  } else {
    const double expaw = std::exp(-aw);
    k0 = specfun::bessel_k0_scaled(aw) * expaw;
    k1 = specfun::bessel_k1_scaled(aw) * expaw;
    t = C * k1;
  }

  double rp = k1 / k0;
  double s = t;

  // Start recursion
  double sp = s;
  for (size_t k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion: r = 1 / rp + 2k / aw
    double r = std::fma(2 * k, oaw, 1.0 / rp);

    // New term
    t *= z / (2 * k + 1) * r;
    s += t;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      if (scaled)
        return std::fma(xmu, std::exp(C + std::log(s)), 0.5);
      else
        return std::fma(xmu, s, 0.5);
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
  const size_t maxiter = 200,
  const double eps = 5e-13
)
{
  // Parameters
  const double xmu = x - mu;

  // Constants
  const double a = xmu / constants::sqrt2;
  const double b = delta / alpha;

  const double z = alpha * alpha / 2.0;
  const double oz = 1.0 / z;
  const double xi = b * z;
  const double C = delta * std::sqrt(z) * constants::oneosqrttwopi;

  const double a2 = a * a;
  const double t0 = 2.0 * xi;
  const double t1 = 2.0 * a2;
  const double t2 = 3.0 * b;
  const double t3 = 2.0 * b * b;
  const double t4 = xi * xi;
  const double t5 = 0.5 - t0;
  const double ot3 = 1.0 / t3;
  const double ob = 1.0 / b;
  const double sqrtob = std::sqrt(ob);
  const double t1mt2 = t1 - t2;
  const double b4 = 4.0 * b;

  // Compute the first two terms of Phi((x-mu) / sqrt(t)) at t=r
  double c0 = specfun::erfc(-a * sqrtob) * 0.5;
  double c1 = -a * 0.5 * std::exp(-a2 * ob) * constants::oneosqrtpi * sqrtob * ob;

  // Compute first elements of binomial sum of Bessel functions recursion
  const double ot0 = 1.0 / t0;
  const double kh = std::sqrt(constants::pihalf * ot0);
  const double khp2 = kh * (t0 + 1) * ot0;

  const double sqrt_xi = std::sqrt(xi);
  double q0 = 2.0 / sqrt_xi * kh;
  double q1 = 0.0;
  double q2 = 2.0 * sqrt_xi * xi * (khp2 - kh);

  double s = c0 * q0;  
  double num = oz;

  // Compute term k = 2 -> n = 0 => np1 = n + 1 = 1
  //      c1 * (t1 - t2)
  // ck = --------------
  //        t3 * 2
  double ck = (c1 * t1mt2) * ot3 * 0.5;
  num *= oz;
  s = std::fma(num * ck, q2, s);

  // Next iteration
  c0 = c1;
  c1 = ck;

  // Start recursion
  double sp = s;
  for (size_t k = 3; k < maxiter; k++)
  {
    // Compute recursion for Phi() at t=r
    size_t n = k - 2;
    size_t np1 = n + 1;

    //      (n + 1) * c1 * (t1 - 4bn - t2) - (2n^2 + n) * c0)
    // ck = -------------------------------------------------
    //                  t3 * (n+1) * (n+2)
    double ck0 = std::fma(n, -b4, t1mt2);
    double ck1 = std::fma(2 * n, n, n);
    double ck2 = std::fma(np1, c1 * ck0, -ck1 * c0);
    double ck = ck2 * ot3 / (np1 * (n + 2));

    // qk = (n + t5) * q2 + xi * (2n + 1/2) * q1 + n * t4 * q0
    double qk0 = std::fma(2, n, 0.5);
    double qk1 = std::fma(n + t5, q2, xi * qk0 * q1);
    double qk = std::fma(n * t4, q0, qk1);

    // New term: 1 / z^k * ck * qk
    num *= oz;
    s = std::fma(num * ck, qk, s);

    if (std::fabs(1.0 - sp / s) < eps) {
      return C * s;
    } else {
      c0 = c1;
      c1 = ck;
      q0 = q1;
      q1 = q2;
      q2 = qk;
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
  const size_t maxiter = 200,
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
  const double caux = delta / constants::pi;
  const double oaw = 1.0 / aw;

  // Check if scaled version is required
  double C, k0, k1, t;
  const bool scaled = (da > 705.342) | (aw > 705.342);

  if (scaled) {
    C = std::fma(alpha, delta - omega, std::log(caux));
    k0 = specfun::bessel_k0_scaled(aw);
    k1 = specfun::bessel_k1_scaled(aw);
    t = k0;
  } else {
    const double expaw = std::exp(-aw);
    k0 = specfun::bessel_k0_scaled(aw) * expaw;
    k1 = specfun::bessel_k1_scaled(aw) * expaw;
    C = caux * std::exp(da);
    t = C * k0;
  }
  
  double rp = k0 / k1;
  double s = t;

  // Start recursion
  double sp = s;
  for (size_t k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion: r = 1 / rp + 2(k-1) / aw
    double r = std::fma(2 * (k - 1), oaw, 1.0 / rp);

    // New term
    t *= -z * (2 * k - 1) * r;
    s += t;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      if (scaled)
        return -std::exp(C + std::log(s)) / xmu;
      else
        return -s / xmu;
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
  const double absxmu = std::fabs(xmu);
  const double xmu2 = xmu * xmu;
  const double omega = std::hypot(xmu, delta);
  const double da = delta * alpha;
  const double aow = alpha / omega;

  const bool use_series_c1 = (absxmu <= 5) & (aow <= 0.25) & (delta >= 2.0 * absxmu);
  const bool use_series_c2 = (xmu2 <= 1.25) & (aow <= 1.0);
  const bool use_series_c3 = (delta >= 1.0);

  const bool use_asymp_alpha_c1 = (xmu2 <= 2.5) & (da >= 200.0);
  const bool use_asymp_alpha_c2 = (alpha >= 5.0) & (delta >= 10.0);
  const bool use_asymp_xmu = (xmu2 >= 70) & (aow >= 1.0);

  if ((use_series_c1 | use_series_c2) & use_series_c3) {
    return bessel_series(x, alpha, mu, delta);
  } else if (use_asymp_alpha_c1 & use_asymp_alpha_c2) {
    return asymptotic_alpha(x, alpha, mu, delta);
  } else if (use_asymp_xmu) {
    if (xmu < 0.0)
      return asymptotic_xmu(x, alpha, mu, delta);
    else
      return 1.0 - asymptotic_xmu(-x, alpha, -mu, delta);
  } else {
    return nig_integration(x, alpha, 0.0, mu, delta, 1e-13, 14);
  }
}