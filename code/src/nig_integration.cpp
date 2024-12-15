/*  Normal Inverse Gaussian cumulative distribution computed via tahn-sinh
 *  quadrature.
 *   
 *  Includes functionalities to devise optimal truncation point and magnitude
 *  estimation to obtain a robust implementation.
 *  
 *  The tanh-sinh quadrature code is based on the implementation in
 *  https://github.com/sigma-py/tanh-sinh with multiple optimizations.
 *  
 *  Guillermo Navas-Palencia <g.navas.palencia@gmail.com>
 *  Copyright (C) 2024
 */

#include <algorithm>
#include <cmath>

#include <constants.hpp>
#include <nig.hpp>


double norm_cdf(const double x)
{
  return 0.5 * std::erfc(-x * constants::osqrt2);
}


double norm_pdf(const double x)
{
  return std::exp(-x*x * 0.5) * constants::oneosqrttwopi;
}


void integrand_and_deriv(
  const double t,
  const double beta,
  const double xmu,
  const double gamma2,
  const double delta2,
  double &f,
  double &fp
)
{
  double sqrtt = std::sqrt(t);
  double t2 = t * t;
  double t32 = t * sqrtt;
  double ot = 1.0 / t;
  double ot2 = 1.0 / t2;
  double ot32 = 1.0 / t32;
  double betat = beta * t;
  double z1aux = xmu - betat;
  double z2aux = xmu + betat;

  double z = z1aux / sqrtt;

  // exp integrand at t
  double A = z2aux * ot32;
  double B = norm_pdf(z) / norm_cdf(z);

  f = 0.5 * (delta2 * ot2 - gamma2 - 3.0 * ot - A * B);

  // derivative exp integrand at t
  double d1 = -0.5 * (beta * ot32 - 1.5 * z2aux * ot * ot32) * B;
  const double AB = A * B;
  double d2 = -z * A * AB;
  double d3 = -AB * AB;

  fp = -delta2 * ot * ot2 + 1.5 * ot2 + d1 + d2 + d3;
}


double saddle_point(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double gamma,
  const double tol = 1e-4,
  const int maxiter = 10
)
{
  // Estimate magnitude to accurately calculate N to achieve eps rel error

  // Constants
  const double alpha2 = alpha * alpha;
  const double delta2 = delta * delta;
  const double gamma2 = gamma * gamma;
  const double xmu = x - mu;
  const double xmu2 = xmu * xmu;

  double a, c;
  if (x - mu < 0.0)
  {
    a = alpha2;
    c = delta2 + xmu2;
  } else {
    a = gamma2;
    c = delta2;
  }

  double x0 = (-1.5 + std::sqrt(2.25 + a * c)) / a;

  // Newton's iteration
  double f, fp;
  for (int k = 0; k < maxiter; k++)
  {
    integrand_and_deriv(x0, beta, xmu, gamma2, delta2, f, fp);
    x0 -= f / fp;

    if (std::abs(f) < tol)
      break;
  }

  return x0;
}


double estimate_magnitude(
  const double x0,
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double gamma,  
  const double C
)
{
  // Constants
  const double delta2 = delta * delta;
  const double gamma2 = gamma * gamma;
  const double xmu = x - mu;

  // Maximum integrand contribution
  const double sqrt_x0 = std::sqrt(x0);
  double z0 = xmu / sqrt_x0 - beta * sqrt_x0;
  double g1 = -0.5 * (delta2 / x0 + gamma2 * x0);
  double g2 = -1.5 * std::log(x0);
  double g3 = std::log(norm_cdf(z0));
  double gx0 = g1 + g2 + g3 + delta * gamma;
  
  return C * std::exp(gx0);
}


int truncation(
  const double delta,
  const double gamma,
  const double eps = 5e-16
)
{
  // Solve: 2 e^{-gamma^2 / 2 N} / N^(3/2) / gamma^2 = eps / C
  // N = 3/gamma^2 W0(gamma^2/(3u)), u = (gamma^2 eps / 2 / C) ^(2/3)

  // Use upper bound of the W0(x): W0(x) < log(x)^(log(x) / (1 + log(x)))
  // To avoid overflow/underflow perform computation using logarithms

  // Constants
  const double gamma2 = gamma * gamma;
  const double loggamma2 = std::log(gamma2);

  // log C = log(delta) + delta * gamma - 1/2 * log(2 pi)
  const double logC = std::log(delta) + delta * gamma - 0.9189385332046727;

  // log y = log(gamma^2) - log(3) - 2/3 (log(gamma^2) + log(eps) - log(2) - logC)
  const double logy =  1./3 * loggamma2 - 2./3 * (std::log(eps) - logC) - 0.636514168294813;
  const double lambertwy = std::pow(logy, logy / (1.0 + logy));

  return (int) std::ceil(3.0 / gamma2 * lambertwy);
}


double estimate_h(
  const double tau,
  const double tol = 1e-10,
  const int maxiter = 10
)
{
  // Solve log(pi / tau) = pi / 2 - exp(x) - x - log(x)

  // Initial guess x0
  const double piotau = constants::pi / tau;
  double x = std::log(2.0 * constants::oneopi * std::log(piotau));

  // Newton's iteration
  double aux = constants::pihalf * std::exp(x);
  double fx = aux - x - std::log(x * piotau);

  for (int k = 0; k < maxiter; k++)
  {
    double fxp = aux - 1.0 - 1.0 / x;
    x -= fx / fxp;

    aux = constants::pihalf * std::exp(x);
    fx = aux - x - std::log(x * piotau);

    if (std::abs(fx) < tol)
      return x;
  }

  return x;
}


double lambertwm1(const double x)
{
  // D. A. Barry, L. Li, and D.-S. Jeng, “Comments on “Numerical evaluation
  // of the Lambert W function and application to generation of generalized
  // Gaussian noise with exponent 1/2”,” IEEE Trans. Signal Process.,
  // vol. 52, no. 5, pp. 1456–1458, May 2004.
  const double a = 0.3205;
  const double logmx = std::log(-x);

  return logmx - 2.0 / a * (1.0 - 1.0 / (1.0 + a * std::sqrt(-0.5 * (1.0 + logmx))));
}


double nig_integration(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double eps = 1e-15,
  const int maxlevel = 10
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constants
  const double xmu = x - mu;
  const double C = delta * constants::oneosqrttwopi;

  // Estimate integrand saddle point
  double x0 = saddle_point(x, alpha, beta, mu, delta, gamma);

  // Estimate integrand magnitude to achieve eps relative error
  double magnitude = estimate_magnitude(x0, x, alpha, beta, mu, delta, gamma, C);

  // Estimate the truncation point N
  int N = truncation(delta, gamma, eps * std::min(magnitude, 1.0));

  // Estimate initial step size h
  const double eps2 = eps * eps;
  double h = estimate_h(eps2);

  // Start integration
  const double alpha2 = N * 0.5;

  // Level 0
  double y = alpha2;
  double weight = h * alpha2 * constants::pihalf;

  // Eval f(y0) and f(N - y0)
  double xl = y;
  double oxl = 1.0 / xl;
  double sqrt_oxl = std::sqrt(oxl);
  double rl = delta - gamma * xl;
  double phil = norm_cdf((xmu - beta * xl) * sqrt_oxl);
  double fl = phil * std::exp(-0.5 * rl*rl * oxl) * sqrt_oxl * oxl;

  double estimate = fl * weight;
  h *= 0.5;  

  for (int level = 1; level <= maxlevel; level++)
  {
    // Estimate j using approximation W_{-1}(x)
    double japprox = std::log(-2.0 / constants::pi * lambertwm1(-eps2 / h * 0.5)) / h;
    unsigned int j = (int) std::ceil(japprox);

    double sum = 0.0;

    for (unsigned int i = 1; i < j + 1; i += 2)
    {
      double t = h * i;
      double sinh_t = constants::pihalf * std::sinh(t);
      double cosh_t = constants::pihalf * std::cosh(t);
      double cosh_sinh_t = 1.0 / std::cosh(sinh_t);
      double exp_sinh_t = std::exp(sinh_t);

      double y = alpha2 / exp_sinh_t * cosh_sinh_t;
      double weight = h * alpha2 * cosh_t * (cosh_sinh_t * cosh_sinh_t);

      // Eval f(y0) and f(N - y0)
      double xl = y;
      double rl = delta - gamma * xl;
      double oxl = 1.0 / xl;
      double sqrt_oxl = std::sqrt(oxl);      
      double phil = norm_cdf((xmu - beta * xl) * sqrt_oxl);
      double fl = phil * std::exp(-0.5 * rl*rl * oxl) * sqrt_oxl * oxl;

      double xr = N - xl;
      double oxr = 1.0 / xr;
      double sqrt_oxr = std::sqrt(oxr);
      double rr = delta - gamma * xr;
      double phir = norm_cdf((xmu - beta * xr) * sqrt_oxr);
      double fr = phir * std::exp(-0.5 * rr*rr * oxr) * sqrt_oxr * oxr;

      sum += (fl + fr) * weight;
    }

    double f = sum + estimate * 0.5;

    if (std::abs(1.0 - estimate / f) < eps)
    {
      return C * f;
    } else {
      estimate = f;
      h *= 0.5;
    }
  }

  double integral = C * estimate;

  // Check if order of magnitude of the computed integral is close to the
  // estimate. If the magnitudes differ significantly, an error during the
  // computation occurred. Then, return the estimated magnitude.
  if (std::abs(std::log10(integral) - std::log(magnitude)) > 5)
    return estimate;
  else
    return integral;
}