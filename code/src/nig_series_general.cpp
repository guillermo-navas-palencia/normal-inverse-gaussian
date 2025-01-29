#include <cmath>
#include <unordered_map>

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
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);  
  const double xmu = x - mu;
  const double omega = std::hypot(xmu, delta);

  // Constants
  const double xmu2 = xmu * xmu;
  const double aw = alpha * omega;
  const double oaw = 2.0 / aw;
  const double z = -omega * beta / alpha / xmu;
  const double v = xmu2 * alpha / omega;

  const double C = delta * xmu * alpha / omega * constants::oneopi * std::exp(
    delta * gamma + xmu * beta - aw);

  // Bessel recursion: First three terms
  const double k0 = specfun::bessel_k0_scaled(aw);
  const double k1 = specfun::bessel_k1_scaled(aw);

  std::unordered_map<int, double> bessel_map;
  bessel_map.insert({0, k0});
  bessel_map.insert({1, k1});

  // Start recursion
  double t = C;
  double s = 0.0;
  double sp = s;

  for (size_t k = 0; k < maxiter; k++)
  {
    // Compute polynomial A(k). First iteration
    if (bessel_map.find(k+1) == bessel_map.end()) {
      double cached = bessel_map[k - 1] + k * oaw * bessel_map[k];
      bessel_map.insert({k + 1, cached});
    }

    double sA = bessel_map[k + 1];
    double u = 1.0;
    double r = 1.0;

    for (size_t j = 1; j <= 2 * k + 1; j++)
    {
      int m = k + 1 - j;
      int am = std::fabs(m);

      u *= z;
      r *= (2.0 * k + 2 - j) / j;
      sA += r * u * bessel_map[am];
    }

    // New term
    t /= (2 * k + 1);
    s += t * sA;

    // Check convergence
    if (std::fabs(1.0 - s / sp) < eps) {
      return 0.5 + s;
    } else {
      t *= v;
      sp = s;
    }
  }

  return -1.0;
}


double incgamma_series_xmu(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const size_t maxiter = 100,
  const double eps = 5e-15
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);  
  const double xmu = x - mu;
  const double omega = std::hypot(xmu, delta);

  // Constants
  const double gw = gamma * omega;
  const double gow = gamma / omega;
  const double ogw = 1.0 / gw;
  const double z = std::sqrt(2.0 * gow) * xmu;

  // Bessel k0, k1, k_1/2, k_3/2
  const double k0 = specfun::bessel_k0_scaled(gw) * std::exp(-gw);
  const double k1 = specfun::bessel_k1_scaled(gw) * std::exp(-gw);

  const double k1o2 = std::sqrt(constants::pihalf / gw) * std::exp(-gw);
  const double k3o2 = k1o2 * (gw + 1.0) / gw;

  // Series 1: Bessel series
  double ki0 = k0;
  double ki1 = k1;
  double kh0 = k1o2;
  double kh1 = k3o2;

  double u = 1.0;
  double g_even = 1.0;
  double g_odd = constants::sqrtpi;

  // Iteration k=0
  double s1 = u * kh0 * g_odd;

  // Iteration k=1
  u *= z;
  s1 += u * ki1 * g_even;

  // Iteration k=2
  u *= z * 0.5;
  g_odd *= 0.5;
  s1 += u * kh1 * g_odd;

  // Start recursion
  double s1p = s1;

  for (size_t k = 3; k < maxiter; k++)
  {
    double g;
    double kn;

    if (k % 2 != 0)
    {
      // Gamma recursion
      g_even *= (k - 1) * 0.5;
      g = g_even;

      // Bessel recursion: k/2 + 1/2 is integer order
      kn = ki0 + (k - 1) * ogw * ki1;
      ki0 = ki1;
      ki1 = kn;

    } else 
    {
      // Gamma recursion
      g_odd *= (k - 1) * 0.5;
      g = g_odd;

      // Bessel recursion: k/2 + 1/2 is half order
      kn = kh0 + (k - 1) * ogw * kh1;
      kh0 = kh1;
      kh1 = kn;
    }

    // New term
    u *= z / k;
    s1 += u * g * kn;

    // Check convergence
    if (std::fabs(1.0 - s1 / s1p) < eps) {
      return 2.0 * std::sqrt(gow) * s1;
      // break;
    } else {
      s1p = s1;
    }
  }
 
  return s1;
 
  // Series 2: Incgamma and Bessel series
}


double hermite_series_xmu(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const size_t maxiter = 100,
  const double eps = 5e-15
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constants
  const double da = delta * alpha;
  const double xmu = x - mu;
  const double xmub = xmu * beta;
  const double oda = 2.0 / da;
  const double z = -0.5 * alpha / (beta * beta * delta);
  const double C = alpha / xmu * constants::oneopi * std::exp(delta * gamma - da);

  // Bessel recursion: First three terms
  const double k0 = specfun::bessel_k0_scaled(da);
  const double k1 = specfun::bessel_k1_scaled(da);

  std::unordered_map<int, double> bessel_map;
  bessel_map.insert({0, k0});
  bessel_map.insert({1, k1});
  bessel_map.insert({2, k0 + oda * k1});

  // Start recursion
  double v = 1.0;
  double s = 0.0;
  double sp = s;

  for (size_t k = 0; k < maxiter; k++)
  {
    // Computer polynomial A(k). First iteration   
    double u = 1.0;
    double r = 1.0 / std::tgamma(k + 1);
    double sA = r * k1;

    for (size_t j = 1; j <= floor(k / 2); j++)
    {
      if (bessel_map.find(j+1) == bessel_map.end())
      {
        // Apply Bessel recursion
        double cached = std::fma(j * oda, bessel_map[j], bessel_map[j-1]);
        bessel_map.insert({j+1, cached});
      }
      
      u *= z;
      r *= (2.0 * j - k - 1) * (2.0 * j - k - 2) / j;
      sA += r * u * bessel_map[j+1];
    }

    // New term
    s += sA * v / (k + 1);

    // Check convergence
    if (std::fabs(1.0 - s / sp) < eps) {
      return std::fma(C, s, nig_x_eq_mu(alpha, beta, delta));
    } else {
      v *= xmub;
      sp = s;
    }
      
  }

  return -1.0;
}


double hermite_series_beta(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const size_t maxiter = 100,
  const double eps = 5e-15
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);  
  const double xmu = x - mu;
  const double omega = std::hypot(xmu, delta);

  // Constants
  const double xmu2 = xmu * xmu;
  const double xmub = xmu * beta;
  const double gw = gamma * omega;
  const double ogw = 2.0 / gw;
  const double z = -0.5 * omega / (gamma * xmu2);
  const double C = -beta * delta * constants::oneopi * std::exp(gamma * delta - gw);

  // Bessel recursion: First three terms
  const double k0 = specfun::bessel_k0_scaled(gw);
  const double k1 = specfun::bessel_k1_scaled(gw);

  std::unordered_map<int, double> bessel_map;
  bessel_map.insert({-1, k1});
  bessel_map.insert({0, k0});
  bessel_map.insert({1, k1});

  // Start recursion
  double v = 1.0;
  double s = 0.0;
  double sp = s;

  for (size_t k = 0; k < maxiter; k++)
  {
    // Compute polynomial B(k). First iteration
    double u = 1.0;
    double r = 1.0 / std::tgamma(k + 1);
    double sB = r * k0;

    for (size_t j = 1; j <= floor(k / 2); j++)
    {
      if (bessel_map.find(j) == bessel_map.end())
      {
        // Apply Bessel recursion
        double cached = std::fma((j - 1) * ogw, bessel_map[j-1], bessel_map[j-2]);
        bessel_map.insert({j, cached});
      }
      
      u *= z;
      r *= (2.0 * j - k - 1) * (2.0 * j - k - 2) / j;
      sB += r * u * bessel_map[j];
    }

    // New term
    s += sB * v / (k + 1);

    // Check convergence
    if (std::fabs(1.0 - s / sp) < eps) {
      return std::fma(C, s, nig_beta_eq_zero(x, gamma, mu, delta));
    } else {
      v *= xmub;
      sp = s;
    }
  }

  return -1.0;
}


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
  // return asymptotic_xmu(x, alpha, beta, mu, delta);
  // return hermite_series_xmu(x, alpha, beta, mu, delta);
  // return hermite_series_beta(x, alpha, beta, mu, delta);
  // return bessel_series_xmu(x, alpha, beta, mu, delta);
  return incgamma_series_xmu(x, alpha, beta, mu, delta);
}