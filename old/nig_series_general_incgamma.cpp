#include <cmath>
#include <unordered_map>

// #include <nig.hpp>
// #include <constants.hpp>
// #include <specfun.hpp>


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