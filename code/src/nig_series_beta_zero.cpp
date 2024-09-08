#include <cmath>

#include <nig.hpp>
#include <constants.hpp>

#include <iostream>
#include <iomanip>

double series_x_mu_alternating(
  double x,
  double alpha,
  double mu,
  double delta,
  int maxiter = 200,
  double eps = 5e-16
)
{
  // Parameters
  double xmu = x - mu;
  double xmu2 = xmu * xmu;
  double ad = alpha * delta;
  double aod = alpha / delta;

  // Constant
  double C = delta * std::exp(ad) / constants::pi;

  // Series: compute the first two terms with Bessel K_1 and K_2 recursion
  double kn = std::cyl_bessel_k(1, ad);
  double knp1 = std::cyl_bessel_k(2, ad);

  double num = C * xmu * aod;
  double s = num * kn;

  double aux = -xmu2 * aod / 2.0;

  num *= aux;
  s += num * knp1 / 3.0;

  // Start recursion
  double sp = s;
  for (int k = 2; k < maxiter; k++)
  {
    // Bessel K recursion
    double knn = kn + 2.0 * k / ad * knp1;

    // New term
    num *= aux / k;
    s += num * knn / (2.0 * k + 1.0);

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      return 0.5 + s;
    }
    else {
      kn = knp1;
      knp1 = knn;
      sp = s;
    }
  }

  return -1.0;
}


double series_x_mu_pos(
  double x,
  double alpha,
  double mu,
  double delta,
  int maxiter = 200,
  double eps = 5e-16
)
{
  // Parameters
  double xmu = x - mu;
  double xmu2 = xmu * xmu;

  double omega = std::hypot(xmu, delta);
  double aw = alpha * omega;
  double aow = alpha / omega;

  // Constant
  double C = delta * std::exp(alpha * delta) / constants::pi;

  // Series: compute the first two terms with Bessel K_1 and K_2 recursion
  double kn = std::cyl_bessel_k(1, aw);
  double knp1 = std::cyl_bessel_k(2, aw);

  double num = C * xmu * aow;
  double s = num * kn;

  num *= xmu2 * aow / 3.0;
  s += num * knp1;

  // Start recursion
  double sp = s;
  for (int k = 2; k < maxiter; k++)
  {
    // Bessel recursion
    double knn = kn + 2.0 * k / aw * knp1;

    // New term
    num *= xmu2 * aow / (2.0 * k + 1.0);
    s += num * knn;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      return 0.5 + s;
    }
    else {
      kn = knp1;
      knp1 = knn;
      sp = s;
    }
  }

  return -1.0;
}


double series_x_mu_pos_robust(
  double x,
  double alpha,
  double mu,
  double delta,
  int maxiter = 1000,
  double eps = 5e-16
)
{
  // Parameters
  double xmu = x - mu;
  double xmu2 = xmu * xmu;

  double omega = std::hypot(xmu, delta);
  double aw = alpha * omega;
  double aow = alpha / omega;

  // Constants
  double z = xmu2 * aow;  
  double z2 = z * z;

  double C = delta * std::exp(alpha * delta) / constants::pi * xmu * aow;

  // Series: compute first two partial sums
  double y0 = 0.0;
  double y1 = std::cyl_bessel_k(1, aw);
  double y2 = y1 + std::cyl_bessel_k(2, aw) * z / 3.0;

  // Start recursion
  double sp = y2;
  for (int k = 0; k < maxiter; k++)
  {
    double aux0 = 2.0 * k;
    double aux1 = 3.0 + aux0;
    double aux2 = 5.0 + aux0;
    double aux3 = aux1 * aux2;

    double n0 = -z2 / aux3 * y0;
    
    double n11 = z* (-2.0 * (2.0 + k)) / (aux2 * aw);
    double n12 = z2 / aux3;
    double n1 = (n11 + n12) * y1;

    double n2 = (1.0 + (2.0 * (2.0 + k) * z) / (aux2 * aw)) * y2;

    double s = n0 + n1 + n2;

    if (std::fabs(1.0 - sp / s) < eps)
    {
      return 0.5 + C * s;
    }
    else {
      y0 = y1;
      y1 = y2;
      y2 = s;
      sp = s;
    }
  }

  return -1.0;
}


double asymptotic_x_mu_neg(
  double x,
  double alpha,
  double mu,
  double delta,
  int maxiter = 200,
  double eps = 5e-16
)
{
  // Parameters
  double xmu = x - mu;
  double xmu2 = xmu * xmu;

  double omega = std::hypot(xmu, delta);
  double aw = alpha * omega;
  double woa = omega / alpha;

  // Constant
  double C = delta * std::exp(alpha * delta) / constants::pi;

  // Series: compute the first two terms with K_0 and K_1
  double kn = std::cyl_bessel_k(0, aw);
  double knp1 = std::cyl_bessel_k(1, aw);

  double oxmu2 = 1.0 / xmu2;
  double num = -1.0 / xmu;
  double s = num * kn;
  double aux = -woa * oxmu2;

  num *= aux;
  s += num * knp1;

  double sp = s;
  for (int k = 2; k < maxiter; k++)
  {
    // Bessel recursion
    double knn = kn + 2.0 * (k - 1.0) / aw * knp1;

    // New term
    num *= aux * (2.0 * k - 1.0);
    s += num * knn;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps)
    {
      return C * s;
    }
    else {
      kn = knp1;
      knp1 = knn;
      sp = s;
    }
  }

  return -1.0;
}


double asymptotic_alpha(
  double x,
  double alpha,
  double mu,
  double delta,
  int maxiter = 200,
  double eps = 5e-16
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

  const double C = delta * std::exp(alpha * delta) * std::sqrt(z) * constants::oneosqrttwopi;

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
  const double kh = std::cyl_bessel_k(0.5, t0);
  const double khp2 = std::cyl_bessel_k(1.5, t0);

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
    double ck = (np1 * c1 * (t1 - 4.0 * b * n - t2) - (2 * n * n + n) * c0) / (t3 * np1 * (n + 2));

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


double nig_series_beta_zero(double x, double alpha, double mu, double delta)
{
  return asymptotic_alpha(x, alpha, mu, delta);
}