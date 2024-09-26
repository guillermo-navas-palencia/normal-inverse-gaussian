#include <cmath>

#include <nig.hpp>
#include <constants.hpp>

#include <iostream>
#include <iomanip>

double series_pos(
  double alpha,
  double beta,
  double delta,
  int maxiter = 200,
  double eps = 5e-16
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constant
  const double ad = alpha * delta;
  const double doa = delta / alpha;
  const double aux = beta * beta * doa;

  const double C = delta * std::exp(delta * gamma) / constants::pi;

  // Series: compute the first two terms with Bessel K_0 and K_1 recursion
  double kn = std::cyl_bessel_k(0, ad);
  double knp1 = std::cyl_bessel_k(1, ad);

  double num = - C * beta;
  double s = num * kn;

  num *= aux / 3.0;
  s += num * knp1;

  // Start recursion
  double sp = s;
  for (int k = 2; k < maxiter; k++)
  {
    // Bessel recursion
    double knn = kn + 2 * (k - 1) / ad * knp1;

    // New term
    num *= aux / (2 * k + 1);
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

double series_pos_robust(
  double alpha,
  double beta,
  double delta,
  int maxiter = 2000,
  double eps = 5e-16
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constant
  const double ad = alpha * delta;
  const double doa = delta / alpha;
  const double z = beta * beta * doa;
  const double z2 = z * z;

  // Constants
  const double C = -beta * delta * std::exp(delta * gamma) / constants::pi;

  // Series: compute first two partial sums
  double y0 = 0.0;
  double y1 = std::cyl_bessel_k(0, ad);
  double y2 = y1 + std::cyl_bessel_k(1, ad) * z / 3.0;

  // Start recursion
  double sp = y2;
  for (int k = 0; k < maxiter; k++)
  {
    int aux0 = 2 * k;
    int aux1 = 3 + aux0;
    int aux2 = 5 + aux0;
    int aux3 = aux1 * aux2;

    double n0 = -z2 / aux3 * y0;

    double n11 = -z * 2 * k / (aux1 * ad);
    double n12 = (-6 + ad * z) * z / (aux3 * ad);
    double n1 = (n11 + n12) * y1;

    double n2 = (1 + (2 * (1 + k) * z) / (aux2 * ad)) * y2;

    double s = n0 + n1 + n2;

    if (std::fabs(1.0 - sp / s) < eps)
    {
      std::cout << k << std::endl;

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


double nig_series_x_eq_mu(double alpha, double beta, double delta)
{
  return series_pos(alpha, beta, delta);
}