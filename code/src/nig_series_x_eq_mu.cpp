#include <cmath>

#include <nig.hpp>
#include <constants.hpp>


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


double nig_series_x_eq_mu(double alpha, double beta, double delta)
{
  return series_pos(alpha, beta, delta);
}