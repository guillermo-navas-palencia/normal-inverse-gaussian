#include <cmath>

#include "nig.hpp"


double nig_series_x_eq_mu(double alpha, double beta, double delta)
{
  return std::cyl_bessel_k(2, 0.5);
}