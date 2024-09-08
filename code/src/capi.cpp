#include "nig.hpp"
#include "nig.h"


double cpp_nig_cdf(double x, double alpha, double beta, double mu, double delta)
{
  return nig_cdf(x, alpha, beta, mu, delta);
}