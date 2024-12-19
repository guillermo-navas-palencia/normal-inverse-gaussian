#ifndef nig_hpp
#define nig_hpp

#include <cstddef>


double nig_integration(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double eps,
  const size_t maxlevel
);

double nig_x_eq_mu(
  const double alpha,
  const double beta,
  const double delta
);

double nig_beta_eq_zero(
  const double x,
  const double alpha,
  const double mu,
  const double delta
);

double nig_general(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta
);

double nig_cdf(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta
);

#endif // nig_hpp