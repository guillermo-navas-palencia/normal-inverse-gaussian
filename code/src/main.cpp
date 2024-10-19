#include <iostream>
#include <iomanip>

#include <nig.hpp>


void test_case_x_eq_mu()
{
  double alpha = 0.005;
  double beta = 0.004;
  double delta = 50.001;

  double nig_cdf = nig_x_eq_mu(alpha, beta, delta);

  std::cout << std::setprecision(16) << nig_cdf << std::endl;
}


void test_case_beta_zero()
{
  double alpha = 1.0;
  double mu = 0.5;
  double delta = 3.0;
  double x = 1.1;

  double nig_cdf = nig_series_beta_zero(x, alpha, mu, delta);

  std::cout << std::setprecision(16) << nig_cdf << std::endl;
}


int main()
{
  // main_normal_distribution();
  // test_case_beta_zero();

  double alpha = 10.005;
  double beta = 0.0025;
  double mu = 0.5;
  double delta = 0.1;
  double x = 0.5;
  // double eps = 5e-16;

  // std::cout << truncation(alpha, beta, mu, delta, eps) << std::endl;
  // double tau = 1e-32;
  // std::cout << estimate_h(tau, 1e-10, 10) << std::endl;
  // double result = nig_integration(x, alpha, beta, mu, delta, 1e-15, 10);
  double result = nig_cdf(x, alpha, beta, mu, delta);
  std::cout << std::setprecision(16) << result << std::endl;

  return 0;
}