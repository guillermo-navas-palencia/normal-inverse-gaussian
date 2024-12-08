#include <iostream>
#include <iomanip>

#include <chrono>

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

  double nig_cdf = nig_beta_eq_zero(x, alpha, mu, delta);

  std::cout << std::setprecision(16) << nig_cdf << std::endl;
}


int main()
{
  // main_normal_distribution();
  // test_case_beta_zero();

  // double alpha = 40.05;
  // double beta = 0;
  // double mu = 0.04;
  // double delta = 26.5;
  // double x = 1.1;

  // double x = -1.8816854798951448;
  // double beta = 0.0;
  // double alpha = 23.184287219731825;
  // double mu = -0.6779473810531575;
  // double delta = 46.86837947038307;

  // double x = -1200;
  // double alpha = 1;
  // double beta = -1./2;
  // double mu = 1;
  // double delta = 10;

  // double x = -1.9103180092440883;
  // double alpha = 22.41626678065764;
  // double beta = 20.0;
  // double mu = 2.9491781940248964;
  // double delta = 20.14817937115534;  
  // double eps = 5e-16;

  // double x = -20;
  // double alpha = 50;
  // double beta = 12;
  // double mu = 1;
  // double delta = 45;

  // double x = 4.219899000274793;
  // double alpha = 28.04845308531407;
  // double beta = 10.0;
  // double mu = -0.3756820526459652;
  // double delta = 40.4430694917478499;

  // std::cout << truncation(alpha, beta, mu, delta, eps) << std::endl;
  // double tau = 1e-32;
  // std::cout << estimate_h(tau, 1e-10, 10) << std::endl;
  // double result = nig_integration(x, alpha, beta, mu, delta, 1e-13, 14);
  // double result = nig_cdf(x, alpha, 0, mu, delta);

  // double x = -10./8;
  // double alpha = 9.0;
  // double beta = -5.0;
  // double mu = 1.0;
  // double delta = 24.0;

  double x = -0.4519899000274793;
  double alpha = 0.54845308531407;
  double beta = 0.3;
  double mu = 0.3756820526459652;
  double delta = 0.4430694917478499;

  double result = nig_integration(x, alpha, beta, mu, delta, 1e-13, 14);
  std::cout << std::setprecision(16) << result << std::endl;

  // double result;
  int N = 1000;
  auto start_time = std::chrono::high_resolution_clock::now();
  for(int count = 0; count < N; count++)
  {
    nig_integration(x, alpha, beta, mu, delta, 1e-13, 14);
  }

  // Record end time
  auto finish_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish_time - start_time;
  std::cout << "Elapsed time for gcc erfc " << elapsed.count() * 1000000 / N << " microseconds\n";

  return 0;
}