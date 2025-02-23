#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>

#include <nig.hpp>


int main()
{
  double x = 2.0;
  double alpha = 2.0;
  double beta = -0.4;
  double mu = 1.75;
  double delta = 2.0;

  double result1 = nig_integration(x, alpha, beta, mu, delta, 1e-13, 14);
  std::cout << std::setprecision(16) << result1 << std::endl;

  double result2 = nig_cdf(x, alpha, beta, mu, delta);
  std::cout << std::setprecision(16) << result2 << std::endl;

  int N = 10000;
  auto start_time = std::chrono::high_resolution_clock::now();
  for(int count = 0; count < N; count++)
    nig_cdf(x, alpha, beta, mu, delta);

  // Record end time
  auto finish_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish_time - start_time;
  std::cout << "Elapsed time " << elapsed.count() << " seconds\n";
  std::cout << "Elapsed time " << elapsed.count() * 1000000 / N << " microseconds\n";

  return 0;
}