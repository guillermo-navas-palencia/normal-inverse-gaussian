#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>

#include <cmath>
#include <unordered_map>

#include <nig.hpp>
#include <constants.hpp>
#include <specfun.hpp>

#include <nig.hpp>


void test_cash_or_nothing_beta_zero()
{
  double alpha = 40.0;
  double beta = 0.0;
  double delta = 25.0;
  double mu = 0.0;

  double St = 3000;
  double K = 5000;
  double r = 0.01;
  double q = 0.0;
  // double tau = 1.0 / 12.0;
  double tau = 2.0;

  double gamma = std::sqrt(alpha*alpha - beta*beta);
  double omega = delta * (std::sqrt(alpha*alpha - (beta + 1)*(beta + 1)) - gamma);
  double k = std::log(St / K) + (r - q + omega) * tau;

  double cdf = nig_general(k, alpha, -beta, mu, delta * tau);
  double call = St * std::exp(-q * tau) * cdf;

  std::cout << std::setprecision(16) << "CALL " << call << std::endl;
  std::cout << "k = " << k << std::endl;
}


double bessel_series_xmu_2(
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
  const double xmu2 = xmu * xmu;
  const double aw = alpha * omega;
  const double oaw = 2.0 / aw;
  const double z = -omega * beta / alpha / xmu;
  const double v = xmu2 * alpha / omega;

  const double C = delta * xmu * alpha / omega * constants::oneopi * std::exp(
    delta * gamma + xmu * beta - aw);

  // Bessel recursion: First three terms
  const double k0 = specfun::bessel_k0_scaled(aw);
  const double k1 = specfun::bessel_k1_scaled(aw);

  std::unordered_map<int, double> bessel_map;
  bessel_map.insert({0, k0});
  bessel_map.insert({1, k1});

  // Start recursion
  double t = C;
  double s = 0.0;
  double sp = s;

  for (size_t k = 0; k < maxiter; k++)
  {
    // Compute polynomial A(k). First iteration
    if (bessel_map.find(k+1) == bessel_map.end()) {
      double cached = bessel_map[k - 1] + k * oaw * bessel_map[k];
      bessel_map.insert({k + 1, cached});
    }

    double sA = bessel_map[k + 1];
    double u = 1.0;
    double r = 1.0;

    for (size_t j = 1; j <= 2 * k + 1; j++)
    {
      int m = k + 1 - j;
      int am = std::fabs(m);

      u *= z;
      r *= (2.0 * k + 2 - j) / j;
      sA += r * u * bessel_map[am];
    }

    // New term
    t /= (2 * k + 1);
    s += t * sA;

    // Check convergence
    if (std::fabs(1.0 - s / sp) < eps) {
      // std::cout << "iteration: " << k << std::endl;
      return 0.5 + s;
    } else {
      t *= v;
      sp = s;
    }
  }

  return -1.0;
}


double bessel_series_2(
  const double x,
  const double alpha,
  const double mu,
  const double delta,
  const size_t maxiter = 10000,
  const double eps = 5e-15
)
{
  // Parameters
  const double xmu = x - mu;
  const double xmu2 = xmu * xmu;
  const double omega = std::hypot(xmu, delta);

  // Constants
  const double aw = alpha * omega;
  const double aow = alpha / omega;
  const double da = delta * alpha;
  const double z = xmu2 * aow;
  const double oaw = 1.0 / aw;

  // Check if scaled version is required
  double C, k0, k1, t;
  const bool scaled = da > 705.342;

  // Constant C = delta * alpha / omega / pi * e^(delta * alpha)
  const double caux = delta * aow * constants::oneopi;

  if (scaled)
    C = std::fma(alpha, delta - omega, std::log(caux));
  else 
    C = caux * std::exp(da);

  // Series: compute the ratio of Bessel K_1 / K_2 for recursion
  if (scaled) {
    k0 = specfun::bessel_k0_scaled(aw);
    k1 = specfun::bessel_k1_scaled(aw);
    t = k1;
  } else {
    const double expaw = std::exp(-aw);
    k0 = specfun::bessel_k0_scaled(aw) * expaw;
    k1 = specfun::bessel_k1_scaled(aw) * expaw;
    t = C * k1;
  }

  double rp = k1 / k0;
  double s = t;

  // Start recursion
  double sp = s;
  for (size_t k = 1; k < maxiter; k++)
  {
    // Ratio Bessel recursion: r = 1 / rp + 2k / aw
    double r = std::fma(2 * k, oaw, 1.0 / rp);

    // New term
    t *= z / (2 * k + 1) * r;
    s += t;

    // Check convergence
    if (std::fabs(1.0 - sp / s) < eps){
      // std::cout << "bessel series beta=0 iterations: " << k << std::endl;
      if (scaled)
        return std::fma(xmu, std::exp(C + std::log(s)), 0.5);
      else
        return std::fma(xmu, s, 0.5);
    }
    else {
      rp = r;
      sp = s;
    }
  }

  return -1.0;
}


void analytical_european_option_call()
{
  double alpha = 40.0;
  double beta = 0.0;
  double delta = 25.0;
  double mu = 0.0;

  double St = 3000;
  double r = 0.01;
  double q = 0.0;
  double tau = 1 / 12.0;

  auto start_time = std::chrono::high_resolution_clock::now();

  double gamma = std::sqrt(alpha*alpha - beta*beta);
  double omega = delta * (std::sqrt(alpha*alpha - (beta + 1)*(beta + 1)) - gamma);

  for (int K = 2000; K < 5000; K++)
  {
      double k = std::log(St / K) + (r - q + omega) * tau;
      double cdf1 = bessel_series_2(k, alpha, mu, delta * tau, 1000, 1e-8);
      double cdf2 = bessel_series_xmu_2(k, alpha, -beta-1, mu, delta * tau, 1000, 1e-8);

      // double cdf1 = nig_general(k, alpha, -beta, mu, delta * tau);
      // double cdf2 = nig_general(k, alpha, -beta-1, mu, delta * tau);

      double call = St * std::exp(-q * tau) * cdf2 - K * std::exp(-r * tau) * cdf1;
      std::cout << std::setprecision(16) << call << std::endl;
  }
  
  auto finish_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish_time - start_time;
  std::cout << "Elapsed time " << elapsed.count() << " seconds\n";
}


int test()
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

int main()
{
  // test_cash_or_nothing_beta_zero();
  // analytical_european_option_call();

  return 0;
}