#include <cmath>

#include <nig.hpp>
#include <constants.hpp>


double normal_cdf(const double x)
{
  return 0.5 * std::erfc(-x / constants::sqrt2);
}


double normal_pdf(const double x)
{
  return std::exp(-x*x * 0.5) * constants::oneosqrttwopi;
}


void integrand_and_deriv(
  double t,
  const double beta,
  const double xmu,
  const double gamma2,
  const double delta2,
  double &f,
  double &fp
)
{
  double sqrtt = std::sqrt(t);
  double t2 = t * t;
  double t32 = t * sqrtt;
  double betat = beta * t;
  double z1aux = xmu - betat;
  double z2aux = xmu + betat;

  double z = z1aux / sqrtt;

  // exp integrand at t
  double A = z2aux / t32;
  double B = normal_pdf(z) / normal_cdf(z);

  f = 0.5 * (delta2 / t2 - gamma2 - 3.0 / t - A * B);

  // derivative exp integrand at t
  double d1 = -0.5 * (beta / t32 - 1.5 * z2aux / t / t32) * B;
  double d2 = -z * A * A * B;
  double d3 = -B * B * A * A;

  fp = -delta2 / t / t2 + 1.5 / t2 + d1 + d2 + d3;
}


double estimate_magnitude(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double gamma,
  const double C,
  const double tol = 1e-4,
  const int maxiter = 10
)
{
  // Estimate magnitude to accurately calculate N to achieve eps rel error

  // Constants
  const double alpha2 = alpha * alpha;
  const double delta2 = delta * delta;
  const double gamma2 = gamma * gamma;
  const double xmu = x - mu;
  const double xmu2 = xmu * xmu;

  double a, c;
  if (x - mu < 0.0)
  {
    a = alpha2;
    c = delta2 + xmu2;
  } else {
    a = gamma2;
    c = delta2;
  }

  double x0 = (-1.5 + std::sqrt(2.25 + a * c)) / a;

  // Newton's iteration
  double f, fp;
  for (int k = 0; k < maxiter; k++)
  {
    integrand_and_deriv(x0, beta, xmu, gamma2, delta2, f, fp);
    x0 -= f / fp;

    if (std::abs(f) < tol)
      break;
  }

  // Maximum integrand contribution
  double z0 = xmu / std::sqrt(x0) - beta * std::sqrt(x0);
  double g1 = -0.5 * (delta2 / x0 + gamma2 * x0);
  double g2 = -1.5 * std::log(x0);
  double g3 = std::log(normal_cdf(z0));
  double gx0 = g1 + g2 + g3 + delta * gamma;
  
  return C * std::exp(gx0);
}


int truncation(
  const double delta,
  const double gamma,
  const double eps = 5e-16
)
{
  // Solve: 2 e^{-gamma^2 / 2 N} / N^(3/2) / gamma^2 = eps / C
  // N = 3/gamma^2 W0(gamma^2/(3u)), u = (gamma^2 eps / 2 / C) ^(2/3)

  // Use upper bound of the W0(x): W0(x) < log(x)^(log(x) / (1 + log(x)))
  // To avoid overflow/underflow perform computation using logarithms

  // Constants
  const double gamma2 = gamma * gamma;
  const double loggamma2 = std::log(gamma2);

  // log C = log(delta) + delta * gamma - 1/2 * log(2 pi)
  const double logC = std::log(delta) + delta * gamma - 0.9189385332046727;

  // log y = log(gamma^2) - log(3) - 2/3 (log(gamma^2) + log(eps) - log(2) - logC)
  const double logy =  1./3 * loggamma2 - 2./3 * (std::log(eps) - logC) - 0.636514168294813;
  const double lambertwy = std::pow(logy, logy / (1.0 + logy));

  return (int) std::ceil(3.0 / gamma2 * lambertwy);
}


double estimate_h(
  const double tau,
  const double tol = 1e-10,
  const int maxiter = 10
)
{
  // Solve log(pi / tau) = pi / 2 - exp(x) - x - log(x)

  // Initial guess x0
  const double piotau = constants::pi / tau;
  double x = std::log(2.0 / constants::pi * std::log(piotau));

  // Newton's iteration
  double aux = constants::pihalf * std::exp(x);
  double fx = aux - x - std::log(x * piotau);

  for (int k = 0; k < maxiter; k++)
  {
    double fxp = aux - 1.0 - 1.0 / x;
    x -= fx / fxp;

    aux = constants::pihalf * std::exp(x);
    fx = aux - x - std::log(x * piotau);

    if (std::abs(fx) < tol)
      return x;
  }

  return x;
}


double lambertwm1(const double x)
{
  // D. A. Barry, L. Li, and D.-S. Jeng, “Comments on “Numerical evaluation
  // of the Lambert W function and application to generation of generalized
  // Gaussian noise with exponent 1/2”,” IEEE Trans. Signal Process.,
  // vol. 52, no. 5, pp. 1456–1458, May 2004.
  const double a = 0.3205;
  const double logmx = std::log(-x);

  return logmx - 2.0 / a * (1.0 - 1.0 / (1.0 + a * std::sqrt(-0.5 * (1.0 + logmx))));
}


double nig_integration(
  const double x,
  const double alpha,
  const double beta,
  const double mu,
  const double delta,
  const double eps = 1e-15,
  const int maxlevel = 10
)
{
  // Parameters
  const double gamma = std::sqrt(alpha * alpha - beta * beta);

  // Constants
  const double xmu = x - mu;
  const double C = delta * constants::oneosqrttwopi;

  // Estimate integrand magnitude to achieve eps relative error
  double magnitude = estimate_magnitude(x, alpha, beta, mu, delta, gamma, C);

  // Estimate the truncation point N
  int N = truncation(delta, gamma, eps * std::min(magnitude, 1.0));

  // Estimate initial step size h
  const double eps2 = eps * eps;
  double h = estimate_h(eps2);

  // Start integration
  const double alpha2 = N / 2.0;

  double estimate = 0.0;
  double f, y0, y1, weight, xl, rl, phil, fl, xr, rr, phir, fr;

  for (int level = 0; level <= maxlevel; level++)
  {
    // Estimate j using approximation W_{-1}(x)

    double japprox = std::log(-2.0 / constants::pi * lambertwm1(-eps2 / h / 2.0)) / h;
    int j = (int) std::ceil(japprox);

    if (level == 0.0)
    {
      y0 = alpha2;
      y1 = -alpha2 * constants::pihalf;
      weight = -h * y1;

      // Eval f(y0) and f(N - y0)
      xl = y0;
      rl = delta - gamma * xl;
      phil = normal_cdf((xmu - beta * xl) / std::sqrt(xl));
      fl = phil * std::exp(-rl*rl / xl / 2.0) * std::pow(xl, -1.5);

      xr = N - xl;
      rr = delta - gamma * xr;
      phir = normal_cdf((xmu - beta * xr) / std::sqrt(xr));
      fr = phir * std::exp(-rr*rr / xr / 2.0) * std::pow(xr, -1.5);

      estimate = fl * weight;
      h /= 2.0;
    }
    else
    {
      double suml = 0.0;
      double sumr = 0.0;
      double t, sinh_t, cosh_t, cosh_sinh_t, exp_sinh_t;

      for (int i = 1; i < j + 1; i += 2)
      {
        t = h * i;
        sinh_t = constants::pihalf * std::sinh(t);
        cosh_t = constants::pihalf * std::cosh(t);
        cosh_sinh_t = std::cosh(sinh_t);
        exp_sinh_t = std::exp(sinh_t);

        y0 = alpha2 / exp_sinh_t / cosh_sinh_t;
        y1 = -alpha2 * cosh_t / (cosh_sinh_t * cosh_sinh_t);
        weight = -h * y1;

        // Eval f(y0) and f(N - y0)
        xl = y0;
        rl = delta - gamma * xl;
        phil = normal_cdf((xmu - beta * xl) / std::sqrt(xl));
        fl = phil * std::exp(-rl*rl / xl / 2.0) * std::pow(xl, -1.5);

        xr = N - xl;
        rr = delta - gamma * xr;
        phir = normal_cdf((xmu - beta * xr) / std::sqrt(xr));
        fr = phir * std::exp(-rr*rr / xr / 2.0) * std::pow(xr, -1.5);

        suml += fl * weight;
        sumr += fr * weight;
      }

      f = suml + sumr + estimate / 2.0;

      if (std::abs(1.0 - estimate / f) < eps)
      {
        return C * f;
      } else {
        estimate = f;
        h /= 2.0;
      }
    }
  }

  double integral = C * estimate;

  // Check if order of magnitude of the computed integral is close to the
  // estimate. If the magnitudes differ significantly, an error during the
  // computation occurred. Then, return the estimated magnitude.
  if (std::abs(std::log10(integral) - std::log(magnitude)) > 5)
    return estimate;
  else
    return integral;
}