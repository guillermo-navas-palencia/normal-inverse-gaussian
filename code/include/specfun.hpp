#ifndef specfun_hpp
#define specfun_hpp

namespace specfun 
{
  double bessel_k0_scaled(const double x);
  double bessel_k1_scaled(const double x);
  double erfc(const double x);
  double norm_cdf(const double x);
  double norm_pdf(const double x);
  double norm_cdf_nag(const double x);
  double norm_cdf_std(const double x);
}


#endif // specfun_hpp