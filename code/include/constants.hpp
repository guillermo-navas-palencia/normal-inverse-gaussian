#ifndef constants_hpp
#define constants_hpp

#include <limits>

namespace constants {
  // machine precision
  const double epsilon       = std::numeric_limits<double>::epsilon();

  // limits
  const double inf           = std::numeric_limits<double>::infinity();
  const double maxdouble     = std::numeric_limits<double>::max();
  const double mindouble     = std::numeric_limits<double>::min();

  // common constants
  const double pi            = 3.1415926535897932385; // pi
  const double pihalf        = 1.5707963267948966192; // pi/2
  const double piquart       = 0.7853981633974483096; // pi/4
  const double pithreequart  = 2.3561944901923449288; // 3pi/4
  const double sqrtpi        = 1.7724538509055160272; // sqrt(pi)
  const double lnpi          = 1.1447298858494001741; // log(pi)
  const double twopi         = 6.2831853071795864769; // 2*pi
  const double oneopi        = 0.3183098861837906715; // 1/pi
  const double oneotwopi     = 0.1591549430918953358; // 1/(2*pi)
  const double oneosqrtpi    = 0.5641895835477562869; // 1/sqrt(pi)
  const double twoosqrtpi    = 1.1283791670955125739; // 2/sqrt(pi)
  const double sqrttwopi     = 2.5066282746310005024; // sqrt(2*pi)
  const double sqrttwoopi    = 0.7978845608028653559; // sqrt(2/pi)
  const double oneosqrttwopi = 0.3989422804014326779; // 1/sqrt(2*pi)
  const double lnsqrttwopi   = 0.9189385332046727418; // log(sqrt(2*pi))
  const double onethird      = 0.3333333333333333333; // 1/3
  const double twothird      = 0.6666666666666666667; // 2/3
  const double onesix        = 0.1666666666666666667; // 1/6
  const double twoexp14      = 1.1892071150027210667; // 2**(1/4)
  const double twoexp13      = 1.2599210498948731648; // 2**(1/3)
  const double sqrt2         = 1.4142135623730950488; // sqrt(2)
  const double sqrt3         = 1.7320508075688772935; // sqrt(3)
  const double log_2         = 0.6931471805599453094; // log(2)
  const double log_10        = 2.3025850929940456840; // log(10)    
}

#endif // constants_hpp