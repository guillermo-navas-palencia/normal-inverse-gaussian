/*  The modified Bessel function K0(x) for real x.
 *  
 *  Implementation based on SPECFUN CALCK0: https://www.netlib.org/specfun/k0
 *  
 *  Authors: W. J. Cody and Laura Stoltz
 *    Mathematics and Computer Science Division
 *    Argonne National Laboratory
 *    Argonne, IL 60439
 * 
 * 
 *  Guillermo Navas-Palencia <g.navas.palencia@gmail.com>
 *  Copyright (C) 2024
 */

#include <cmath>

#include <constants.hpp>
#include <specfun.hpp>


double specfun::bessel_k0_scaled(const double x)
{
  if (x < constants::epsilon)
    return 0.11593151565841245 - std::log(x);

  if (x <= 1.0) {
    const double t = x * x;

    const double sp = (
      (((5.8599221412826100000E-04 * t + 1.3166052564989571850E-01) * t +
        1.1999463724910714109E+01) * t + 4.6850901201934832188E+02) * t +
        5.9169059852270512312E+03) * t + 2.4708152720399552679E+03;

    const double sq = (
      (-2.4994418972832303646E+02 + t) * t + 2.1312714303849120380E+04);

    const double sf = (
      (-1.6414452837299064100E+00 * t -2.9601657892958843866E+02) * t +
       -1.7733784684952985886E+04) * t + -4.0320340761145482298E+05;

    const double sg = (
      (-2.5064972445877992730E+02 + t) * t +
      2.9865713163054025489E+04) * t + -1.6128136304458193998E+06;

    const double logx = std::log(x);
    return std::exp(x) * (sp / sq - t * sf * logx / sg - logx);

  } else {
    const double t = 1.0 / x;

    const double sp = ((((((((1.1394980557384778174E+02 * t +
      3.6832589957340267940E+03) * t + 3.1075408980684392399E+04) * t +
      1.0577068948034021957E+05) * t + 1.7398867902565686251E+05) * t + 
      1.5097646353289914539E+05) * t + 7.1557062783764037541E+04) * t +
      1.8321525870183537725E+04) * t + 2.3444738764199315021E+03) * t +
      1.1600249425076035558E+02;

    const double sq = (((((((((t + 2.0013443064949242491E+02) * t +
      4.4329628889746408858E+03) * t + 3.1474655750295278825E+04) * t +
      9.7418829762268075784E+04) * t + 1.5144644673520157801E+05) * t +
      1.2689839587977598727E+05) * t + 5.8824616785857027752E+04) * t +
      1.4847228371802360957E+04) * t + 1.8821890840982713696E+03) * t +
      9.2556599177304839811E+01;

    return sp / sq / std::sqrt(x);
  }
}