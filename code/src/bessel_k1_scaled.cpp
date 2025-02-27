/*  The modified Bessel function K1(x) for real x.
 *  
 *  Implementation based on SPECFUN CALCK1: https://www.netlib.org/specfun/k1
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


double specfun::bessel_k1_scaled(const double x)
{
  if (x < constants::epsilon)
    return 1.0 / x;

  if (x <= 1.0) {
    const double t = x * x;

    const double sp = (
      (((4.8127070456878442310E-1 * t + 9.9991373567429309922E+1) * t +
        7.1885382604084798576E+3) * t + 1.7733324035147015630E+5) * t +
        7.1938920065420586101E+5) * t -2.2149374878243304548E+6;

    const double sq = ((t + -2.8143915754538725829E+2) * t +
      3.7264298672067697862E+4) * t -2.2149374878243304548E+6;

    const double sf = (
      ((-2.2795590826955002390E-1 * t + -5.3103913335180275253E+1) * t +
        -4.5051623763436087023E+3) * t + -1.4758069205414222471E+5) * t
        -1.3531161492785421328E+6;

    const double sg = ((t - 3.0507151578787595807E+2) * t +
      4.3117653211351080007E+4) * t - 2.7062322985570842656E+6;

    return std::exp(x) * (t * std::log(x) * sf/sg + sp/sq) / x;

  } else {
    const double t = 1.0 / x;

    const double sp = (((((((((6.4257745859173138767E-2 * t +
      7.5584584631176030810E+0) * t + 1.3182609918569941308E+2) * t +
      8.1094256146537402173E+2) * t + 2.3123742209168871550E+3) * t + 
      3.4540675585544584407E+3) * t + 2.8590657697910288226E+3) * t +
      1.3319486433183221990E+3) * t + 3.4122953486801312910E+2) * t +
      4.4137176114230414036E+1) * t + 2.2196792496874548962E+0;

    const double sq = ((((((((t + 3.6001069306861518855E+1) * t +
      3.3031020088765390854E+2) * t + 1.2082692316002348638E+3) * t +
      2.1181000487171943810E+3) * t + 1.9448440788918006154E+3) * t +
      9.6929165726802648634E+2) * t + 2.5951223655579051357E+2) * t +
      3.4552228452758912848E+1) * t + 1.7710478032601086579E+0;

    return sp / sq / std::sqrt(x);
  }
}