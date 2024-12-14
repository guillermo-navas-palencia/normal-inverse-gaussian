/*  The complementary error function erfc(x) for real x.
 *  
 *  Implementation based on SPECFUN CALERF: https://www.netlib.org/specfun/erf
 *  
 *  Authors: W. J. Cody
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


double specfun::erfc(const double x)
{
    const double y = std::fabs(x);
    const double y2 = y * y;

    if (y <= 0.46875) {
        double t = (y > constants::epsilon) ? y2 : 0.0;

        const double sa = (((1.85777706184603153E-1 * t +
            3.16112374387056560E00) * t + 1.13864154151050156E02) * t +
            3.77485237685302021E02) * t + 3.20937758913846947E03;

        const double sb = (((t + 2.36012909523441209E01) * t +
            2.44024637934444173E02) * t + 1.28261652607737228E03) * t +
            2.84423683343917062E03;

        return 1.0 - x * sa / sb;

    } else if (y <= 4.0) {
        const double sc = (((((((2.15311535474403846E-8 * y +
            5.64188496988670089E-1) * y + 8.88314979438837594E00) * y +
            6.61191906371416295E01) * y + 2.98635138197400131E02) * y +
            8.81952221241769090E02) * y + 1.71204761263407058E03) * y +
            2.05107837782607147E03) * y + 1.23033935479799725E03;

        const double sd = (((((((y + 1.57449261107098347E01) * y +
            1.17693950891312499E02) * y + 5.37181101862009858E02) * y +
            1.62138957456669019E03) * y + 3.29079923573345963E03) * y +
            4.36261909014324716E03) * y + 3.43936767414372164E03) * y +
            1.23033935480374942E03;

        double result = std::exp(-y2) * sc / sd;
        return (x < 0.0) ? 2.0 - result : result;

    } else {
        const double t = 1.0 / y2;

        const double sp = ((((1.63153871373020978E-2 * t +
            3.05326634961232344E-1) * t + 3.60344899949804439E-1) * t +
            1.25781726111229246E-1) * t + 1.60837851487422766E-2) * t +
            6.58749161529837803E-4;

        const double sq = ((((t + 2.56852019228982242E00) * t +
            1.87295284992346047E00) * t + 5.27905102951428412E-1) * t +
            6.05183413124413191E-2) * t + 2.33520497626869185E-3;

        double result = (constants::oneosqrtpi - t * sp / sq) / y * std::exp(-y2);
        return (x < 0.0) ? 2.0 - result : result;
    }
}


double specfun::norm_cdf(const double x)
{
    return 0.5 * specfun::erfc(-x * constants::osqrt2);
}


double specfun::norm_cdf_std(const double x)
{
    return 0.5 * std::erfc(-x * constants::osqrt2);
}


double specfun::norm_cdf_nag(const double x)
{
// We only want to handle the case of x<=0.0
    const double p = std::exp(-0.5*x*x);
    const double z = -std::fabs(x);
    // Make exp(-z^2/2) *P[z]/Q[z] using Horner forms
    const double numerator = 0.5000000000000000 +
        z * (-0.7726229273322750 +
          z * (0.5907342938552783 +
            z * (-0.2869333616209955 +
              z * (0.0966067226123471 +
                z * (-0.02328333056701721 +
                  z * (0.004016012307190090 +
                    z * (-0.0004795982820969556 +
                      z * (0.00003624685524802775 + 
                        z * (-1.339331478644128e-6)))))))));
    
    const double denominator = 1.000000000000000 +
        z *(-2.343130415467415 +
          z * (2.551016170159611 +
            z * (-1.703679452304189 +
              z * (0.7752274191855893 +
                z * (-0.2520420647830403 +
                  z * (0.05955811540207907 + 
                    z * (-0.01015750738006275 +
                      z * (0.001205531820751269 +
                        z * (-0.0000908573922286600 +
                          z * (3.357206153485080e-6))))))))));
    
    const double h = p * (numerator / denominator);
    // Compute ND(x) taking care if x is +ve or -ve
    const double nd = (x <= 0.0 ? h : 1.0 - h);
    return(nd);
}